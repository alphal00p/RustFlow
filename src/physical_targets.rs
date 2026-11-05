//! Exact target reductions applied to cached Laurent master coefficients.
use crate::transport_cache::{CachedBoundary, CachedPoint, EpsilonRange};
use crate::{ComplexFloat, Error, Precision, PreparedPhysicalFamily, Result};
use std::{collections::BTreeMap, sync::Arc};
use symbolica::domains::float::{ComplexBall, FloatField, RoundingDirection};
use symbolica::poly::series::{Series, SeriesDepth};
use symbolica::prelude::*;

/// A target value projected from an independently verified master boundary.
/// Errors include inherited master uncertainty and precision refinement of the
/// projection. No new epsilon fit or independent master calculation is implied.
#[derive(Clone, Debug)]
pub struct ProjectedLaurentExpansion {
    pub coefficients: BTreeMap<i32, ComplexFloat>,
    pub absolute_errors: BTreeMap<i32, Float>,
    pub verified_digits: u32,
    pub working_bits: u32,
}

type ExactSeries = Series<symbolica::domains::atom::AtomField>;

#[derive(Clone, Copy)]
struct ExactTargetFraction {
    numerator: usize,
    denominator: usize,
}

/// Symbolica owns both rational normalization and the shared evaluation graph.
/// Only polynomials enter the graph: specialize every denominator/guard before
/// any inversion, including at complex points where the formal i satisfies i²=-1.
struct ExactTargetBinding {
    parameters: Vec<Symbol>,
    evaluator: ExpressionEvaluator<Complex<Rational>>,
    outputs: usize,
    coefficients: BTreeMap<Atom, ExactTargetFraction>,
    guards: Vec<ExactTargetFraction>,
}

impl ExactTargetBinding {
    fn new(rows: &[Vec<SeriesWeight>], guards: &[Atom]) -> Result<Option<Self>> {
        // A zero output also makes the all-zero target table an ordinary graph.
        let mut polynomials = vec![Atom::new()];
        let mut indices = BTreeMap::from([(Atom::new(), 0)]);
        let mut parameters = std::collections::BTreeSet::new();
        let mut fractions = BTreeMap::new();
        for coefficient in rows
            .iter()
            .flatten()
            .flat_map(|w| w.series.terms().map(|(_, a)| a))
            .chain(guards)
        {
            if fractions.contains_key(coefficient) {
                continue;
            }
            let rational: RationalPolynomial<IntegerRing, u16> =
                crate::family::encode_complex(coefficient)
                    .try_to_rational_polynomial(&Q, &Z, None)
                    .map_err(|e| {
                        Error::Unsupported(format!("rational target coefficient required: {e}"))
                    })?;
            for variable in rational.get_variables().iter() {
                let PolyVariable::Symbol(symbol) = variable else {
                    // Nonpolynomial atoms retain the legacy specialization path.
                    return Ok(None);
                };
                parameters.insert(*symbol);
            }
            let mut intern = |polynomial: MultivariatePolynomial<IntegerRing, u16>| {
                let expression = polynomial.to_expression();
                *indices.entry(expression.clone()).or_insert_with(|| {
                    let index = polynomials.len();
                    polynomials.push(expression);
                    index
                })
            };
            fractions.insert(
                coefficient.clone(),
                ExactTargetFraction {
                    numerator: intern(rational.numerator),
                    denominator: intern(rational.denominator),
                },
            );
        }
        let parameters = parameters.into_iter().collect::<Vec<_>>();
        let atoms = parameters
            .iter()
            .copied()
            .map(Atom::var)
            .collect::<Vec<_>>();
        let evaluator = Atom::evaluator_multiple(&polynomials, &atoms)
            .horner_iterations(0)
            .cpe_iterations(Some(1))
            .build()
            .map_err(|e| Error::Unsupported(format!("exact target binding evaluator: {e}")))?;
        Ok(Some(Self {
            parameters,
            evaluator,
            outputs: polynomials.len(),
            guards: guards.iter().map(|a| fractions[a]).collect(),
            coefficients: fractions,
        }))
    }

    /// None requests legacy preparation, never a zero coefficient or a rounded
    /// coordinate. Normalize coordinates once; all remaining arithmetic is Q(i).
    fn parameters_at_point(
        &self,
        point: &BTreeMap<Symbol, Atom>,
    ) -> Option<Vec<Complex<Rational>>> {
        let mut parameters = Vec::with_capacity(self.parameters.len());
        for symbol in &self.parameters {
            if *symbol == crate::family::imaginary_parameter() {
                parameters.push(Complex::new(Rational::from(0), Rational::from(1)));
                continue;
            }
            let Some(value) = point.get(symbol) else {
                return None;
            };
            let value = value.together().cancel();
            let AtomView::Num(n) = value.as_view() else {
                return None;
            };
            let Coefficient::Complex(value) = n.get_coeff_view().to_owned() else {
                return None;
            };
            parameters.push(value);
        }
        Some(parameters)
    }

    fn at_point(&self, point: &BTreeMap<Symbol, Atom>) -> Result<Option<BTreeMap<Atom, Atom>>> {
        let Some(parameters) = self.parameters_at_point(point) else {
            return Ok(None);
        };
        let ring = FloatField::from_rep(Complex::new(Rational::from(1), Rational::from(0)));
        let mut values = vec![ring.zero(); self.outputs];
        // Each point owns its scratch stack; the retained graph remains reusable.
        self.evaluator
            .clone()
            .try_evaluate_in_ring(&parameters, &mut values, &ring)
            .map_err(|e| Error::Unsupported(format!("exact target binding evaluation: {e}")))?;
        if self
            .guards
            .iter()
            .any(|g| ring.is_zero(&values[g.numerator]) || ring.is_zero(&values[g.denominator]))
            || self
                .coefficients
                .values()
                .any(|c| ring.is_zero(&values[c.denominator]))
        {
            return Ok(None);
        }
        let mut inverses = BTreeMap::new();
        let mut coefficients = BTreeMap::new();
        for (coefficient, fraction) in &self.coefficients {
            let inverse = inverses
                .entry(fraction.denominator)
                .or_insert_with(|| ring.try_inv(&values[fraction.denominator]).unwrap());
            let value = ring.mul(&values[fraction.numerator], &*inverse);
            coefficients.insert(
                coefficient.clone(),
                crate::family::encode_complex(&Atom::num(value)),
            );
        }
        Ok(Some(coefficients))
    }
}

/// Only exact parameters are retained at a point. Both numerical evaluations
/// map the shared native polynomial graph afresh at the owner's precision.
struct NumericPointBinding {
    graph: Arc<ExactTargetBinding>,
    parameters: Vec<Complex<Rational>>,
}

struct NumericCoefficient {
    low: ComplexFloat,
    high: ComplexFloat,
    // Refined coefficient roundoff, including exact-coordinate conversion.
    error: Float,
}

type NumericCoefficients = BTreeMap<Atom, NumericCoefficient>;

impl NumericPointBinding {
    fn evaluate(&self, p: Precision) -> Result<Option<BTreeMap<Atom, ComplexBall>>> {
        let convert =
            |c: &Complex<Rational>| ComplexBall::from_rational_ball(c, &Rational::from(0), p.bits);
        let parameters = self.parameters.iter().map(&convert).collect::<Vec<_>>();
        let mut evaluator = self
            .graph
            .evaluator
            .clone()
            .map_coeff_with_prec(&convert, p.bits);
        let zero = convert(&Complex::new(Rational::from(0), Rational::from(0)));
        let mut values = vec![zero; self.graph.outputs];
        evaluator
            .try_evaluate(&parameters, &mut values)
            .map_err(|e| Error::Numerical(format!("numerical target binding evaluation: {e}")))?;
        let mut coefficients = BTreeMap::new();
        for (coefficient, fraction) in &self.graph.coefficients {
            let denominator = &values[fraction.denominator];
            if denominator.contains_zero() {
                // Exact guards passed, but this precision cannot resolve the
                // denominator. Reuse legacy exact specialization, not its
                // rounded inverse or a silently vanishing weight.
                return Ok(None);
            }
            let value = values[fraction.numerator].clone() / denominator;
            if !value.re.center.is_finite()
                || !value.im.center.is_finite()
                || !value.re.radius.is_finite()
                || !value.im.radius.is_finite()
            {
                return Ok(None);
            }
            coefficients.insert(coefficient.clone(), value);
        }
        Ok(Some(coefficients))
    }

    fn coefficients(
        &self,
        p: Precision,
        refined: Precision,
    ) -> Result<Option<NumericCoefficients>> {
        let Some(low) = self.evaluate(p)? else {
            return Ok(None);
        };
        let Some(high) = self.evaluate(refined)? else {
            return Ok(None);
        };
        Ok(Some(
            low.into_iter()
                .map(|(coefficient, low)| {
                    let high = &high[&coefficient];
                    let error = high.re.radius.add_round(
                        &high.im.radius,
                        refined.bits,
                        RoundingDirection::Up,
                    );
                    (
                        coefficient,
                        NumericCoefficient {
                            low: ComplexFloat::new(low.re.center, low.im.center),
                            high: ComplexFloat::new(high.re.center.clone(), high.im.center.clone()),
                            error,
                        },
                    )
                })
                .collect(),
        ))
    }
}

fn expansion(coefficient: &Atom, epsilon: Symbol, depth: SeriesDepth) -> Result<ExactSeries> {
    crate::family::scalar_symbols(coefficient.as_view(), &mut Default::default())?;
    let coefficient = crate::family::encode_complex(coefficient);
    let _: RationalPolynomial<IntegerRing, u16> = coefficient
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::Unsupported(format!("rational target coefficient required: {e}")))?;
    coefficient
        .series(epsilon, 0, depth)
        .map_err(|e| Error::Unsupported(format!("target coefficient series: {e}")))
}

fn valuation(coefficient: &Atom, epsilon: Symbol) -> Result<Option<i32>> {
    if coefficient.is_zero() {
        return Ok(None);
    }
    let series = expansion(coefficient, epsilon, SeriesDepth::relative(1))?;
    if series.is_zero() {
        return Err(Error::Numerical(
            "nonzero target weight has no leading series coefficient".into(),
        ));
    }
    Ok(Some(series_leading(&series)?))
}

fn series_leading(series: &ExactSeries) -> Result<i32> {
    series
        .get_trailing_exponent()
        .to_string()
        .parse()
        .map_err(|_| Error::Limit("target weight epsilon valuation exceeds i32".into()))
}

/// Target weights expanded once over symbolic kinematics, in the owner's cached
/// master basis (including epsilon shearing). This owned table can be stored
/// beside its family without a self-reference. No master values are retained.
///
/// Bind a new point with `at_point` before transporting and projecting masters.
/// Away from specialization singularities this evaluates a retained native
/// exact polynomial graph; it does not re-expand specialized rational weights.
pub struct PreparedKinematicTargetProjection {
    identity_key: String,
    targets: Vec<crate::reduction::LinearCombination>,
    range: EpsilonRange,
    required: EpsilonRange,
    rows: Arc<Vec<Vec<SeriesWeight>>>,
    binding: Option<Arc<ExactTargetBinding>>,
    numeric_guards: Option<ExactTargetBinding>,
}

impl PreparedKinematicTargetProjection {
    /// Generic demand, valid when the retained specialization guards hold.
    /// At exceptional points use the range returned by `at_point` instead:
    /// denominator valuations can increase the required master depth.
    pub fn required_master_range(&self) -> EpsilonRange {
        self.required
    }

    /// Validate the original family identity and target map, then bind retained
    /// coefficients, preserving exact input and source accuracy. The returned
    /// point plan borrows only `family`, so this table can be moved or dropped.
    /// Zeros can lower demand. A vanishing denominator leading coefficient or
    /// removable series-coefficient singularity uses legacy exact preparation;
    /// generic valuations are never assumed at such a point. Rational and
    /// rational-complex coordinates use the same exact evaluator. Unsupported
    /// atoms or coordinates use legacy preparation without rounding.
    pub fn at_point<'a>(
        &self,
        family: &'a PreparedPhysicalFamily,
        point: &BTreeMap<Symbol, Atom>,
    ) -> Result<PreparedTargetProjection<'a>> {
        if family.flow().identity().key() != self.identity_key.as_str()
            || family.target_reductions() != self.targets.as_slice()
        {
            return Err(Error::InvalidInput(
                "kinematic target projection has a different physical basis identity or target map"
                    .into(),
            ));
        }
        let cached_point = CachedPoint::Exact(point.clone());
        family.flow().identity().validate_point(&cached_point)?;
        let Some(binding) = &self.binding else {
            return family.prepare_target_projection(point, self.range);
        };
        let Some(coefficients) = binding.at_point(point)? else {
            return family.prepare_target_projection(point, self.range);
        };
        let mut last = self.required.leading;
        let mut rows = Vec::with_capacity(self.rows.len());
        for row in self.rows.iter() {
            let mut bound = Vec::with_capacity(row.len());
            for weight in row {
                let series = weight
                    .series
                    .map_coeff(|coefficient| coefficients[coefficient].clone());
                if series.is_zero() {
                    // Only the retained orders vanish. Earlier master poles
                    // require legacy re-expansion, handled by the point plan.
                    continue;
                }
                let leading = series_leading(&series)?;
                last = last.max(required_last(leading, self.range)?);
                bound.push(SeriesWeight {
                    column: weight.column,
                    leading,
                    series,
                });
            }
            rows.push(bound);
        }
        Ok(PreparedTargetProjection {
            family,
            point: cached_point,
            range: self.range,
            required: EpsilonRange::new(self.required.leading, last)?,
            weights: None,
            rows: Arc::new(rows),
            numeric: None,
        })
    }

    /// Opt in to numerical evaluation of the retained symbolic coefficients.
    /// The exact point is preserved. Original denominator/epsilon guards are
    /// checked exactly without specializing coefficient numerators. Generic
    /// master demand is conservative, including weights that vanish here.
    ///
    /// Projection evaluates the shared Symbolica graph independently at the
    /// boundary's working precision and its refined precision. Native ball
    /// radii account for coefficient/coordinate roundoff in addition to the
    /// usual refinement estimate and inherited master uncertainty. No source
    /// credit is added. Exceptional/unsupported specializations, unresolved
    /// numerical denominators and earlier master poles use the exact route.
    /// Like `at_point`, the result borrows only `family`.
    pub fn at_point_numeric<'a>(
        &self,
        family: &'a PreparedPhysicalFamily,
        point: &BTreeMap<Symbol, Atom>,
    ) -> Result<PreparedTargetProjection<'a>> {
        if family.flow().identity().key() != self.identity_key.as_str()
            || family.target_reductions() != self.targets.as_slice()
        {
            return Err(Error::InvalidInput(
                "kinematic target projection has a different physical basis identity or target map"
                    .into(),
            ));
        }
        let cached_point = CachedPoint::Exact(point.clone());
        family.flow().identity().validate_point(&cached_point)?;
        let (Some(graph), Some(guards)) = (&self.binding, &self.numeric_guards) else {
            return family.prepare_target_projection(point, self.range);
        };
        let Some(parameters) = graph.parameters_at_point(point) else {
            return family.prepare_target_projection(point, self.range);
        };
        if guards.at_point(point)?.is_none() {
            return family.prepare_target_projection(point, self.range);
        }
        Ok(PreparedTargetProjection {
            family,
            point: cached_point,
            range: self.range,
            required: self.required,
            weights: None,
            rows: Arc::clone(&self.rows),
            numeric: Some(NumericPointBinding {
                graph: Arc::clone(graph),
                parameters,
            }),
        })
    }
}

/// Target weights bound to one exact physical point and target epsilon range.
///
/// Borrow the retained family while transporting the requested master range,
/// then project the resulting boundary with this same plan. No numerical values
/// or source accuracy are retained or upgraded by preparation. A boundary at a
/// different point or with a different physical basis identity is rejected.
///
/// Exact plans substitute coordinates, cancel and determine epsilon valuation
/// at preparation. Opt-in numerical plans retain symbolic coefficients and
/// conservative generic demand, evaluating only during projection. Absolute
/// weight series are retained through the ordinary master pole bound's depth. A boundary starting below that
/// bound is supported by expanding full point-bound weights more deeply. Plans
/// made by `prepare_target_projection` or `prepare_target_combinations` retain
/// those weights (including exact target combinations); plans made by
/// `PreparedKinematicTargetProjection::at_point` rebind them on that exceptional
/// path through the legacy API.
pub struct PreparedTargetProjection<'a> {
    family: &'a PreparedPhysicalFamily,
    point: CachedPoint,
    range: EpsilonRange,
    required: EpsilonRange,
    // Legacy plans keep full point-bound rational weights for deeper poles.
    // Kinematic plans bind only retained coefficients and rebind on that path.
    weights: Option<Vec<Vec<BoundWeight>>>,
    rows: Arc<Vec<Vec<SeriesWeight>>>,
    numeric: Option<NumericPointBinding>,
}

struct BoundWeight {
    column: usize,
    coefficient: Atom,
    leading: i32,
}

struct SeriesWeight {
    column: usize,
    leading: i32,
    series: ExactSeries,
}

impl PreparedTargetProjection<'_> {
    /// Master orders required by these weights (conservative for numerical plans).
    pub fn required_master_range(&self) -> EpsilonRange {
        self.required
    }

    /// Project a boundary at the plan's point using its retained weight series.
    /// Kinematic plans rebind for earlier-than-ordinary master poles; numerical
    /// plans also fall back when generic demand or denominators are inconclusive.
    /// The same accuracy, missing-order and pole guards as `project_targets`
    /// on the retained family apply, including support for earlier master poles.
    pub fn project_targets(
        &self,
        boundary: &CachedBoundary,
        digits: u32,
    ) -> Result<Vec<ProjectedLaurentExpansion>> {
        let master_leading = self
            .family
            .validate_target_projection(boundary, self.range, digits)?;
        if !self.point.same_coordinates(&boundary.point)? {
            return Err(Error::InvalidInput(
                "target projection plan has different physical point coordinates".into(),
            ));
        }
        if boundary.range.leading < master_leading && self.weights.is_none() {
            return self.family.project_targets(boundary, self.range, digits);
        }
        if self.numeric.is_some() && boundary.range.last < self.required.last {
            // A specialized zero may lower demand. Only exact preparation may
            // establish that; a numerical zero cannot invent missing orders.
            return self.family.project_targets(boundary, self.range, digits);
        }
        let (p, refined) = projection_precision(boundary)?;
        let numeric = if let Some(binding) = &self.numeric {
            let Some(coefficients) = binding.coefficients(p, refined)? else {
                return self.family.project_targets(boundary, self.range, digits);
            };
            Some(coefficients)
        } else {
            None
        };
        let mut result = Vec::with_capacity(self.rows.len());
        for (index, rows) in self.rows.iter().enumerate() {
            if let Some(weights) = &self.weights {
                validate_master_depth(
                    weights[index].iter().map(|w| w.leading),
                    self.range,
                    boundary,
                )?;
            } else {
                validate_master_depth(rows.iter().map(|w| w.leading), self.range, boundary)?;
            }
            // The old API uses the boundary's actual leading order. Do not
            // truncate earlier coefficients merely because the ordinary bound
            // used to prepare the default series is less negative.
            let deeper;
            let rows = if let Some(weights) = &self.weights
                && boundary.range.leading < master_leading
            {
                deeper = expand_row(
                    &weights[index],
                    self.family.family().epsilon,
                    self.range,
                    boundary.range.leading,
                )?;
                &deeper
            } else {
                rows
            };
            result.push(project_row(
                rows,
                self.family.family().epsilon,
                boundary,
                self.range,
                digits,
                p,
                refined,
                numeric.as_ref(),
            )?);
        }
        Ok(result)
    }
}

impl PreparedPhysicalFamily {
    /// Retain exact Laurent target weights over symbolic kinematics once for
    /// this target range. The generic demand can be conservative; always use
    /// the point plan's demand after specialization. Full exact preparation is
    /// repeated only at exceptional denominator specializations or for masters
    /// with poles earlier than the ordinary bound.
    pub fn prepare_kinematic_target_projection(
        &self,
        range: EpsilonRange,
    ) -> Result<PreparedKinematicTargetProjection> {
        EpsilonRange::new(range.leading, range.last)?;
        let leading = self.ordinary_master_leading()?;
        let weights = self.bind_target_weights(&BTreeMap::new())?;
        let required = master_range(&weights, leading, range)?;
        let mut rows = weights
            .iter()
            .map(|row| expand_row(row, self.family().epsilon, range, leading))
            .collect::<Result<Vec<_>>>()?;
        let mut guards = Vec::new();
        // Even a weight beyond the retained orders can acquire a pole after
        // specialization, e.g. eps^3/(s+eps^4) at s=0.
        for weight in weights.iter().flatten() {
            let denominator = target_denominator(&weight.coefficient)?;
            let denominator = expansion(
                &denominator,
                self.family().epsilon,
                SeriesDepth::relative(1),
            )?;
            guards.push(
                denominator
                    .coefficient(denominator.get_trailing_exponent())
                    .ok_or_else(|| {
                        Error::Numerical("missing target denominator leading coefficient".into())
                    })?,
            );
        }
        for row in &mut rows {
            for weight in row {
                weight.series = weight.series.map_coeff(|a| a.together().cancel());
                // Native series expressions may have removable coefficient
                // poles beyond those of the original rational denominator.
                for (_, coefficient) in weight.series.terms() {
                    guards.push(target_denominator(coefficient)?);
                }
            }
        }
        guards.retain(|a| !matches!(a.as_view(), AtomView::Num(_)));
        guards.sort_by_cached_key(AtomCore::to_canonical_string);
        guards.dedup();
        let binding = ExactTargetBinding::new(&rows, &guards)?.map(Arc::new);
        // A separate graph contains only original leading-denominator guards
        // and unique series denominators, never the coefficient numerators.
        let numeric_guards = ExactTargetBinding::new(&[], &guards)?;
        Ok(PreparedKinematicTargetProjection {
            identity_key: self.flow().identity().key().to_owned(),
            targets: self.target_reductions().to_vec(),
            range,
            required,
            rows: Arc::new(rows),
            binding,
            numeric_guards,
        })
    }

    pub(crate) fn ordinary_master_leading(&self) -> Result<i32> {
        let family = self.family();
        let mut symbols = std::collections::BTreeSet::new();
        for coefficient in family.external_gram.iter().flatten().chain(
            family
                .propagators
                .iter()
                .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products)),
        ) {
            crate::family::scalar_symbols(coefficient.as_view(), &mut symbols)?;
        }
        if symbols.contains(&Atom::var(family.epsilon)) {
            return Err(Error::Unsupported("automatic master pole bounds require epsilon-independent propagators and kinematics".into()));
        }
        i32::try_from(family.loops.len())
            .ok()
            .and_then(|n| n.checked_mul(-2))
            .ok_or_else(|| Error::Limit("master Laurent bound overflow".into()))
    }

    /// Bind exact target reductions once at this point for both master demand
    /// and subsequent numerical projection. Retain one plan per new point and
    /// target range; the family can be reused across requests. Preparation uses
    /// a relative leading-term series for valuation and an absolute series for
    /// projection; neither is repeated when projecting from the ordinary bound.
    pub fn prepare_target_projection(
        &self,
        point: &BTreeMap<Symbol, Atom>,
        range: EpsilonRange,
    ) -> Result<PreparedTargetProjection<'_>> {
        EpsilonRange::new(range.leading, range.last)?;
        let cached_point = CachedPoint::Exact(point.clone());
        self.flow().identity().validate_point(&cached_point)?;
        let leading = self.ordinary_master_leading()?;
        let weights = self.bind_target_weights(point)?;
        let required = master_range(&weights, leading, range)?;
        let rows = weights
            .iter()
            .map(|row| expand_row(row, self.family().epsilon, range, leading))
            .collect::<Result<_>>()?;
        Ok(PreparedTargetProjection {
            family: self,
            point: cached_point,
            range,
            required,
            weights: Some(weights),
            rows: Arc::new(rows),
            numeric: None,
        })
    }

    /// Prepare exact linear combinations of the retained targets at one point.
    /// `weights[h][t]` multiplies target `t` in `target_reductions()` order;
    /// the matrix must have at least one row and exactly one column per target.
    /// Weights must be exact rational/Gaussian rational functions of epsilon
    /// and the supplied kinematics. Numerical factors must be lifted exactly
    /// by the caller; this does not grant them additional source accuracy.
    ///
    /// Bind reductions and weights once, then sum and cancel by ordered cached
    /// master column before valuation, Laurent expansion or uncertainty norms.
    /// Thus targets sharing a master retain that common dependency; different
    /// master coefficients still use the ordinary error bounds without any
    /// assumed covariance.
    /// The returned plan has one output per matrix row and uses the existing
    /// point, basis, depth, source-floor and arithmetic-refinement checks.
    /// `required_master_range()` reports the ordinary master pole bound and
    /// the depth needed after cancellation; it does not assume omitted poles
    /// vanish or certify a stronger physical pole bound for these combinations.
    /// For a nonzero combined row its natural output bound is the ordinary
    /// master leading order plus its smallest surviving coefficient valuation.
    /// Request those earlier output orders explicitly to check stronger poles.
    pub fn prepare_target_combinations(
        &self,
        point: &BTreeMap<Symbol, Atom>,
        weights: &[Vec<Atom>],
        range: EpsilonRange,
    ) -> Result<PreparedTargetProjection<'_>> {
        EpsilonRange::new(range.leading, range.last)?;
        let targets = self.target_reductions().len();
        if weights.is_empty() || targets == 0 || weights.iter().any(|row| row.len() != targets) {
            return Err(Error::InvalidInput(
                "target combination weights require a nonempty matrix with one column per retained target".into(),
            ));
        }
        let cached_point = CachedPoint::Exact(point.clone());
        self.flow().identity().validate_point(&cached_point)?;
        let leading = self.ordinary_master_leading()?;
        // Do not determine individual target valuations: poles and whole
        // master dependencies can disappear only after the exact sum below.
        let reductions = self.bind_target_coefficients(point)?;
        let substitutions = point
            .iter()
            .map(|(&s, a)| (Atom::var(s), a.clone()))
            .collect();
        let variables = point
            .keys()
            .copied()
            .chain(std::iter::once(self.family().epsilon))
            .collect();
        let mut combined = Vec::with_capacity(weights.len());
        for row in weights {
            // Reuse the owner's exact rational-domain checks, including
            // denominators of weights multiplying an empty or cancelled target.
            let conditions =
                crate::physical_conditions::rational_denominator_conditions(row, &variables)?;
            crate::physical_conditions::validate_conditions_at(
                &conditions,
                self.family().epsilon,
                point,
            )?;
            let mut columns = BTreeMap::<usize, Atom>::new();
            for (weight, target) in row.iter().zip(&reductions) {
                let weight = crate::family::substitute(weight, &substitutions)
                    .together()
                    .cancel();
                for (column, reduction) in target {
                    *columns.entry(*column).or_default() += &weight * reduction;
                }
            }
            combined.push(bind_row(
                columns
                    .into_iter()
                    .map(|(column, coefficient)| (column, coefficient.together().cancel())),
                self.family().epsilon,
            )?);
        }
        let required = master_range(&combined, leading, range)?;
        let rows = combined
            .iter()
            .map(|row| expand_row(row, self.family().epsilon, range, leading))
            .collect::<Result<_>>()?;
        Ok(PreparedTargetProjection {
            family: self,
            point: cached_point,
            range,
            required,
            weights: Some(combined),
            rows: Arc::new(rows),
            numeric: None,
        })
    }

    /// Determine how far masters must be known before reducing targets. A pole
    /// in a reduction coefficient requires additional positive epsilon orders;
    /// a cache truncated at the target's last order is not sufficient in general.
    pub fn required_master_range(
        &self,
        point: &BTreeMap<Symbol, Atom>,
        targets: EpsilonRange,
    ) -> Result<EpsilonRange> {
        EpsilonRange::new(targets.leading, targets.last)?;
        self.flow()
            .identity()
            .validate_point(&CachedPoint::Exact(point.clone()))?;
        let leading = self.ordinary_master_leading()?;
        master_range(&self.bind_target_weights(point)?, leading, targets)
    }

    fn bind_target_weights(&self, point: &BTreeMap<Symbol, Atom>) -> Result<Vec<Vec<BoundWeight>>> {
        self.bind_target_coefficients(point)?
            .into_iter()
            .map(|row| bind_row(row, self.family().epsilon))
            .collect()
    }

    fn bind_target_coefficients(
        &self,
        point: &BTreeMap<Symbol, Atom>,
    ) -> Result<Vec<Vec<(usize, Atom)>>> {
        let substitutions = point
            .iter()
            .map(|(&s, a)| (Atom::var(s), a.clone()))
            .collect();
        let columns = self
            .basis()
            .iter()
            .enumerate()
            .map(|(column, integral)| (integral, column))
            .collect::<BTreeMap<_, _>>();
        self.target_reductions()
            .iter()
            .map(|target| {
                let mut row = Vec::with_capacity(target.len());
                for (integral, weight) in target {
                    let column = *columns.get(integral).ok_or_else(|| {
                        Error::IncompleteReduction(
                            "target reduction leaves the physical master basis".into(),
                        )
                    })?;
                    let weight = self.cache_target_weight(column, weight);
                    let coefficient = crate::family::substitute(&weight, &substitutions)
                        .together()
                        .cancel();
                    row.push((column, coefficient));
                }
                Ok(row)
            })
            .collect()
    }

    fn validate_target_projection(
        &self,
        boundary: &CachedBoundary,
        range: EpsilonRange,
        digits: u32,
    ) -> Result<i32> {
        EpsilonRange::new(range.leading, range.last)?;
        boundary.validate()?;
        if boundary.identity.key() != self.flow().identity().key() {
            return Err(Error::InvalidInput(
                "target projection has a different physical basis identity".into(),
            ));
        }
        if digits == 0 || digits > boundary.accuracy.verified_digits() {
            return Err(Error::Accuracy(
                "target projection requires sufficient inherited accuracy".into(),
            ));
        }
        let master_leading = self.ordinary_master_leading()?;
        if boundary.range.leading > master_leading {
            return Err(Error::InvalidInput(format!(
                "target projection requires master coefficients from epsilon order {master_leading}; omitted pole coefficients cannot be assumed zero"
            )));
        }
        Ok(master_leading)
    }

    /// Apply the retained exact target reduction at this boundary's coordinates.
    /// Insufficient master orders and excessive propagated uncertainty are typed
    /// errors; missing positive orders are never silently replaced by zero.
    pub fn project_targets(
        &self,
        boundary: &CachedBoundary,
        range: EpsilonRange,
        digits: u32,
    ) -> Result<Vec<ProjectedLaurentExpansion>> {
        self.validate_target_projection(boundary, range, digits)?;
        let coordinates = boundary.point.restart_coordinates()?;
        let (p, refined) = projection_precision(boundary)?;
        let weights = self.bind_target_weights(&coordinates)?;
        let mut result = Vec::with_capacity(weights.len());
        for row in &weights {
            validate_master_depth(row.iter().map(|w| w.leading), range, boundary)?;
            let rows = expand_row(row, self.family().epsilon, range, boundary.range.leading)?;
            result.push(project_row(
                &rows,
                self.family().epsilon,
                boundary,
                range,
                digits,
                p,
                refined,
                None,
            )?);
        }
        Ok(result)
    }
}

fn bind_row(
    coefficients: impl IntoIterator<Item = (usize, Atom)>,
    epsilon: Symbol,
) -> Result<Vec<BoundWeight>> {
    let mut row = Vec::new();
    for (column, coefficient) in coefficients {
        if let Some(leading) = valuation(&coefficient, epsilon)? {
            row.push(BoundWeight {
                column,
                coefficient,
                leading,
            });
        }
    }
    Ok(row)
}

fn target_denominator(coefficient: &Atom) -> Result<Atom> {
    let rational: RationalPolynomial<IntegerRing, u16> = crate::family::encode_complex(coefficient)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::Unsupported(format!("rational target coefficient required: {e}")))?;
    Ok(rational.denominator.to_expression())
}

fn master_range(
    weights: &[Vec<BoundWeight>],
    leading: i32,
    range: EpsilonRange,
) -> Result<EpsilonRange> {
    let mut last = leading;
    for weight in weights.iter().flatten() {
        last = last.max(required_last(weight.leading, range)?);
    }
    EpsilonRange::new(leading, last)
}

fn required_last(leading: i32, range: EpsilonRange) -> Result<i32> {
    range
        .last
        .checked_sub(leading)
        .ok_or_else(|| Error::Limit("target Laurent order overflow".into()))
}

fn validate_master_depth(
    weights: impl IntoIterator<Item = i32>,
    range: EpsilonRange,
    boundary: &CachedBoundary,
) -> Result<()> {
    for weight in weights {
        let required_last = required_last(weight, range)?;
        if required_last > boundary.range.last {
            return Err(Error::InvalidInput(format!(
                "target reduction requires master epsilon order {required_last}, but boundary ends at {}",
                boundary.range.last
            )));
        }
    }
    Ok(())
}

fn expand_row(
    weights: &[BoundWeight],
    epsilon: Symbol,
    range: EpsilonRange,
    master_leading: i32,
) -> Result<Vec<SeriesWeight>> {
    let mut rows = Vec::with_capacity(weights.len());
    for weight in weights {
        let last_weight = range
            .last
            .checked_sub(master_leading)
            .ok_or_else(|| Error::Limit("target weight order overflow".into()))?;
        if weight.leading > last_weight {
            continue;
        }
        let depth = i64::from(last_weight) + 1;
        rows.push(SeriesWeight {
            column: weight.column,
            leading: weight.leading,
            series: expansion(&weight.coefficient, epsilon, SeriesDepth::absolute(depth))?,
        });
    }
    Ok(rows)
}

fn projection_precision(boundary: &CachedBoundary) -> Result<(Precision, Precision)> {
    let p = Precision {
        bits: boundary.accuracy.working_bits(),
    };
    let refined = Precision {
        bits: p
            .bits
            .checked_add(64)
            .ok_or_else(|| Error::Limit("target projection precision overflow".into()))?,
    };
    Ok((p, refined))
}

fn project_row(
    rows: &[SeriesWeight],
    epsilon: Symbol,
    boundary: &CachedBoundary,
    range: EpsilonRange,
    digits: u32,
    p: Precision,
    refined: Precision,
    numeric: Option<&NumericCoefficients>,
) -> Result<ProjectedLaurentExpansion> {
    let count = (i64::from(range.last) - i64::from(range.leading) + 1) as usize;
    let mut values = vec![p.zero(); count];
    let mut high_values = vec![refined.zero(); count];
    let mut input_errors = vec![refined.zero(); count];
    for SeriesWeight {
        column,
        leading,
        series,
    } in rows
    {
        // Shift both Laurent series to ordinary series, then delegate
        // the products to Symbolica. Keep the absolute truncation length:
        // trailing zeros must not turn unknown higher terms into zeros.
        let offset = i64::from(boundary.range.leading) + i64::from(*leading);
        let length = i64::from(range.last) - offset + 1;
        if length <= 0 {
            continue;
        }
        let length = usize::try_from(length)
            .map_err(|_| Error::Limit("target product length overflow".into()))?;
        if length > boundary.coefficients.len() {
            return Err(Error::InvalidInput(
                "target product needs additional master coefficients".into(),
            ));
        }
        let mut weights_low = Vec::with_capacity(length);
        let mut weights_high = Vec::with_capacity(length);
        let mut weights_norm = Vec::with_capacity(length);
        let mut coefficient_errors = Vec::with_capacity(length);
        let mut master_norms = Vec::with_capacity(length);
        let mut masters = Vec::with_capacity(length);
        let mut errors = Vec::with_capacity(length);
        for index in 0..length {
            let power = i64::from(*leading) + index as i64;
            let weight = series.coefficient(Rational::from(power)).ok_or_else(|| {
                Error::Accuracy(
                    "target weight coefficient lies in the unknown series remainder".into(),
                )
            })?;
            let (low, high, coefficient_error) = if let Some(coefficients) = numeric {
                if weight.is_zero() {
                    (p.zero(), refined.zero(), refined.real(0))
                } else {
                    let coefficient = &coefficients[&weight];
                    (
                        coefficient.low.clone(),
                        coefficient.high.clone(),
                        coefficient.error.clone(),
                    )
                }
            } else {
                (
                    p.eval(&weight, &Default::default())?,
                    refined.eval(&weight, &Default::default())?,
                    refined.real(0),
                )
            };
            let mut norm = refined.norm(&high);
            if numeric.is_some() {
                norm = norm.add_round(&coefficient_error, refined.bits, RoundingDirection::Up);
                coefficient_errors.push(ComplexFloat::new(coefficient_error, refined.real(0)));
                master_norms.push(ComplexFloat::new(
                    refined.norm(&boundary.coefficients[index][*column]),
                    refined.real(0),
                ));
            }
            weights_norm.push(ComplexFloat::new(norm, refined.real(0)));
            weights_low.push(low);
            weights_high.push(high);
            let master = &boundary.coefficients[index][*column];
            let scale = maximum(refined.norm(master), refined.real(1));
            let floor = refined.tolerance(boundary.accuracy.verified_digits()) * scale;
            let error = maximum(
                boundary.accuracy.comparison_errors()[index][*column].clone(),
                floor,
            );
            errors.push(ComplexFloat::new(error, refined.real(0)));
            masters.push(master.clone());
        }
        let multiply = |p, left: &[ComplexFloat], right: &[ComplexFloat]| {
            let a = crate::fixed_series::fixed_series(p, epsilon, left)?;
            let b = crate::fixed_series::fixed_series(p, epsilon, right)?;
            crate::fixed_series::coefficients(&(&a * &b), length)
        };
        let low = multiply(p, &weights_low, &masters)?;
        let high = multiply(refined, &weights_high, &masters)?;
        let mut propagated = multiply(refined, &weights_norm, &errors)?;
        if numeric.is_some() {
            // |w|*delta_m + delta_w*(|m|+delta_m). The first
            // product already uses |w_high|+delta_w, so only
            // delta_w*|m| remains. Keep the same native convolution.
            let roundoff = multiply(refined, &coefficient_errors, &master_norms)?;
            for (error, roundoff) in propagated.iter_mut().zip(roundoff) {
                *error = refined.add(error, &roundoff);
            }
        }
        for (index, power) in (range.leading..=range.last).enumerate() {
            let source = i64::from(power) - offset;
            if source < 0 {
                continue;
            }
            let source = usize::try_from(source)
                .map_err(|_| Error::Limit("target product index overflow".into()))?;
            values[index] = p.add(&values[index], &low[source]);
            high_values[index] = refined.add(&high_values[index], &high[source]);
            input_errors[index] = refined.add(&input_errors[index], &propagated[source]);
        }
    }
    let mut coefficients = BTreeMap::new();
    let mut absolute_errors = BTreeMap::new();
    for (index, power) in (range.leading..=range.last).enumerate() {
        let value = &values[index];
        let high = &high_values[index];
        if !p.finite(value) || !refined.finite(high) {
            return Err(Error::Numerical("nonfinite projected target".into()));
        }
        let error = input_errors[index].re.add_round(
            &refined.norm(&refined.sub(high, value)),
            p.bits,
            symbolica::domains::float::RoundingDirection::Nearest,
        );
        let scale = maximum(p.norm(value), p.real(1));
        if !error.is_finite() || error > p.tolerance(digits) * scale {
            return Err(Error::Accuracy("target reduction amplifies cached master uncertainty beyond the requested accuracy".into()));
        }
        coefficients.insert(power, value.clone());
        absolute_errors.insert(power, error);
    }
    Ok(ProjectedLaurentExpansion {
        coefficients,
        absolute_errors,
        verified_digits: digits,
        working_bits: p.bits,
    })
}

fn maximum(a: Float, b: Float) -> Float {
    if a > b { a } else { b }
}
