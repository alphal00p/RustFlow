//! Exact target reductions applied to cached Laurent master coefficients.
use crate::transport_cache::{CachedBoundary, CachedPoint, EpsilonRange};
use crate::{ComplexFloat, Error, Precision, PreparedPhysicalFamily, Result};
use std::collections::BTreeMap;
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
/// Away from specialization singularities this only substitutes exact Laurent
/// coefficients; it does not expand the full specialized rational weights.
pub struct PreparedKinematicTargetProjection {
    identity_key: String,
    targets: Vec<crate::reduction::LinearCombination>,
    range: EpsilonRange,
    required: EpsilonRange,
    rows: Vec<Vec<SeriesWeight>>,
    specialization_guards: Vec<Atom>,
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
    /// generic valuations are never assumed at such a point.
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
        let mut substitutions = point
            .iter()
            .map(|(&s, a)| (Atom::var(s), a.clone()))
            .collect::<BTreeMap<_, _>>();
        // Test specialization in Q(i), not with i as an independent symbol.
        substitutions.insert(
            Atom::var(crate::family::imaginary_parameter()),
            Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
        );
        for guard in &self.specialization_guards {
            if crate::family::substitute(guard, &substitutions)
                .together()
                .cancel()
                .is_zero()
            {
                return family.prepare_target_projection(point, self.range);
            }
        }
        let mut last = self.required.leading;
        let mut rows = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            let mut bound = Vec::with_capacity(row.len());
            for weight in row {
                let series = weight.series.map_coeff(|coefficient| {
                    let coefficient = crate::family::substitute(coefficient, &substitutions)
                        .together()
                        .cancel();
                    crate::family::encode_complex(&coefficient)
                });
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
            rows,
        })
    }
}

/// Exact target weights bound to one physical point and target epsilon range.
///
/// Borrow the retained family while transporting the requested master range,
/// then project the resulting boundary with this same plan. No numerical values
/// or source accuracy are retained or upgraded by preparation. A boundary at a
/// different point or with a different physical basis identity is rejected.
///
/// Coordinate substitution, cancellation and epsilon valuation happen only at
/// preparation. Absolute weight series are also retained through the depth
/// needed for the ordinary master pole bound. A boundary starting below that
/// bound is supported by expanding full point-bound weights more deeply. Plans
/// made by `prepare_target_projection` retain those weights; plans made by
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
    rows: Vec<Vec<SeriesWeight>>,
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
    /// Master orders required by these exact point-bound target weights.
    pub fn required_master_range(&self) -> EpsilonRange {
        self.required
    }

    /// Project a boundary at the plan's point using its retained weight series.
    /// A kinematic plan rebinds only for earlier-than-ordinary master poles.
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
        let (p, refined) = projection_precision(boundary)?;
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
        Ok(PreparedKinematicTargetProjection {
            identity_key: self.flow().identity().key().to_owned(),
            targets: self.target_reductions().to_vec(),
            range,
            required,
            rows,
            specialization_guards: guards,
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
            rows,
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
                    if let Some(leading) = valuation(&coefficient, self.family().epsilon)? {
                        row.push(BoundWeight {
                            column,
                            coefficient,
                            leading,
                        });
                    }
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
            )?);
        }
        Ok(result)
    }
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
        let mut masters = Vec::with_capacity(length);
        let mut errors = Vec::with_capacity(length);
        for index in 0..length {
            let power = i64::from(*leading) + index as i64;
            let weight = series.coefficient(Rational::from(power)).ok_or_else(|| {
                Error::Accuracy(
                    "target weight coefficient lies in the unknown series remainder".into(),
                )
            })?;
            let low = p.eval(&weight, &Default::default())?;
            let high = refined.eval(&weight, &Default::default())?;
            weights_norm.push(ComplexFloat::new(refined.norm(&high), refined.real(0)));
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
        let propagated = multiply(refined, &weights_norm, &errors)?;
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
