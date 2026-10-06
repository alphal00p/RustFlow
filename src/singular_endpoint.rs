//! Supplied-boundary limits in an explicitly admitted local endpoint chart.
//!
//! Endpoint values are terminal evidence, not initial data for a singular ODE.
//! The retained regular matching boundary supplies the reusable initial data.
use crate::kinematics::KinematicPath;
use crate::transport_cache::{
    BoundaryAccuracy, BoundaryIdentity, CachedBoundary, CachedPoint, EpsilonRange, RootGerm,
};
use crate::{ComplexFloat as C, Error, FlowOptions, Precision, Result, RunContext};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::domains::float::RoundingDirection::Up;
use symbolica::prelude::*;

/// Exact physical approach, with its endpoint at parameter zero. Rational maps
/// such as s=1/z also describe infinity without inventing a finite coordinate.
#[derive(Clone, Debug)]
pub struct EndpointChart {
    pub path: KinematicPath,
    pub matching_parameter: Atom,
    /// Discrete root signs at the regular matching point, never at a vanishing root.
    pub matching_germ: Option<RootGerm>,
    /// The same log(z)+2*pi*i*winding convention is used for powers and logs.
    pub winding: i32,
    /// Caller-declared terminal homotopy, in the flow's existing branch domain.
    pub homotopy: String,
}

impl EndpointChart {
    pub fn matching_point(&self) -> Result<CachedPoint> {
        self.path.validate()?;
        // An exact numerical constant is checked by the ordinary point validator.
        CachedPoint::Exact(BTreeMap::from([(
            self.path.parameter,
            self.matching_parameter.clone(),
        )]))
        .validate()?;
        if self.matching_parameter.together().cancel().is_zero() || self.homotopy.trim().is_empty()
        {
            return Err(Error::InvalidInput(
                "endpoint chart requires a nonzero matching parameter and a homotopy declaration"
                    .into(),
            ));
        }
        let rules = BTreeMap::from([(
            Atom::var(self.path.parameter),
            self.matching_parameter.clone(),
        )]);
        let point = CachedPoint::Exact(
            self.path
                .coordinates
                .iter()
                .map(|(&s, a)| (s, crate::family::substitute(a, &rules).together().cancel()))
                .collect(),
        );
        match &self.matching_germ {
            Some(germ) => point.with_root_germ(germ.clone()),
            None => {
                point.validate()?;
                Ok(point)
            }
        }
    }

    pub(crate) fn validate(&self, identity: &BoundaryIdentity) -> Result<()> {
        identity.validate_point(&self.matching_point()?)
    }

    pub(crate) fn key(&self, identity: &BoundaryIdentity) -> Result<String> {
        self.validate(identity)?;
        let coordinates = self
            .path
            .coordinates
            .iter()
            .map(|(&s, a)| (Atom::var(s).to_canonical_string(), a.to_canonical_string()))
            .collect::<BTreeMap<_, _>>();
        let data = serde_json::json!([
            "finite-epsilon-coefficient-endpoint-v1",
            identity.key(),
            Atom::var(self.path.parameter).to_canonical_string(),
            coordinates,
            self.matching_parameter.to_canonical_string(),
            self.matching_point()?.key()?,
            self.winding,
            self.homotopy,
        ]);
        Ok(
            blake3::hash(&serde_json::to_vec(&data).map_err(|e| Error::Cache(e.to_string()))?)
                .to_hex()
                .to_string(),
        )
    }
}

/// Explicit resource limits; these do not declare achieved accuracy.
#[derive(Clone, Debug)]
pub struct EndpointOptions {
    pub max_lift_dimension: usize,
    pub series_order: usize,
}
impl Default for EndpointOptions {
    fn default() -> Self {
        Self {
            max_lift_dimension: 256,
            series_order: 64,
        }
    }
}
impl EndpointOptions {
    pub(crate) fn validate(&self) -> Result<()> {
        if self.max_lift_dimension == 0
            || self.max_lift_dimension > 4096
            || self.series_order < 8
            || self.series_order > 10000
        {
            return Err(Error::InvalidInput(
                "endpoint limits require dimension 1..4096 and order 8..10000".into(),
            ));
        }
        Ok(())
    }
}

/// The regular route's admission remains separate. This callback admits the
/// final radial approach on the stated logarithm branch into the physical domain.
pub trait EndpointAdmission {
    fn admit(&self, matching: &CachedBoundary, chart: &EndpointChart) -> Result<bool>;
}
impl<F> EndpointAdmission for F
where
    F: Fn(&CachedBoundary, &EndpointChart) -> Result<bool>,
{
    fn admit(&self, matching: &CachedBoundary, chart: &EndpointChart) -> Result<bool> {
        self(matching, chart)
    }
}

/// A finite coefficientwise endpoint, with its regular anchor and input evidence.
/// It is never eligible as an ordinary cached initial condition.
#[derive(Clone, Debug)]
pub struct EndpointBoundary {
    pub identity: BoundaryIdentity,
    pub chart: EndpointChart,
    pub range: EpsilonRange,
    pub coefficients: Vec<Vec<C>>,
    pub accuracy: BoundaryAccuracy,
    pub matching_boundary: CachedBoundary,
    /// Explicit asserted physical relations, separate from numerical evidence.
    pub constraints: Option<EndpointConstraints>,
}
impl EndpointBoundary {
    pub fn validate(&self) -> Result<()> {
        self.chart.validate(&self.identity)?;
        if let Some(constraints) = &self.constraints {
            constraints.preflight(&self.identity, self.range, &RunContext::default())?;
        }
        self.matching_boundary.validate()?;
        if self.matching_boundary.identity.key() != self.identity.key()
            || self.matching_boundary.point.key()? != self.chart.matching_point()?.key()?
            || self.matching_boundary.range != self.range
        {
            return Err(Error::InvalidInput(
                "endpoint matching boundary does not match its chart, identity or epsilon range"
                    .into(),
            ));
        }
        // Reuse the ordinary evidence/dimension checker at the regular anchor.
        CachedBoundary {
            identity: self.identity.clone(),
            point: self.matching_boundary.point.clone(),
            kind: self.matching_boundary.kind,
            range: self.range,
            coefficients: self.coefficients.clone(),
            accuracy: self.accuracy.clone(),
        }
        .validate()?;
        if self.accuracy.verified_digits() > self.matching_boundary.accuracy.verified_digits()
            || self.accuracy.input_verified_digits()
                > self.matching_boundary.accuracy.input_verified_digits()
        {
            return Err(Error::InvalidInput(
                "endpoint evidence exceeds its regular input evidence".into(),
            ));
        }
        Ok(())
    }
    pub(crate) fn key(&self) -> Result<String> {
        Ok(format!(
            "{}:{}:{}:{}",
            self.chart.key(&self.identity)?,
            self.range.leading,
            self.range.last,
            self.constraints
                .as_ref()
                .map(EndpointConstraints::key)
                .transpose()?
                .unwrap_or_default(),
        ))
    }
    pub(crate) fn truncated(&self, range: EpsilonRange) -> Result<Self> {
        if self.constraints.is_some() && range != self.range {
            return Err(Error::Unsupported("constrained endpoint evidence requires the exact epsilon range; dropping higher constraints requires a separate projection proof".into()));
        }
        let count = EpsilonRange::new(range.leading, range.last)?;
        if count.leading != self.range.leading || count.last > self.range.last {
            return Err(Error::InvalidInput(
                "endpoint epsilon range is not covered".into(),
            ));
        }
        let n = (i64::from(range.last) - i64::from(range.leading) + 1) as usize;
        let mut result = self.clone();
        result.range = range;
        result.coefficients.truncate(n);
        result.accuracy = self.accuracy.reindexed(
            self.accuracy.comparison_errors()[..n].to_vec(),
            "compatible endpoint cache hit",
        )?;
        result.matching_boundary.range = range;
        result.matching_boundary.coefficients.truncate(n);
        result.matching_boundary.accuracy = self.matching_boundary.accuracy.reindexed(
            self.matching_boundary.accuracy.comparison_errors()[..n].to_vec(),
            "endpoint regular anchor",
        )?;
        result.validate()?;
        Ok(result)
    }
}

pub struct EndpointResult {
    pub boundary: EndpointBoundary,
    pub matching_transport: Option<crate::diffexp::EpsilonSolution>,
    /// Regular records retained by the matching transport. Anchor rebinding
    /// and cache deduplication are not counted as new transport checkpoints.
    pub inserted_regular_points: usize,
    pub cache_hit: bool,
    pub boundary_attempts: Vec<crate::physical_transport::BoundaryAttempt>,
}

/// A terminal request is explicitly separate from an ordinary physical point.
#[derive(Clone, Debug)]
pub struct EndpointRequest {
    pub chart: EndpointChart,
    pub range: EpsilonRange,
    pub options: EndpointOptions,
}

/// Exact asserted physical relations and bounded exact coefficient work.
/// The relations may correlate epsilon coefficients in flattened order
/// `epsilon_offset * physical_dimension + component`.
#[derive(Clone, Debug)]
pub struct EndpointConstraints {
    pub asymptotic: crate::asymptotic::ExactAsymptoticConstraints,
    pub limits: crate::frobenius::ExactFrobeniusLimits,
}
impl EndpointConstraints {
    pub(crate) fn preflight(
        &self,
        identity: &BoundaryIdentity,
        range: EpsilonRange,
        context: &RunContext,
    ) -> Result<()> {
        EpsilonRange::new(range.leading, range.last)?;
        let count = usize::try_from(i64::from(range.last) - i64::from(range.leading) + 1)
            .map_err(|_| Error::Limit("constrained endpoint epsilon count overflow".into()))?;
        let dimension = identity
            .dimension()
            .checked_mul(count)
            .ok_or_else(|| Error::Limit("constrained endpoint dimension overflow".into()))?;
        self.asymptotic.validate(dimension, &self.limits, context)
    }
    pub(crate) fn key(&self) -> Result<String> {
        self.asymptotic.key()
    }
}
impl EndpointRequest {
    /// Bound the epsilon hierarchy before pulling back or allocating its
    /// coefficient matrices. The later root lift enforces this same cap as
    /// additional monomial blocks are discovered.
    pub(crate) fn preflight(&self, identity: &BoundaryIdentity) -> Result<usize> {
        self.options.validate()?;
        EpsilonRange::new(self.range.leading, self.range.last)?;
        let count = usize::try_from(i64::from(self.range.last) - i64::from(self.range.leading) + 1)
            .map_err(|_| Error::Limit("endpoint epsilon count overflow".into()))?;
        let physical = identity
            .dimension()
            .checked_mul(count)
            .ok_or_else(|| Error::Limit("endpoint physical hierarchy dimension overflow".into()))?;
        if physical > self.options.max_lift_dimension {
            return Err(Error::Limit(format!(
                "endpoint physical hierarchy dimension {physical} exceeds lift limit {} before source pullback",
                self.options.max_lift_dimension
            )));
        }
        Ok(count)
    }
}

pub(crate) struct PreparedEndpoint {
    prepared: crate::algebraic::PreparedAlgebraicFrobenius,
    request: EndpointRequest,
    radius: Rational,
    guard_count: usize,
    constraints: Option<EndpointConstraints>,
}

impl PreparedEndpoint {
    pub(crate) fn new(
        mut source: crate::algebraic::AlgebraicSystem,
        identity: &BoundaryIdentity,
        request: &EndpointRequest,
        constraints: Option<&EndpointConstraints>,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        context.cancellation.check()?;
        let count = request.preflight(identity)?;
        request.chart.validate(identity)?;
        let z = request.chart.path.parameter;
        source
            .nonzero_conditions
            .extend(identity.prescribed_path_conditions(&request.chart.path)?);
        source.nonzero_conditions.extend(
            crate::physical_conditions::rational_denominator_conditions(
                &request
                    .chart
                    .path
                    .coordinates
                    .values()
                    .cloned()
                    .collect::<Vec<_>>(),
                &BTreeSet::from([z]),
            )?,
        );
        if let Some(constraints) = constraints {
            constraints.preflight(identity, request.range, context)?;
            if !source.roots.is_empty() {
                return Err(Error::Unsupported("constrained registered-root endpoints require an exact sheet-consistency proof".into()));
            }
        }
        let prepared = source.prepare_frobenius(request.options.max_lift_dimension, context)?;
        if constraints.is_none() {
            prepared.rational_preparation().admit_endpoint_support(
                identity.dimension().checked_mul(count).ok_or_else(|| {
                    Error::Limit("endpoint physical hierarchy dimension overflow".into())
                })?,
                context,
            )?;
        }
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        let mut guards = prepared.endpoint_guards(p)?;
        let rational = prepared.rational_preparation();
        for matrix in [
            &rational.normalized_system().matrix,
            rational.basis_transformation(),
        ] {
            guards.extend(crate::physical_conditions::rational_denominator_conditions(
                &matrix.iter().flatten().cloned().collect::<Vec<_>>(),
                &BTreeSet::from([z]),
            )?);
        }
        let guards =
            crate::physical_conditions::canonical_conditions(&guards, &BTreeSet::from([z]))?;
        let exact = crate::ode::source::ExactSpecialization::new(p, &Default::default())?;
        let mut radius = Rational::one();
        for guard in &guards {
            context.cancellation.check()?;
            let coefficients = exact.polynomial(guard, z)?;
            // Only the endpoint's zero factor is removed. Every other source,
            // normalization and chart hole must lie outside this open disk.
            let first = coefficients
                .iter()
                .position(|c| !c.is_zero())
                .ok_or_else(|| {
                    Error::InvalidInput("identically zero endpoint source guard".into())
                })?;
            let constant = &coefficients[first];
            let lower = constant.re.clone().abs().max(constant.im.clone().abs());
            let sum = coefficients[first + 1..]
                .iter()
                .fold(Rational::zero(), |sum, c| {
                    sum + c.re.clone().abs() + c.im.clone().abs()
                });
            if !sum.is_zero() {
                radius = radius.min(lower / (Rational::from(2) * sum));
            }
        }
        let matching = exact.polynomial(&request.chart.matching_parameter, z)?;
        if matching.len() != 1
            || matching[0].is_zero()
            || matching[0].re.clone().abs() + matching[0].im.clone().abs() >= radius
        {
            return Err(Error::InvalidInput(format!(
                "endpoint matching parameter must lie strictly inside the certified punctured disk |z| < {radius} (the admission check uses |Re z|+|Im z|)"
            )));
        }
        context.cancellation.check()?;
        Ok(Self {
            prepared,
            request: request.clone(),
            radius,
            guard_count: guards.len(),
            constraints: constraints.cloned(),
        })
    }

    pub(crate) fn evaluate(
        &self,
        mut matching: CachedBoundary,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<EndpointBoundary> {
        matching.validate()?;
        // A regular cache hit may have a Derived coordinate record at exactly
        // this requested point. Rebind only after the regular engine's equality
        // check; the matching values and all evidence remain unchanged.
        matching.point = self.request.chart.matching_point()?;
        let size = matching.identity.dimension();
        let count = matching.coefficients.len();
        let mut source = matching
            .coefficients
            .iter()
            .flatten()
            .cloned()
            .collect::<Vec<_>>();
        if self.constraints.is_some() {
            source.insert(
                0,
                Precision::decimal(options.digits + options.guard_digits)?.i(1),
            );
        }
        let seeds = self
            .request
            .chart
            .matching_germ
            .as_ref()
            .map(RootGerm::seeds)
            .unwrap_or_default();
        let mut digits = options.digits + options.guard_digits;
        let mut previous: Option<EndpointProfile> = None;
        let mut profiles = Vec::new();
        let mut last =
            "endpoint precision/order profiles did not establish requested accuracy".to_owned();
        for attempt in 0..=options.max_precision_attempts {
            context.cancellation.check()?;
            let p = Precision::decimal(digits)?;
            let order = self
                .request
                .options
                .series_order
                .checked_add(
                    attempt
                        .checked_mul(32)
                        .ok_or_else(|| Error::Limit("endpoint refinement order overflow".into()))?,
                )
                .filter(|order| *order <= 10000)
                .ok_or_else(|| Error::Limit("endpoint refinement order exceeds 10000".into()))?;
            let parameter = CachedPoint::Exact(BTreeMap::from([(
                self.request.chart.path.parameter,
                self.request.chart.matching_parameter.clone(),
            )]))
            .evaluate(p)?[&self.request.chart.path.parameter]
                .clone();
            let (map, reconstruction) = if let Some(constraints) = &self.constraints {
                match self.prepared.constrained_endpoint_map(
                    &parameter,
                    p,
                    order,
                    self.request.chart.winding,
                    &constraints.asymptotic,
                    &constraints.limits,
                    context,
                ) {
                    Ok(Some(map)) => (map.endpoint, Some(map.reconstruction)),
                    Ok(None) => {
                        last = "exact constrained endpoint prefix does not yet cover all requested or nonpositive coefficients".into();
                        digits = digits
                            .checked_add(20)
                            .ok_or_else(|| Error::Limit("endpoint precision overflow".into()))?;
                        continue;
                    }
                    Err(Error::InsufficientPrecision {
                        minimum_bits,
                        context: message,
                    }) => {
                        last = message;
                        digits = Precision::refinement_digits(digits, minimum_bits)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                }
            } else {
                let map = match self.prepared.endpoint_map(
                    &parameter,
                    &seeds,
                    p,
                    order,
                    self.request.chart.winding,
                    context,
                ) {
                    Ok(Some(map)) => map,
                    Ok(None) => {
                        last = "endpoint order does not reach every nonpositive exponent after the normalization's lowest Laurent shift".into();
                        digits = digits
                            .checked_add(20)
                            .ok_or_else(|| Error::Limit("endpoint precision overflow".into()))?;
                        continue;
                    }
                    Err(Error::InsufficientPrecision {
                        minimum_bits,
                        context: message,
                    }) => {
                        last = message;
                        digits = Precision::refinement_digits(digits, minimum_bits)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                };
                (map, None)
            };
            profiles.push((p.bits, order));
            let values = contract(&map.values, &source, p);
            let reconstruction_values = reconstruction
                .as_ref()
                .map(|map| contract(&map.values, &source, p));
            if let Some(prior) = &previous {
                let mut input = matching
                    .accuracy
                    .comparison_errors()
                    .iter()
                    .flatten()
                    .zip(matching.coefficients.iter().flatten())
                    .map(|(error, value)| {
                        let scale = unit_scale(value, p);
                        let floor = p
                            .tolerance(matching.accuracy.verified_digits())
                            .mul_round(&scale, p.bits, Up);
                        error.add_round(&floor, p.bits, Up)
                    })
                    .collect::<Vec<_>>();
                if self.constraints.is_some() {
                    input.insert(0, p.real(0));
                }
                let errors = contraction_errors(
                    &map,
                    &prior.map,
                    &values,
                    &prior.values,
                    &source,
                    &input,
                    p,
                    digits,
                );
                if let (
                    Some(reconstruction),
                    Some(reconstructed),
                    Some((prior_map, prior_values)),
                ) = (
                    &reconstruction,
                    &reconstruction_values,
                    &prior.reconstruction,
                ) {
                    let reconstruction_errors = contraction_errors(
                        reconstruction,
                        prior_map,
                        reconstructed,
                        prior_values,
                        &source,
                        &input,
                        p,
                        digits,
                    );
                    let numerical_errors = contraction_errors(
                        reconstruction,
                        prior_map,
                        reconstructed,
                        prior_values,
                        &source,
                        &vec![p.real(0); input.len()],
                        p,
                        digits,
                    );
                    let stable =
                        reconstructed
                            .iter()
                            .zip(&numerical_errors)
                            .all(|(value, error)| {
                                *error
                                    <= p.tolerance(options.digits).mul_round(
                                        &unit_scale(value, p),
                                        p.bits,
                                        Up,
                                    )
                            });
                    if !stable {
                        last = "constrained matching reconstruction has not converged across independent profiles".into();
                        previous = Some(EndpointProfile {
                            map: map.values,
                            values,
                            reconstruction: Some((
                                reconstruction.values.clone(),
                                reconstructed.clone(),
                            )),
                        });
                        digits = digits
                            .checked_add(20)
                            .ok_or_else(|| Error::Limit("endpoint precision overflow".into()))?;
                        continue;
                    }
                    for (i, (value, error)) in
                        reconstructed.iter().zip(&reconstruction_errors).enumerate()
                    {
                        let residual = p.norm(&p.sub(value, &source[i + 1]));
                        if residual > error.add_round(&input[i + 1], p.bits, Up) {
                            return Err(Error::Accuracy(format!(
                                "supplied regular component {i} is inconsistent with the asserted exact asymptotic relations and propagated input evidence"
                            )));
                        }
                    }
                }
                let values = values.chunks(size).map(<[C]>::to_vec).collect::<Vec<_>>();
                let errors = errors
                    .chunks(size)
                    .map(<[Float]>::to_vec)
                    .collect::<Vec<_>>();
                if crate::physical_transport::errors_meet(p, &values, &errors, options.digits) {
                    let verified = crate::physical_transport::strongest_evidence(
                        p,
                        &values,
                        &errors,
                        options.digits,
                        matching.accuracy.verified_digits().min(digits),
                    );
                    let provenance = format!(
                        "finite epsilon coefficient limit; independent (bits, order) profiles {profiles:?}; one native inverse and all physical unit directions per profile; propagated independent-component matching errors and input evidence floor; {} exact source/chart/normalization guards with punctured disk radius {}; consistency estimates, not a global error certificate",
                        self.guard_count, self.radius
                    );
                    let provenance = if let Some(constraints) = &self.constraints {
                        format!(
                            "{provenance}; exact affine nonpositive-support proof from explicitly asserted physical relations ({}); all regular components checked against constrained reconstruction and input evidence; conservative consistency checks, not a rigorous joint uncertainty certificate",
                            constraints.asymptotic.provenance
                        )
                    } else {
                        provenance
                    };
                    let accuracy = BoundaryAccuracy::endpoint(
                        &matching.accuracy,
                        verified,
                        p.bits,
                        errors,
                        &provenance,
                    )?;
                    let result = EndpointBoundary {
                        identity: matching.identity.clone(),
                        chart: self.request.chart.clone(),
                        range: self.request.range,
                        coefficients: values,
                        accuracy,
                        matching_boundary: matching,
                        constraints: self.constraints.clone(),
                    };
                    result.validate()?;
                    context.cancellation.check()?;
                    return Ok(result);
                }
                last = format!(
                    "endpoint profile comparison or propagated input uncertainty exceeds {} digits",
                    options.digits
                );
            }
            debug_assert_eq!(values.len(), size * count);
            previous = Some(EndpointProfile {
                map: map.values,
                values,
                reconstruction: reconstruction
                    .zip(reconstruction_values)
                    .map(|(map, values)| (map.values, values)),
            });
            digits = digits
                .checked_add(20)
                .ok_or_else(|| Error::Limit("endpoint precision overflow".into()))?;
        }
        Err(Error::Accuracy(last))
    }
}

struct EndpointProfile {
    map: Vec<Vec<C>>,
    values: Vec<C>,
    reconstruction: Option<(Vec<Vec<C>>, Vec<C>)>,
}
fn contract(map: &[Vec<C>], source: &[C], p: Precision) -> Vec<C> {
    map.iter()
        .map(|row| {
            row.iter()
                .zip(source)
                .fold(p.zero(), |sum, (a, b)| p.add(&sum, &p.mul(a, b)))
        })
        .collect()
}
#[allow(clippy::too_many_arguments)]
fn contraction_errors(
    map: &crate::algebraic::EndpointLinearMap,
    prior: &[Vec<C>],
    values: &[C],
    prior_values: &[C],
    source: &[C],
    input: &[Float],
    p: Precision,
    digits: u32,
) -> Vec<Float> {
    map.values
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut error = p.norm(&p.sub(&values[i], &prior_values[i]));
            for (j, transfer) in row.iter().enumerate() {
                let disagreement = p.norm(&p.sub(transfer, &prior[i][j]));
                let transfer_error =
                    disagreement.add_round(&map.arithmetic_errors[i][j], p.bits, Up);
                let transfer_upper = p.norm(transfer).add_round(&transfer_error, p.bits, Up);
                error =
                    error.add_round(&transfer_upper.mul_round(&input[j], p.bits, Up), p.bits, Up);
                error = error.add_round(
                    &transfer_error.mul_round(&p.norm(&source[j]), p.bits, Up),
                    p.bits,
                    Up,
                );
            }
            let magnitude = row.iter().zip(source).fold(p.real(0), |sum, (a, b)| {
                sum.add_round(&p.norm(a).mul_round(&p.norm(b), p.bits, Up), p.bits, Up)
            });
            error.add_round(
                &magnitude.mul_round(&p.tolerance(digits.saturating_sub(5)), p.bits, Up),
                p.bits,
                Up,
            )
        })
        .collect()
}

fn unit_scale(value: &C, p: Precision) -> Float {
    let norm = p.norm(value);
    if norm > p.real(1) { norm } else { p.real(1) }
}
