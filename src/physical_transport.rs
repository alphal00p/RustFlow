//! Cache-first transport between regular physical kinematic points.
//!
//! Every request selects a compatible boundary before forming its path. Accepted
//! intermediate physical points are independently checked and retained. The
//! current orchestrator supports regular straight paths; threshold prescriptions
//! and partial singular boundaries are handled by separate solver interfaces.
use crate::algebraic::{
    AlgebraicKinematicSystem, AlgebraicSystem, CanonicalAlgebraicSystem, CompiledAlgebraicSystem,
};
use crate::diffexp::{EpsilonSolution, EpsilonSystem, transport_epsilon};
use crate::kinematics::{KinematicPath, KinematicSystem};
use crate::transport_cache::{
    BoundaryAccuracy, BoundaryIdentity, BoundaryQuery, CachedBoundary, CachedPoint, EpsilonRange,
    PointKind, RootGerm, RustFlowCache, TransportCost,
};
use crate::{Error, FlowOptions, Precision, Prescription, Result, RunContext};
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub struct RustFlow<S = KinematicSystem> {
    system: S,
    identity: BoundaryIdentity,
}

/// Outcome of an actually attempted compatible cached source.
#[derive(Clone, Debug)]
pub enum BoundaryAttemptOutcome {
    Accepted,
    AccuracyRejected { message: String },
}

#[derive(Clone, Debug)]
pub struct BoundaryAttempt {
    pub starting_point: CachedPoint,
    /// The arbitrary-precision cost returned by the caller's policy.
    pub cost: Float,
    pub source_verified_digits: u32,
    pub outcome: BoundaryAttemptOutcome,
}

pub struct PhysicalResult {
    pub boundary: CachedBoundary,
    pub starting_point: CachedPoint,
    /// None for a compatible cache hit at precisely the requested coordinates.
    pub transport: Option<EpsilonSolution>,
    pub inserted_points: usize,
    /// Compatible sources tried in selection order, including an exact hit.
    pub boundary_attempts: Vec<BoundaryAttempt>,
}

impl RustFlow {
    pub fn new(
        system: KinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
    ) -> Result<Self> {
        Self::with_conditions(
            system,
            basis,
            normalization,
            prescription,
            branch_domain,
            &[],
        )
    }

    /// Preserve all exact assumptions used in deriving the physical system.
    pub fn with_conditions(
        system: KinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        let identity = BoundaryIdentity::with_conditions(
            &system,
            basis,
            normalization,
            prescription,
            branch_domain,
            conditions,
        )?;
        Ok(Self { system, identity })
    }

    /// Select a verified cached boundary and transport to this exact point.
    /// The cost policy must validate physical path/sheet compatibility. Nearby
    /// entries are starting values, never replacements for the requested point.
    /// The entire cache update is committed only after successful validation.
    pub fn evaluate_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
    ) -> Result<PhysicalResult> {
        evaluate_physical(
            Connection::Rational(&self.system),
            &self.identity,
            cache,
            CachedPoint::Exact(destination.clone()),
            range,
            options,
            context,
            policy,
        )
    }
}

impl<S> RustFlow<S> {
    /// The exact common-basis physical connection.
    pub fn system(&self) -> &S {
        &self.system
    }
    pub fn identity(&self) -> &BoundaryIdentity {
        &self.identity
    }
}

impl RustFlow<AlgebraicKinematicSystem> {
    pub fn new_algebraic(
        system: AlgebraicKinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
    ) -> Result<Self> {
        Self::with_algebraic_conditions(
            system,
            basis,
            normalization,
            prescription,
            branch_domain,
            &[],
        )
    }
    pub fn with_algebraic_conditions(
        system: AlgebraicKinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        if system.roots.is_empty() {
            return Err(Error::InvalidInput(
                "use the rational RustFlow constructor for a system without registered roots"
                    .into(),
            ));
        }
        let identity = BoundaryIdentity::with_algebraic_conditions(
            &system,
            basis,
            normalization,
            prescription,
            branch_domain,
            conditions,
        )?;
        Ok(Self { system, identity })
    }
    /// Transport within a regular real domain. Every radicand must remain real
    /// and nonzero; the requested local root germ is explicit. The caller's
    /// admissibility policy additionally controls integral/logarithmic monodromy.
    #[allow(clippy::too_many_arguments)] // Ordinary physical inputs plus the required destination root germ.
    pub fn evaluate_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: &RootGerm,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
    ) -> Result<PhysicalResult> {
        evaluate_physical(
            Connection::Algebraic(&self.system),
            &self.identity,
            cache,
            CachedPoint::Exact(destination.clone()).with_root_germ(germ.clone())?,
            range,
            options,
            context,
            policy,
        )
    }
}

impl RustFlow<CanonicalAlgebraicSystem> {
    pub fn new_canonical(
        system: CanonicalAlgebraicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
    ) -> Result<Self> {
        Self::with_canonical_conditions(
            system,
            basis,
            normalization,
            prescription,
            branch_domain,
            &[],
        )
    }
    /// Keep the canonical form separate until a physical path is selected.
    pub fn with_canonical_conditions(
        system: CanonicalAlgebraicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        let identity = BoundaryIdentity::with_canonical_conditions(
            &system,
            basis,
            normalization,
            prescription,
            branch_domain,
            conditions,
        )?;
        Ok(Self { system, identity })
    }
    /// Use an explicit germ when roots are registered, and None otherwise.
    /// The same real regular path and caller monodromy contract as the dense
    /// algebraic interface applies; no dense physical matrix is assembled.
    #[allow(clippy::too_many_arguments)] // Physical transport inputs plus explicit optional root germ.
    pub fn evaluate_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: Option<&RootGerm>,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
    ) -> Result<PhysicalResult> {
        let point = CachedPoint::Exact(destination.clone());
        let point = match (self.system.roots().is_empty(), germ) {
            (true, None) => point,
            (false, Some(germ)) => point.with_root_germ(germ.clone())?,
            _ => {
                return Err(Error::InvalidInput(
                    "canonical root germ must be supplied exactly when roots are registered".into(),
                ));
            }
        };
        evaluate_physical(
            Connection::Canonical(&self.system),
            &self.identity,
            cache,
            point,
            range,
            options,
            context,
            policy,
        )
    }
}

#[derive(Clone, Copy)]
enum Connection<'a> {
    Rational(&'a KinematicSystem),
    Algebraic(&'a AlgebraicKinematicSystem),
    Canonical(&'a CanonicalAlgebraicSystem),
}
impl Connection<'_> {
    fn pullback(&self, path: &KinematicPath, order: usize) -> Result<PreparedConnection> {
        Ok(match self {
            Self::Rational(system) => {
                PreparedConnection::Rational(EpsilonSystem::from_differential_system(
                    &system.pullback(path)?,
                    system.epsilon,
                    order,
                )?)
            }
            Self::Algebraic(system) => PreparedConnection::Algebraic(system.pullback(path, order)?),
            Self::Canonical(system) => PreparedConnection::Algebraic(system.pullback(path, order)?),
        })
    }
}
enum PreparedConnection {
    Rational(EpsilonSystem),
    Algebraic(AlgebraicSystem),
}
impl PreparedConnection {
    fn compile(&self, p: Precision) -> Result<CompiledConnection> {
        Ok(match self {
            Self::Rational(system) => {
                CompiledConnection::Rational(system.compile(p, &Default::default())?)
            }
            Self::Algebraic(system) => CompiledConnection::Algebraic(system.compile(p)?),
        })
    }
    fn transport(
        &self,
        source: &CachedBoundary,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<EpsilonSolution> {
        match self {
            Self::Rational(system) => transport_epsilon(
                system,
                |p| source.as_epsilon_boundary(&Atom::new(), range, p),
                &[Atom::one()],
                options,
                context,
                true,
            ),
            Self::Algebraic(system) => {
                let seeds = if system.roots.is_empty() {
                    BTreeMap::new()
                } else {
                    source
                        .point
                        .root_germ()
                        .ok_or_else(|| Error::InvalidInput("missing source root germ".into()))?
                        .seeds()
                };
                crate::diffexp::refine_epsilon_transport(options, context, true, |p, refined| {
                    let compiled = system.compile(p)?;
                    let result = compiled.transport(
                        &source.as_epsilon_boundary(&Atom::new(), range, p)?,
                        &[p.i(1)],
                        &seeds,
                        refined,
                        context,
                        true,
                    )?;
                    let expected = compiled.branch_state_at(&p.i(1), &seeds)?;
                    if !result
                        .branches
                        .roots
                        .iter()
                        .zip(&expected.roots)
                        .all(|((a, x), (b, y))| a == b && p.close(x, y, options.digits + 8))
                    {
                        return Err(Error::Accuracy(
                            "regular physical transport changed its root germ".into(),
                        ));
                    }
                    Ok(result.solution)
                })
            }
        }
    }
}
enum CompiledConnection {
    Rational(crate::diffexp::CompiledEpsilonSystem),
    Algebraic(CompiledAlgebraicSystem),
}
impl CompiledConnection {
    fn poles(&self) -> &[crate::ComplexFloat] {
        match self {
            Self::Rational(c) => c.poles(),
            Self::Algebraic(c) => c.singularities(),
        }
    }
    fn error_amplification_weighted(
        &self,
        start: &crate::ComplexFloat,
        end: &crate::ComplexFloat,
        weights: &[Float],
    ) -> Result<Float> {
        match self {
            Self::Rational(c) => c.error_amplification_weighted(start, end, weights),
            Self::Algebraic(c) => c.error_amplification_weighted(start, end, weights),
        }
    }
}

#[allow(clippy::too_many_arguments)] // Shared orchestration for both exact connection kinds.
fn evaluate_physical(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    cache: &mut RustFlowCache,
    target: CachedPoint,
    range: EpsilonRange,
    options: &FlowOptions,
    context: &RunContext,
    policy: &dyn TransportCost,
) -> Result<PhysicalResult> {
    options.validate()?;
    context.cancellation.check()?;
    let p = Precision::decimal(options.digits + options.guard_digits)?;
    let query = BoundaryQuery::new(identity, &target, range, options.digits)?;
    // Reduction domains are enforced independently of the caller's sheet
    // policy. Reject unsafe sources before ranking, so another valid source
    // can be selected rather than failing after a nearest-source choice.
    let guarded_policy = GuardedCost {
        identity,
        policy,
        digits: options.digits,
        context,
    };
    let mut excluded = std::collections::BTreeSet::new();
    let mut attempts = Vec::new();
    let mut last_accuracy = None;
    for _ in 0..options.max_boundary_attempts {
        context.cancellation.check()?;
        let Some((index, matched)) = cache.best_excluding(&query, &guarded_policy, p, &excluded)?
        else {
            return Err(if let Some(message) = last_accuracy {
                Error::Accuracy(format!(
                    "none of {} compatible cached boundaries reached the target accuracy; last attempt: {message}",
                    attempts.len()
                ))
            } else {
                Error::IncompleteReduction("no compatible cached physical boundary with the requested epsilon range, verified accuracy and regular reduction domain".into())
            });
        };
        let source = matched.boundary.clone();
        let mut diagnostic = BoundaryAttempt {
            starting_point: source.point.clone(),
            cost: matched.cost,
            source_verified_digits: source.accuracy.verified_digits(),
            outcome: BoundaryAttemptOutcome::Accepted,
        };
        excluded.insert(index);
        context.cancellation.check()?;
        match evaluate_from_source(
            system,
            identity,
            cache,
            target.clone(),
            source,
            range,
            options,
            context,
        ) {
            Ok(mut result) => {
                attempts.push(diagnostic);
                result.boundary_attempts = attempts;
                return Ok(result);
            }
            Err(Error::Accuracy(message)) => {
                diagnostic.outcome = BoundaryAttemptOutcome::AccuracyRejected {
                    message: message.clone(),
                };
                attempts.push(diagnostic);
                last_accuracy = Some(message);
            }
            Err(error) => return Err(error),
        }
    }
    Err(Error::Limit(format!(
        "physical boundary attempt budget {} exhausted after accuracy failures; last attempt: {}",
        options.max_boundary_attempts,
        last_accuracy.unwrap_or_default()
    )))
}

#[allow(clippy::too_many_arguments)] // One already selected source and the shared physical transport inputs.
fn evaluate_from_source(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    cache: &mut RustFlowCache,
    target: CachedPoint,
    source: CachedBoundary,
    range: EpsilonRange,
    options: &FlowOptions,
    context: &RunContext,
) -> Result<PhysicalResult> {
    let p = Precision::decimal(options.digits + options.guard_digits)?;
    let destination = target.restart_coordinates()?;
    let coordinates = source.point.rounded_coordinates_as_exact()?;
    let count = (i64::from(range.last) - i64::from(range.leading) + 1) as usize;
    if source.point.root_germ() == target.root_germ()
        && coordinates
            .iter()
            .all(|(s, a)| (a - &destination[s]).together().cancel().is_zero())
    {
        let mut boundary = source.clone();
        boundary.range = range;
        boundary.coefficients.truncate(count);
        boundary.accuracy = BoundaryAccuracy::supplied(
            source.accuracy.verified_digits(),
            source.accuracy.working_bits(),
            source.accuracy.comparison_errors()[..count].to_vec(),
            "compatible exact-coordinate cache hit",
        )?;
        return Ok(PhysicalResult {
            boundary,
            starting_point: source.point,
            transport: None,
            inserted_points: 0,
            boundary_attempts: Vec::new(),
        });
    }
    let path = KinematicPath::straight_line(
        identity.path_parameter("physical_path_parameter")?,
        &coordinates,
        &destination,
    )?;
    if !identity.conditions_admit_straight_path(&source.point, &target, p, options.digits)? {
        return Err(Error::Unsupported(
            "selected physical path leaves the reduction domain".into(),
        ));
    }
    let system = system.pullback(&path, count - 1)?;
    let prepared = system.compile(p)?;
    if prepared.poles().iter().any(|pole| {
        pole.re >= p.real(0)
            && pole.re <= p.real(1)
            && p.norm(&crate::ComplexFloat::new(p.real(0), pole.im.clone()))
                <= p.tolerance(options.digits + 8)
    }) {
        return Err(Error::Unsupported("selected physical straight path meets a singularity; choose a compatible regular boundary/path or use explicit contour transport".into()));
    }
    let mut solution = system.transport(&source, range, options, context)?;
    let p = Precision {
        bits: solution.diagnostics.working_bits,
    };
    let prepared = system.compile(p)?;
    // A second transport from the same cached values checks integration
    // error only. Carry the supplied boundary uncertainty forward too.
    // Fix one scale per coefficient for the whole trajectory. A large
    // higher epsilon coefficient then contributes only through actual
    // couplings, rather than imposing its absolute error on every order.
    let weights = source.coefficients[..count]
        .iter()
        .flatten()
        .map(|value| {
            let norm = p.norm(value);
            if norm > p.real(1) { norm } else { p.real(1) }
        })
        .collect::<Vec<_>>();
    let mut input_relative_error = p.real(0);
    for (error, weight) in source.accuracy.comparison_errors()[..count]
        .iter()
        .flatten()
        .zip(&weights)
    {
        let relative = error.clone() / weight;
        if relative > input_relative_error {
            input_relative_error = relative;
        }
    }
    input_relative_error += p.tolerance(source.accuracy.verified_digits());
    let evidence_cap = source
        .accuracy
        .verified_digits()
        .min((options.digits + options.guard_digits).saturating_sub(10));
    let mut amplification = p.real(1);
    let mut accepted = Vec::new();
    for (index, segment) in solution.segments.iter().enumerate() {
        context.cancellation.check()?;
        amplification *=
            prepared.error_amplification_weighted(&segment.center, &segment.end, &weights)?;
        if let Some(checkpoint) = solution.checkpoints.iter().find(|c| c.segment == index) {
            let mut checkpoint = checkpoint.clone();
            for (error, weight) in checkpoint
                .comparison_errors
                .iter_mut()
                .flatten()
                .zip(&weights)
            {
                *error += weight.clone() * &input_relative_error * &amplification;
            }
            if errors_meet(
                p,
                &checkpoint.coefficients,
                &checkpoint.comparison_errors,
                options.digits,
            ) {
                accepted.push(checkpoint);
            }
        }
    }
    for (error, weight) in solution
        .comparison_errors
        .iter_mut()
        .flatten()
        .zip(&weights)
    {
        *error += weight.clone() * &input_relative_error * &amplification;
    }
    if !errors_meet(
        p,
        &solution.coefficients,
        &solution.comparison_errors,
        options.digits,
    ) {
        return Err(Error::Accuracy("cached boundary uncertainty grows beyond the requested target accuracy; provide a more accurate or closer boundary".into()));
    }
    solution.verified_digits = Some(strongest_evidence(
        p,
        &solution.coefficients,
        &solution.comparison_errors,
        options.digits,
        evidence_cap,
    ));
    solution.checkpoints = accepted;
    let boundary = CachedBoundary {
        identity: identity.clone(),
        point: target,
        kind: PointKind::Physical,
        range,
        coefficients: solution.coefficients.clone(),
        accuracy: BoundaryAccuracy::from_solution(
            &solution,
            source.accuracy.verified_digits(),
            "independently refined physical transport including propagated cached input uncertainty",
        )?,
    };
    boundary.validate()?;
    let mut pending = RustFlowCache::trajectory_boundaries_with_germ(
        identity,
        &path,
        &solution,
        source.point.root_germ(),
        |index, _| {
            let Some(checkpoint) = solution.checkpoints.iter().find(|c| c.segment == index) else {
                return Ok(None);
            };
            let checked = strongest_evidence(
                p,
                &checkpoint.coefficients,
                &checkpoint.comparison_errors,
                options.digits,
                evidence_cap,
            );
            let accuracy = BoundaryAccuracy::supplied(
                checked,
                p.bits,
                checkpoint.comparison_errors.clone(),
                "independent precision/order check at trajectory endpoint, with propagated cached input uncertainty",
            )?;
            Ok(Some((PointKind::Physical, accuracy)))
        },
    )?;
    let inserted = pending.len() + 1;
    pending.push(boundary.clone());
    cache.insert_many(pending)?;
    Ok(PhysicalResult {
        boundary,
        starting_point: source.point,
        transport: Some(solution),
        inserted_points: inserted,
        boundary_attempts: Vec::new(),
    })
}

fn errors_meet(
    p: Precision,
    values: &[Vec<crate::ComplexFloat>],
    errors: &[Vec<Float>],
    digits: u32,
) -> bool {
    values
        .iter()
        .flatten()
        .zip(errors.iter().flatten())
        .all(|(v, e)| {
            let norm = p.norm(v);
            let scale = if norm > p.real(1) { norm } else { p.real(1) };
            e.is_finite() && *e <= p.tolerance(digits) * scale
        })
}

// Stronger cache evidence requires explicitly passing tighter comparisons and
// propagated input errors. Arithmetic precision supplies only a conservative
// ceiling. Keep two decimal digits beyond any accuracy advertised to the bank.
fn strongest_evidence(
    p: Precision,
    values: &[Vec<crate::ComplexFloat>],
    errors: &[Vec<Float>],
    requested: u32,
    cap: u32,
) -> u32 {
    let mut checked = requested;
    for digits in requested.saturating_add(1)..=cap.saturating_sub(2) {
        if !errors_meet(p, values, errors, digits.saturating_add(2)) {
            break;
        }
        checked = digits;
    }
    checked
}

/// Compatibility name for the physical transport engine.
pub type PhysicalTransport = RustFlow;

struct GuardedCost<'a> {
    identity: &'a BoundaryIdentity,
    policy: &'a dyn TransportCost,
    digits: u32,
    context: &'a RunContext,
}
impl TransportCost for GuardedCost<'_> {
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        self.context.cancellation.check()?;
        if !self
            .identity
            .conditions_admit_straight_path(&source.point, target, p, self.digits)?
        {
            return Ok(None);
        }
        self.policy.cost(source, target, p)
    }
    fn lower_bound(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        self.context.cancellation.check()?;
        self.policy.lower_bound(source, target, p)
    }
    fn compare_tied_costs(
        &self,
        left: &CachedBoundary,
        right: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<std::cmp::Ordering> {
        self.context.cancellation.check()?;
        self.policy.compare_tied_costs(left, right, target, p)
    }
}
