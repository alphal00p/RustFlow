//! Cache-first transport between regular physical kinematic points.
//!
//! Every request selects a compatible boundary before forming its path. Accepted
//! intermediate physical points are independently checked and retained. The
//! current orchestrator supports regular straight paths; threshold prescriptions
//! and partial singular boundaries are handled by separate solver interfaces.
use crate::diffexp::{EpsilonSolution, EpsilonSystem, transport_epsilon};
use crate::kinematics::{KinematicPath, KinematicSystem};
use crate::transport_cache::{
    BoundaryAccuracy, BoundaryIdentity, BoundaryQuery, CachedBoundary, CachedPoint, EpsilonRange,
    PointKind, RustFlowCache, TransportCost,
};
use crate::{Error, FlowOptions, Precision, Prescription, Result, RunContext};
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub struct RustFlow {
    system: KinematicSystem,
    identity: BoundaryIdentity,
}

pub struct PhysicalResult {
    pub boundary: CachedBoundary,
    pub starting_point: CachedPoint,
    /// None for a compatible cache hit at precisely the requested coordinates.
    pub transport: Option<EpsilonSolution>,
    pub inserted_points: usize,
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

    /// The exact common-basis physical connection.
    pub fn system(&self) -> &KinematicSystem {
        &self.system
    }

    pub fn identity(&self) -> &BoundaryIdentity {
        &self.identity
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
        options.validate()?;
        context.cancellation.check()?;
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        let target = CachedPoint::Exact(destination.clone());
        let query = BoundaryQuery::new(&self.identity, &target, range, options.digits)?;
        // Reduction domains are enforced independently of the caller's sheet
        // policy. Reject unsafe sources before ranking, so another valid source
        // can be selected rather than failing after a nearest-source choice.
        let guarded_policy = GuardedCost {
            identity: &self.identity,
            policy,
            digits: options.digits,
        };
        let matched=cache.best(&query,&guarded_policy,p)?.ok_or_else(||Error::IncompleteReduction("no compatible cached physical boundary with the requested epsilon range, verified accuracy and regular reduction domain".into()))?;
        let source = matched.boundary.clone();
        let coordinates = source.point.rounded_coordinates_as_exact()?;
        let count = (i64::from(range.last) - i64::from(range.leading) + 1) as usize;
        if coordinates
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
            });
        }
        let path = KinematicPath::straight_line(
            symbol!("symbolica_amflow::physical_path_parameter"),
            &coordinates,
            destination,
        )?;
        if !self.identity.conditions_admit_straight_path(
            &source.point,
            &target,
            p,
            options.digits,
        )? {
            return Err(Error::Unsupported(
                "selected physical path leaves the reduction domain".into(),
            ));
        }
        let pulled_back = self.system.pullback(&path)?;
        let system =
            EpsilonSystem::from_differential_system(&pulled_back, self.system.epsilon, count - 1)?;
        let prepared = system.compile(p, &Default::default())?;
        if prepared.poles().iter().any(|pole| {
            pole.re >= p.real(0)
                && pole.re <= p.real(1)
                && p.norm(&crate::ComplexFloat::new(p.real(0), pole.im.clone()))
                    <= p.tolerance(options.digits + 8)
        }) {
            return Err(Error::Unsupported("selected physical straight path meets a singularity; choose a compatible regular boundary/path or use explicit contour transport".into()));
        }
        let mut solution = transport_epsilon(
            &system,
            |p| source.as_epsilon_boundary(&Atom::new(), range, p),
            &[Atom::num(1)],
            options,
            context,
            true,
        )?;
        let p = Precision {
            bits: solution.diagnostics.working_bits,
        };
        let prepared = system.compile(p, &Default::default())?;
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
            identity: self.identity.clone(),
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
        let mut pending = RustFlowCache::trajectory_boundaries(
            &self.identity,
            &path,
            &solution,
            |index, _| {
                let Some(checkpoint) = solution.checkpoints.iter().find(|c| c.segment == index)
                else {
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
        })
    }
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
}
impl TransportCost for GuardedCost<'_> {
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        if !self
            .identity
            .conditions_admit_straight_path(&source.point, target, p, self.digits)?
        {
            return Ok(None);
        }
        self.policy.cost(source, target, p)
    }
    fn compare_tied_costs(
        &self,
        left: &CachedBoundary,
        right: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<std::cmp::Ordering> {
        self.policy.compare_tied_costs(left, right, target, p)
    }
}
