//! Bounded, native guarded discovery of a derivative-closed weighted basis.
//!
//! Unreduced integrals are provisional candidates only. A successful result
//! independently audits every target reconstruction and every basis derivative
//! with one source-replayed native program. It makes no minimality claim.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use rustred::solver::guarded::GuardedUnresolvedReason;
use serde::Serialize;
use symbolica::prelude::*;

use super::guarded::{
    GuardedApplicationFailure, GuardedAtomUnresolved, GuardedContext, GuardedDiscovery,
    GuardedDiscoveryOptions, GuardedReductionLimits, GuardedReductionProgram, GuardedUnresolved,
    IndexBounds, IndexDomain, IndexRole, validate_storage_arity, validate_storage_domain,
    validate_storage_label,
};
use crate::reduction::{LinearCombination, ReducedSystem};
use crate::{DifferentialSystem, Error, Integral, Progress, Result, RunContext};

mod active;
mod requested;
#[cfg(test)]
mod requested_tests;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AuxiliaryConvention {
    /// Stored Euclidean factors are rho+t; dI(a)/dt=-a I(a+1).
    EuclideanPlusT,
    /// Stored native Minkowski factors are D-eta; dI(a)/deta=+a I(a+1).
    NativeMinusEta,
}

/// Physical shell/occupation factors are fixed. The family owner must also
/// certify that each selected ordinary slot is an admitted physical quadratic,
/// not a numerator-completion coordinate, and that its source identities use
/// exactly the stated deformation. This algebra service cannot certify a graph.
#[derive(Clone, Debug)]
pub struct FixedShellDeformation<const N: usize> {
    pub variable: Symbol,
    pub convention: AuxiliaryConvention,
    roles: [IndexRole; N],
    shifted: [bool; N],
    admitted: IndexDomain<N>,
    physical_arity: usize,
}

impl<const N: usize> FixedShellDeformation<N> {
    pub fn new(
        variable: Symbol,
        convention: AuxiliaryConvention,
        roles: [IndexRole; N],
        shifted: [bool; N],
    ) -> Result<Self> {
        if !shifted.iter().any(|selected| *selected) {
            return Err(Error::Unsupported(
                "no admitted uncut factor flows; a uniform compact terminal is required".into(),
            ));
        }
        if shifted
            .iter()
            .zip(&roles)
            .any(|(selected, role)| *selected && *role != IndexRole::Ordinary)
        {
            return Err(Error::InvalidInput(
                "fixed-shell flow cannot deform required cuts or occupation factors".into(),
            ));
        }
        Ok(Self {
            variable,
            convention,
            roles,
            shifted,
            admitted: IndexDomain::for_roles(&roles),
            physical_arity: N,
        })
    }

    /// Declare how many coordinates represent actual physical factors. The
    /// remaining storage slots have no factors or deformation and stay zero.
    pub fn with_physical_arity(mut self, physical_arity: usize) -> Result<Self> {
        validate_storage_arity::<N>(physical_arity)?;
        if self.roles[physical_arity..]
            .iter()
            .any(|role| *role != IndexRole::Ordinary)
            || self.shifted[physical_arity..]
                .iter()
                .any(|&shifted| shifted)
        {
            return Err(Error::InvalidInput(
                "weighted deformation uses a storage tail coordinate".into(),
            ));
        }
        if self.admitted.bounds()[physical_arity..]
            .iter()
            .any(|bound| !bound.contains(0))
        {
            return Err(Error::InvalidInput(
                "weighted admitted domain excludes zero storage tail".into(),
            ));
        }
        self.admitted = IndexDomain::new(std::array::from_fn(|axis| {
            if axis < physical_arity {
                self.admitted.bounds()[axis]
            } else {
                IndexBounds::fixed(0)
            }
        }))
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
        self.physical_arity = physical_arity;
        Ok(self)
    }

    pub fn physical_arity(&self) -> usize {
        self.physical_arity
    }

    /// Restrict valid integral labels independently of the native source guards.
    /// Invalid labels must never become provisional or constant terminal terms.
    pub fn with_admitted_domain(mut self, domain: IndexDomain<N>) -> Result<Self> {
        validate_storage_domain(&domain, self.physical_arity)?;
        if !domain.is_subset_of(&IndexDomain::for_roles(&self.roles)) {
            return Err(Error::InvalidInput(
                "weighted admitted domain includes undefined occupation indices".into(),
            ));
        }
        self.admitted = domain;
        Ok(self)
    }

    pub fn admitted_domain(&self) -> &IndexDomain<N> {
        &self.admitted
    }

    fn validate_integral(&self, integral: &[i64; N]) -> Result<()> {
        validate_storage_label(integral, self.physical_arity)?;
        if !self.admitted.contains(integral) {
            return Err(Error::InvalidInput(format!(
                "weighted integral {integral:?} lies outside the admitted index domain"
            )));
        }
        Ok(())
    }

    pub fn shifted(&self) -> &[bool; N] {
        &self.shifted
    }

    pub fn derivative(&self, integral: [i64; N]) -> Result<BTreeMap<[i64; N], Atom>> {
        self.validate_integral(&integral)?;
        let mut result = BTreeMap::new();
        for (slot, (&power, &selected)) in integral.iter().zip(&self.shifted).enumerate() {
            if !selected || power == 0 {
                continue;
            }
            let mut raised = integral;
            raised[slot] = raised[slot].checked_add(1).ok_or_else(|| {
                Error::Limit("weighted auxiliary derivative index overflow".into())
            })?;
            self.validate_integral(&raised)?;
            let coefficient = match self.convention {
                AuxiliaryConvention::EuclideanPlusT => -Atom::num(power),
                AuxiliaryConvention::NativeMinusEta => Atom::num(power),
            };
            result.insert(raised, coefficient);
        }
        Ok(result)
    }
}

/// Subdivide only finite native proof gaps reached by actual reductions. These
/// limits control additional native searches, never an alternative elimination.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct GuardRefinementOptions {
    /// Additional discovery passes per provisional or final-audit program.
    /// Zero preserves discovery without frontend interval subdivision.
    pub max_passes: usize,
    /// Total distinct added faces over the complete closure preparation.
    pub max_added_domains: usize,
    /// Inclusive integer interval width (upper minus lower), not point count.
    pub max_interval_width: u64,
}

impl Default for GuardRefinementOptions {
    fn default() -> Self {
        Self {
            max_passes: 3,
            max_added_domains: 256,
            max_interval_width: 2,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct GuardRefinementRecord {
    pub pass: usize,
    pub parent_bounds: Vec<[Option<i64>; 2]>,
    pub native_failure: String,
    pub axis: usize,
    pub fixed_values: Vec<i64>,
    pub trigger: Vec<i64>,
}

/// One concrete conditional index face submitted to native discovery. This is
/// search evidence, not a replacement for a successful native application.
#[derive(Clone, Debug, Serialize)]
pub struct ConditionalPointRefinement {
    /// Zero-based pass within one bounded discovery call.
    pub pass: usize,
    pub integral: Vec<i64>,
    pub allocated_domains: usize,
}

/// Bounded native searches at each explicitly requested label. A symbolic ray
/// is followed by an exact-point search only when its native application misses.
/// This schedules the existing source-replayed reducer; it changes no identities.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct RequestedRayPointOptions {
    pub max_domains_per_ray: usize,
    pub max_domains_per_point: usize,
}

#[derive(Clone, Debug)]
pub struct WeightedClosureOptions {
    pub max_rounds: usize,
    pub max_frontier: usize,
    pub max_requested: usize,
    pub discovery: GuardedDiscoveryOptions,
    pub application: GuardedReductionLimits,
    pub guard_refinement: GuardRefinementOptions,
    /// Request every provisional frontier label as well as derivatives, so
    /// auxiliary-constant lower sectors receive native source discovery.
    pub search_frontier_sectors: bool,
    /// Keep requested ordinary zero indices on fixed-zero faces, and negative
    /// indices on strictly negative boxes. This narrows native discovery only;
    /// it does not classify a physical integral as zero or alter the ordering.
    pub split_ordinary_zero_faces: bool,
    /// Start each nonzero symbolic ray at the requested integer, with zero
    /// indices fixed. Narrows discovery only; all resulting native rules still
    /// carry their exact domains and must reconstruct the full connection.
    pub requested_index_rays: bool,
    /// Ask native RustRed to explore domains containing requested labels first.
    /// This only schedules bounded work; unvisited domains remain explicit gaps.
    pub prioritize_requested_indices: bool,
    /// Retain replayed native rules across refinement passes and closure
    /// rounds, with fresh rules before older fallbacks. Zero disables reuse;
    /// otherwise this bounds the whole union, including duplicate rules.
    pub max_reused_rules: usize,
    /// Opt-in residual-focused discovery. Each actual NoApplicableRule leaf
    /// receives this many allocated domain visits, sharing discovery.max_domains
    /// across refinement passes. ConditionVanished leaves receive exact-point
    /// searches from the same allocation. Requires requested rays and retained rules.
    /// Existing rules take precedence so adding rules preserves old rewrites.
    /// Zero keeps the original multi-domain search policy.
    pub max_domains_per_residual: usize,
    /// Optional explicit-request ray/point schedule, with empty stopping sets.
    /// Requires active closure and retained rules. Mutually exclusive with
    /// residual allocation, interval refinement and requested-priority hints.
    pub requested_ray_point: Option<RequestedRayPointOptions>,
    /// Rebuild the needed integral module from original target sums after each
    /// native rule update. Historical requests remain bounded evidence, rather
    /// than additional targets. Requires a retained residual or requested-point schedule.
    pub active_target_closure: bool,
    /// Native one-source zero proofs at requested concrete points before
    /// ordinary discovery. A per-call attempt cap; zero disables the prepass.
    /// Requires a retained discovery schedule. Completed point searches are memoized.
    pub max_direct_zero_attempts: usize,
    /// A distinct directory per source/deformation context. Each round stores
    /// its exact native program plus explicitly provisional/closed metadata.
    pub checkpoints: Option<PathBuf>,
}

impl Default for WeightedClosureOptions {
    fn default() -> Self {
        Self {
            max_rounds: 12,
            max_frontier: 256,
            max_requested: 4096,
            discovery: GuardedDiscoveryOptions {
                max_depth: 3,
                max_domains: 8192,
                ..Default::default()
            },
            application: GuardedReductionLimits::default(),
            guard_refinement: GuardRefinementOptions::default(),
            search_frontier_sectors: true,
            split_ordinary_zero_faces: false,
            requested_index_rays: false,
            prioritize_requested_indices: false,
            max_reused_rules: 0,
            max_domains_per_residual: 0,
            requested_ray_point: None,
            active_target_closure: false,
            max_direct_zero_attempts: 0,
            checkpoints: None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct WeightedClosureDiagnostics {
    pub rounds: usize,
    pub requested: usize,
    /// Unique labels submitted across all rounds, including retired requests.
    pub historical_requested: usize,
    /// Historical labels that the final program does not individually reduce;
    /// these are not candidates or masters in an active-target closed result.
    pub retired_unresolved: usize,
    pub direct_zero_calls: usize,
    pub direct_zero_attempts: usize,
    pub direct_zero_rules: usize,
    /// Completed bounded point searches, not zero or master classifications.
    pub direct_zero_completed_points: BTreeSet<Vec<i64>>,
    /// Provisional labels explicitly added for native source discovery, including
    /// auxiliary-constant lower sectors that derivatives cannot request.
    pub native_frontier_requests: usize,
    pub provisional_sizes: Vec<usize>,
    pub native_rules: usize,
    pub native_rule_unions: usize,
    pub residual_search_requests: usize,
    /// Conservative allocations, which can exceed actual native domain visits.
    pub residual_allocated_domains: usize,
    pub requested_ray_point_calls: usize,
    /// Conservatively allocated caps, not measured native visits.
    pub requested_ray_point_allocated_domains: usize,
    pub requested_ray_point_fallbacks: usize,
    pub requested_ray_point_deferred: usize,
    /// Exact-point searches for native vanishing guards, separate from finite
    /// axis subdivisions. Failed searches never install a terminal or zero.
    pub conditional_point_refinements: Vec<ConditionalPointRefinement>,
    pub native_rule_applications: usize,
    /// Discovery need not solve every native index domain to close the finite
    /// requested target and derivative set. Under residual-focused scheduling,
    /// this counts the deduplicated historical gaps, not just the final search.
    pub uncovered_discovery_domains: usize,
    pub historical_discovery_domains: usize,
    pub guard_refinement_passes: usize,
    pub guard_refinement_added_domains: usize,
    pub guard_refinement_budget_exhausted: bool,
    /// Historical parent gaps remain evidence of incomplete symbolic coverage;
    /// adding faces alone is never a claim that the parent has been certified.
    pub guard_refinements: Vec<GuardRefinementRecord>,
}

#[derive(Debug)]
pub struct WeightedReducedSystem<const N: usize> {
    pub reduced: ReducedSystem,
    pub variable: Symbol,
    pub roles: [IndexRole; N],
    pub program: GuardedReductionProgram<N>,
    pub discovery_unresolved: Vec<GuardedUnresolved<N>>,
    pub diagnostics: WeightedClosureDiagnostics,
}

impl<const N: usize> WeightedReducedSystem<N> {
    pub fn physical_arity(&self) -> usize {
        self.program.physical_arity()
    }

    /// An identically zero target collection needs no fictitious ODE component.
    pub fn differential_system(&self) -> Option<DifferentialSystem> {
        (!self.reduced.basis.is_empty()).then(|| DifferentialSystem {
            variable: self.variable,
            matrix: self.reduced.matrix.clone(),
        })
    }
}

#[derive(Debug)]
pub struct UnclosedWeightedSystem<const N: usize> {
    pub reason: String,
    pub provisional_frontier: Vec<[i64; N]>,
    pub unresolved: Vec<GuardedAtomUnresolved<N>>,
    pub discovery_unresolved: Vec<GuardedUnresolved<N>>,
    pub nonzero_conditions: Vec<Atom>,
    pub diagnostics: WeightedClosureDiagnostics,
}

#[derive(Debug)]
pub enum WeightedClosureOutcome<const N: usize> {
    Closed(WeightedReducedSystem<N>),
    Unresolved(UnclosedWeightedSystem<N>),
}

/// Build a closed spanning system from input targets, without a caller-selected
/// terminal basis. Input target coefficients are fixed along the auxiliary flow;
/// eta-dependent coefficients produced by native target reduction are retained.
pub fn prepare_weighted_system<const N: usize>(
    context: &GuardedContext<N>,
    targets: &[BTreeMap<[i64; N], Atom>],
    deformation: &FixedShellDeformation<N>,
    options: WeightedClosureOptions,
    run: &RunContext,
) -> Result<WeightedClosureOutcome<N>> {
    if targets.is_empty()
        || options.max_rounds == 0
        || options.max_frontier == 0
        || options.max_requested == 0
    {
        return Err(Error::InvalidInput(
            "weighted closure needs targets and positive work limits".into(),
        ));
    }
    if options.guard_refinement.max_passes > 0
        && (options.guard_refinement.max_added_domains < 2
            || options.guard_refinement.max_interval_width == 0)
    {
        return Err(Error::InvalidInput(
            "enabled weighted guard refinement needs at least two faces and positive interval width".into(),
        ));
    }
    if options.max_domains_per_residual > 0
        && (options.max_reused_rules == 0 || !options.requested_index_rays)
    {
        return Err(Error::InvalidInput(
            "residual-focused discovery requires retained rules and requested index rays".into(),
        ));
    }
    if let Some(schedule) = options.requested_ray_point {
        if schedule.max_domains_per_ray == 0
            || schedule.max_domains_per_point == 0
            || schedule
                .max_domains_per_ray
                .checked_add(schedule.max_domains_per_point)
                .is_none()
            || options.discovery.max_domains == 0
            || !options.active_target_closure
            || !options.search_frontier_sectors
            || !options.requested_index_rays
            || options.max_reused_rules == 0
            || options.max_domains_per_residual != 0
            || options.guard_refinement.max_passes != 0
            || options.prioritize_requested_indices
        {
            return Err(Error::InvalidInput(
                "requested ray/point discovery needs positive bounded allocations, active frontier closure, requested rays and retained rules; residual allocation, guard refinement and priority hints must be disabled".into(),
            ));
        }
    }
    if options.active_target_closure
        && ((options.max_domains_per_residual == 0 && options.requested_ray_point.is_none())
            || !options.search_frontier_sectors)
    {
        return Err(Error::InvalidInput(
            "active target closure requires a retained discovery schedule and frontier searches"
                .into(),
        ));
    }
    if options.max_direct_zero_attempts > 0
        && options.max_domains_per_residual == 0
        && options.requested_ray_point.is_none()
    {
        return Err(Error::InvalidInput(
            "direct zero discovery requires a retained discovery schedule".into(),
        ));
    }
    if context.physical_arity() != deformation.physical_arity {
        return Err(Error::InvalidInput(
            "weighted source/deformation physical arities differ".into(),
        ));
    }
    if context.sources().roles() != &deformation.roles {
        return Err(Error::InvalidInput(
            "weighted source/deformation index roles differ".into(),
        ));
    }
    if !context
        .sources()
        .native_sources()
        .coefficient_variables()
        .iter()
        .skip(N)
        .any(|variable| variable == &PolyVariable::Symbol(deformation.variable))
    {
        return Err(Error::InvalidInput(
            "auxiliary variable must be a declared physical parameter of the weighted source context".into(),
        ));
    }
    if targets
        .iter()
        .flat_map(|target| target.values())
        .any(|coefficient| coefficient.contains_symbol(deformation.variable))
    {
        return Err(Error::InvalidInput(
            "input target coefficients must stay fixed along auxiliary flow".into(),
        ));
    }
    let physical_variables = context
        .sources()
        .native_sources()
        .coefficient_variables()
        .iter()
        .skip(N)
        .map(|variable| match variable {
            PolyVariable::Symbol(symbol) => Ok(*symbol),
            _ => Err(Error::InvalidInput(
                "weighted coefficient context contains a non-symbol parameter".into(),
            )),
        })
        .collect::<Result<BTreeSet<_>>>()?;
    let mut conditions = BTreeMap::new();
    for target in targets {
        for indices in target.keys() {
            deformation.validate_integral(indices)?;
        }
        // Validate every original coefficient before cancellation and before
        // cut/support zero filtering can remove its integral from the system.
        retain_conditions(
            &mut conditions,
            crate::physical_conditions::rational_denominator_conditions(
                &target.values().cloned().collect::<Vec<_>>(),
                &physical_variables,
            )?,
        );
    }
    if options.active_target_closure {
        return active::prepare(context, targets, deformation, options, conditions, run);
    }
    let mut requested = targets
        .iter()
        .flat_map(|target| target.iter())
        .filter(|(_, coefficient)| !coefficient.is_zero())
        .map(|(indices, _)| *indices)
        .collect::<BTreeSet<_>>();
    let mut diagnostics = WeightedClosureDiagnostics::default();
    let mut last_frontier = Vec::new();
    let mut last_unresolved = Vec::new();
    let mut last_discovery = Vec::new();
    let mut refined_domains = Vec::new();
    let mut retained_program = None;
    let mut historical_gaps = BTreeMap::new();
    for round in 0..options.max_rounds {
        run.cancellation.check()?;
        if requested.len() > options.max_requested {
            return Ok(unclosed(
                "weighted requested-integral budget exhausted",
                last_frontier,
                last_unresolved,
                last_discovery,
                conditions,
                diagnostics,
            ));
        }
        let domains = discovery_domains(
            &requested,
            &deformation.roles,
            options.split_ordinary_zero_faces,
            options.requested_index_rays,
        )?
        .into_iter()
        .filter_map(|domain| domain.intersection(deformation.admitted_domain()))
        .collect::<Vec<_>>();
        let found = discover_with_refinement(
            context,
            &domains,
            &[],
            &requested,
            deformation,
            &options,
            &mut refined_domains,
            &mut retained_program,
            &mut historical_gaps,
            &mut diagnostics,
            &mut conditions,
            run,
        )?;
        diagnostics.rounds = round + 1;
        diagnostics.requested = requested.len();
        diagnostics.historical_requested = requested.len();
        diagnostics.native_rules = found.program.native().rules().len();
        diagnostics.uncovered_discovery_domains = found.unresolved.len();
        let mut frontier = BTreeSet::new();
        let mut failures = Vec::new();
        let mut residuals = Vec::new();
        for target in &requested {
            deformation.validate_integral(target)?;
            run.cancellation.check()?;
            let reduced = found.program.reduce(*target, options.application)?;
            diagnostics.native_rule_applications += reduced.rule_applications;
            retain_conditions(&mut conditions, reduced.nonzero_conditions);
            for indices in reduced.terms.keys() {
                deformation.validate_integral(indices)?;
                frontier.insert(*indices);
            }
            for unresolved in reduced.unresolved {
                deformation.validate_integral(&unresolved.integral)?;
                if unresolved.reason == GuardedApplicationFailure::NoApplicableRule {
                    frontier.insert(unresolved.integral);
                    residuals.push(unresolved);
                } else {
                    failures.push(unresolved);
                }
            }
        }
        let frontier = frontier.into_iter().collect::<Vec<_>>();
        diagnostics.provisional_sizes.push(frontier.len());
        checkpoint(
            &options,
            round,
            false,
            &found.program,
            &requested,
            &frontier,
            None,
            run,
        )?;
        if !failures.is_empty() {
            return Ok(unclosed(
                "native weighted application failed",
                frontier,
                failures,
                found.unresolved,
                conditions,
                diagnostics,
            ));
        }
        if frontier.len() > options.max_frontier {
            return Ok(unclosed(
                "weighted provisional-frontier budget exhausted",
                frontier,
                residuals,
                found.unresolved,
                conditions,
                diagnostics,
            ));
        }
        let derivatives = frontier
            .iter()
            .map(|integral| deformation.derivative(*integral))
            .collect::<Result<Vec<_>>>()?;
        let new_derivatives = derivatives
            .iter()
            .flat_map(|terms| terms.keys())
            .filter(|indices| !requested.contains(*indices))
            .copied()
            .collect::<BTreeSet<_>>();
        run.emit(Progress::DifferentialClosure {
            round,
            requested: requested.len(),
            basis_size: frontier.len(),
            new_derivatives: new_derivatives.len(),
        })?;
        let new_frontier = frontier
            .iter()
            .filter(|indices| options.search_frontier_sectors && !requested.contains(*indices))
            .copied()
            .collect::<BTreeSet<_>>();
        diagnostics.native_frontier_requests += new_frontier.len();
        if !new_frontier.is_empty() {
            run.emit(Progress::Stage {
                name: format!(
                    "searching {} new native weighted frontier labels, including constant sectors",
                    new_frontier.len()
                ),
            })?;
        }
        if new_derivatives.is_empty() && new_frontier.is_empty() {
            if options.max_reused_rules > 0 {
                retained_program = Some(found.program);
            }
            // Rebuild one source-replayed program whose explicit stopping set
            // is now audited against every derivative and target. This is the
            // only point at which provisional candidates become a closed basis.
            let final_found = discover_with_refinement(
                context,
                &domains,
                &frontier,
                &requested,
                deformation,
                &options,
                &mut refined_domains,
                &mut retained_program,
                &mut historical_gaps,
                &mut diagnostics,
                &mut conditions,
                run,
            )?;
            diagnostics.native_rules = final_found.program.native().rules().len();
            diagnostics.uncovered_discovery_domains = final_found.unresolved.len();
            let mut reductions = BTreeMap::new();
            let mut failed = Vec::new();
            for target in &requested {
                deformation.validate_integral(target)?;
                let reduced = final_found.program.reduce(*target, options.application)?;
                diagnostics.native_rule_applications += reduced.rule_applications;
                retain_conditions(&mut conditions, reduced.nonzero_conditions);
                for residual in &reduced.unresolved {
                    deformation.validate_integral(&residual.integral)?;
                }
                failed.extend(reduced.unresolved);
                for indices in reduced.terms.keys() {
                    deformation.validate_integral(indices)?;
                }
                reductions.insert(*target, reduced.terms);
            }
            if !failed.is_empty() {
                return Ok(unclosed(
                    "final weighted closure audit left unresolved terms",
                    frontier,
                    failed,
                    final_found.unresolved,
                    conditions,
                    diagnostics,
                ));
            }
            let mut matrix = vec![vec![Atom::new(); frontier.len()]; frontier.len()];
            for (row, derivative) in derivatives.iter().enumerate() {
                for (target, factor) in derivative {
                    for (index, coefficient) in &reductions[target] {
                        let column = frontier.binary_search(index).map_err(|_| {
                            Error::IncompleteReduction(
                                "audited weighted derivative escaped its basis".into(),
                            )
                        })?;
                        matrix[row][column] += factor * coefficient;
                    }
                }
                for coefficient in &mut matrix[row] {
                    *coefficient = coefficient.together().cancel();
                }
            }
            let target_weights = targets
                .iter()
                .map(|target| {
                    let mut weights = BTreeMap::new();
                    for (indices, weight) in target {
                        if weight.is_zero() {
                            continue;
                        }
                        for (basis, coefficient) in &reductions[indices] {
                            *weights.entry(*basis).or_insert_with(Atom::new) +=
                                weight * coefficient;
                        }
                    }
                    linear_combination(&weights, deformation.physical_arity)
                })
                .collect::<Result<Vec<_>>>()?;
            let candidates = reductions
                .iter()
                .map(|(indices, terms)| {
                    Ok((
                        native_integral(indices, deformation.physical_arity)?,
                        linear_combination(terms, deformation.physical_arity)?,
                    ))
                })
                .collect::<Result<BTreeMap<_, _>>>()?;
            checkpoint(
                &options,
                round,
                true,
                &final_found.program,
                &requested,
                &frontier,
                None,
                run,
            )?;
            return Ok(WeightedClosureOutcome::Closed(WeightedReducedSystem {
                reduced: ReducedSystem {
                    basis: frontier
                        .iter()
                        .map(|indices| native_integral(indices, deformation.physical_arity))
                        .collect::<Result<_>>()?,
                    matrix,
                    targets: target_weights,
                    nonzero_conditions: conditions.into_values().collect(),
                    candidates,
                    transformations: Vec::new(),
                },
                variable: deformation.variable,
                roles: deformation.roles,
                program: final_found.program,
                discovery_unresolved: final_found.unresolved,
                diagnostics,
            }));
        }
        requested.extend(new_derivatives);
        requested.extend(new_frontier);
        last_frontier = frontier;
        last_unresolved = residuals;
        last_discovery = found.unresolved;
        if options.max_reused_rules > 0 {
            retained_program = Some(found.program);
        }
    }
    Ok(unclosed(
        "weighted derivative closure round budget exhausted",
        last_frontier,
        last_unresolved,
        last_discovery,
        conditions,
        diagnostics,
    ))
}

/// Native discovery and application remain the sole proof and reduction owners.
/// The frontend only changes which exact integer boxes are submitted to them.
type HistoricalDiscoveryGaps<const N: usize> =
    BTreeMap<([[Option<i64>; 2]; N], String), GuardedUnresolved<N>>;

#[allow(clippy::too_many_arguments)]
fn discover_direct_zeros<const N: usize>(
    context: &GuardedContext<N>,
    terminals: &[[i64; N]],
    requested: &BTreeSet<[i64; N]>,
    deformation: &FixedShellDeformation<N>,
    options: &WeightedClosureOptions,
    retained: &mut Option<GuardedReductionProgram<N>>,
    historical_gaps: &mut HistoricalDiscoveryGaps<N>,
    diagnostics: &mut WeightedClosureDiagnostics,
    run: &RunContext,
) -> Result<()> {
    let current_terminals = terminals.iter().copied().collect::<BTreeSet<_>>();
    if options.max_direct_zero_attempts > 0 {
        let points = requested
            .iter()
            .filter(|point| {
                !current_terminals.contains(*point)
                    && !diagnostics
                        .direct_zero_completed_points
                        .contains(point.as_slice())
            })
            .copied()
            .collect::<Vec<_>>();
        for point in &points {
            deformation.validate_integral(point)?;
        }
        if !points.is_empty() {
            run.cancellation.check()?;
            let zero = context.discover_direct_zeros(
                &points,
                terminals.iter().copied(),
                options.max_direct_zero_attempts,
            )?;
            run.cancellation.check()?;
            let rule_count = zero.discovery.program.native().rules().len();
            if rule_count > options.max_reused_rules {
                return Err(Error::Limit(
                    "native direct-zero rules exceed the retained-rule budget".into(),
                ));
            }
            let program = if let Some(previous) = retained.take() {
                // A replayed zero can replace an expanding recurrence. All
                // subsequent ordinary rules retain this zero-first precedence.
                diagnostics.native_rule_unions += 1;
                zero.discovery.program.union_verified(
                    previous,
                    terminals.iter().copied(),
                    options.max_reused_rules,
                )?
            } else {
                zero.discovery.program
            };
            *retained = Some(program);
            diagnostics.direct_zero_calls += 1;
            diagnostics.direct_zero_attempts += zero.attempted_rows;
            diagnostics.direct_zero_rules += rule_count;
            diagnostics
                .direct_zero_completed_points
                .extend(zero.completed_points.iter().map(|p| p.to_vec()));
            for gap in zero.discovery.unresolved {
                let bounds = std::array::from_fn(|axis| {
                    let bound = gap.domain.bounds()[axis];
                    [bound.lower(), bound.upper()]
                });
                historical_gaps.insert((bounds, format!("{:?}:{}", gap.reason, gap.detail)), gap);
            }
            diagnostics.historical_discovery_domains = historical_gaps.len();
            if let Some(directory) = &options.checkpoints {
                std::fs::create_dir_all(directory)?;
                let path = directory.join(format!(
                    "direct-zero-{:04}.json",
                    diagnostics.direct_zero_calls
                ));
                let temporary = path.with_extension("json.part");
                std::fs::write(&temporary, serde_json::to_vec_pretty(&serde_json::json!({
                    "schema":1, "attempt_budget":options.max_direct_zero_attempts,
                    "attempted_rows":zero.attempted_rows, "zero_rules":rule_count,
                    "requested_points":points.iter().map(|p| p.to_vec()).collect::<Vec<_>>(),
                    "completed_points":zero.completed_points.iter().map(|p| p.to_vec()).collect::<Vec<_>>(),
                    "completed_search_does_not_imply_zero":true, "skipped_seeds":zero.skipped_seeds,
                })).map_err(|error| Error::Cache(error.to_string()))?)?;
                std::fs::rename(temporary, path)?;
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn discover_with_refinement<const N: usize>(
    context: &GuardedContext<N>,
    domains: &[IndexDomain<N>],
    terminals: &[[i64; N]],
    requested: &BTreeSet<[i64; N]>,
    deformation: &FixedShellDeformation<N>,
    options: &WeightedClosureOptions,
    learned: &mut Vec<IndexDomain<N>>,
    retained: &mut Option<GuardedReductionProgram<N>>,
    historical_gaps: &mut HistoricalDiscoveryGaps<N>,
    diagnostics: &mut WeightedClosureDiagnostics,
    conditions: &mut BTreeMap<String, Atom>,
    run: &RunContext,
) -> Result<GuardedDiscovery<N>> {
    let residual_policy = options.max_domains_per_residual > 0;
    let mut remaining_domains = options.discovery.max_domains;
    // The source corpus and point-search options are immutable during this
    // call. An unsuccessful exact-point search cannot gain coverage by being
    // repeated in every refinement pass.
    let mut attempted_guard_points = BTreeSet::new();
    let current_terminals = terminals.iter().copied().collect::<BTreeSet<_>>();
    discover_direct_zeros(
        context,
        terminals,
        requested,
        deformation,
        options,
        retained,
        historical_gaps,
        diagnostics,
        run,
    )?;
    for pass in 0..=options.guard_refinement.max_passes {
        run.cancellation.check()?;
        let (search_requests, guard_points) = if residual_policy {
            let program = if let Some(previous) = retained.take() {
                if previous.native().terminals() == &current_terminals {
                    previous
                } else {
                    previous.with_terminals_verified(
                        terminals.iter().copied(),
                        options.max_reused_rules,
                    )?
                }
            } else {
                context
                    .discover(Vec::new(), terminals.iter().copied(), options.discovery)?
                    .program
            };
            let (residuals, guard_points, failed) = residual_search_requests(
                &program,
                requested,
                deformation,
                options.application,
                diagnostics,
                conditions,
                run,
            )?;
            if failed {
                // The caller repeats native application and records the actual
                // failure. No failed leaf is silently promoted to a terminal.
                return Ok(GuardedDiscovery {
                    program,
                    unresolved: historical_gaps.values().cloned().collect(),
                });
            }
            *retained = Some(program);
            diagnostics.residual_search_requests += residuals.len() + guard_points.len();
            (residuals, guard_points)
        } else {
            (requested.clone(), BTreeSet::new())
        };
        let priority_points = if options.prioritize_requested_indices {
            search_requests
                .union(&guard_points)
                .copied()
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let new_guard_points = guard_points
            .difference(&attempted_guard_points)
            .copied()
            .collect::<Vec<_>>();
        let mut submitted = if residual_policy {
            // A coupled exceptional locus (for example a_i=a_j) need not be
            // representable by a coordinate ray. Its actual failed integer
            // label is an exact box and must not be broadened back into a ray.
            let mut points = new_guard_points
                .iter()
                .map(|point| {
                    IndexDomain::new(point.map(IndexBounds::fixed))
                        .map_err(|error| Error::InvalidInput(error.to_string()))
                })
                .collect::<Result<Vec<_>>>()?;
            points.extend(
                discovery_domains(&search_requests, &deformation.roles, true, true)?
                    .into_iter()
                    .filter_map(|domain| domain.intersection(deformation.admitted_domain())),
            );
            points
        } else {
            domains.to_vec()
        };
        let initial_domains = submitted.clone();
        submitted.extend(
            learned
                .iter()
                .filter(|face| {
                    !initial_domains.contains(face)
                        && (!residual_policy
                            || search_requests.iter().any(|point| face.contains(point)))
                })
                .cloned(),
        );
        let mut found = if residual_policy {
            let mut guard_budget = remaining_domains;
            for point in &new_guard_points {
                let allocated = options.max_domains_per_residual.min(guard_budget);
                if allocated == 0 {
                    break;
                }
                guard_budget -= allocated;
                attempted_guard_points.insert(*point);
                diagnostics
                    .conditional_point_refinements
                    .push(ConditionalPointRefinement {
                        pass,
                        integral: point.to_vec(),
                        allocated_domains: allocated,
                    });
            }
            let allocated = options
                .max_domains_per_residual
                .saturating_mul(submitted.len())
                .min(remaining_domains);
            let found = context.discover_partitioned(
                submitted,
                terminals.iter().copied(),
                &priority_points,
                GuardedDiscoveryOptions {
                    max_domains: remaining_domains,
                    ..options.discovery
                },
                options.max_domains_per_residual,
            )?;
            remaining_domains -= allocated;
            diagnostics.residual_allocated_domains += allocated;
            found
        } else {
            context.discover_with_priority_points(
                submitted,
                terminals.iter().copied(),
                &priority_points,
                options.discovery,
            )?
        };
        if options.max_reused_rules > 0 {
            if found.program.native().rules().len() > options.max_reused_rules {
                return Err(Error::Limit(
                    "fresh native rule program exceeds the retained-rule budget".into(),
                ));
            }
            if let Some(previous) = retained.take() {
                found.program = if residual_policy {
                    // Fill uncovered leaves, including exact conditional
                    // faces, without diverting successful old reductions.
                    previous.union_verified(
                        found.program,
                        terminals.iter().copied(),
                        options.max_reused_rules,
                    )?
                } else {
                    found.program.union_verified(
                        previous,
                        terminals.iter().copied(),
                        options.max_reused_rules,
                    )?
                };
                diagnostics.native_rule_unions += 1;
            }
        }
        if residual_policy {
            for gap in found.unresolved {
                let bounds = std::array::from_fn(|axis| {
                    let bound = gap.domain.bounds()[axis];
                    [bound.lower(), bound.upper()]
                });
                historical_gaps.insert((bounds, format!("{:?}:{}", gap.reason, gap.detail)), gap);
            }
            found.unresolved = historical_gaps.values().cloned().collect();
            diagnostics.historical_discovery_domains = found.unresolved.len();
            if let Some(directory) = options.checkpoints.as_deref() {
                std::fs::create_dir_all(directory)?;
                let gaps = found
                    .unresolved
                    .iter()
                    .map(|gap| {
                        serde_json::json!({
                            "bounds":gap.domain.bounds().iter().map(|bound|
                                [bound.lower(),bound.upper()]).collect::<Vec<_>>(),
                            "reason":format!("{:?}",gap.reason), "detail":gap.detail,
                        })
                    })
                    .collect::<Vec<_>>();
                std::fs::write(
                    directory.join("historical-discovery-gaps.json"),
                    serde_json::to_vec_pretty(&serde_json::json!({
                        "schema":1, "scope":"historical native gaps, not current target coverage",
                        "whole_domain_coverage_certified":false, "gaps":gaps,
                    }))
                    .map_err(|error| Error::Cache(error.to_string()))?,
                )?;
                std::fs::write(
                    directory.join("conditional-point-refinements.json"),
                    serde_json::to_vec_pretty(&serde_json::json!({
                        "schema":1,
                        "scope":"bounded native exact-point searches for ConditionVanished leaves",
                        "search_does_not_imply_coverage":true,
                        "refinements":diagnostics.conditional_point_refinements,
                    }))
                    .map_err(|error| Error::Cache(error.to_string()))?,
                )?;
            }
        }
        if options.guard_refinement.max_passes == 0 && !residual_policy {
            return Ok(found);
        }
        let mut residuals = BTreeSet::new();
        let mut conditional_points = BTreeSet::new();
        for target in requested {
            run.cancellation.check()?;
            let reduced = found.program.reduce(*target, options.application)?;
            diagnostics.native_rule_applications += reduced.rule_applications;
            retain_conditions(conditions, reduced.nonzero_conditions);
            for integral in reduced.terms.keys() {
                deformation.validate_integral(integral)?;
            }
            for residual in reduced.unresolved {
                deformation.validate_integral(&residual.integral)?;
                if residual_policy
                    && matches!(
                        residual.reason,
                        GuardedApplicationFailure::ConditionVanished { .. }
                    )
                {
                    conditional_points.insert(residual.integral);
                    continue;
                }
                if residual.reason != GuardedApplicationFailure::NoApplicableRule {
                    // Let the ordinary caller report native application failure.
                    return Ok(found);
                }
                // Constant residuals need their own native lower-sector search.
                // They remain provisional until all targets and derivatives are
                // reconstructed by the final source-replayed program.
                if options.search_frontier_sectors
                    || !terminals.is_empty()
                    || !deformation.derivative(residual.integral)?.is_empty()
                {
                    residuals.insert(residual.integral);
                }
            }
        }
        let new_conditional_points = conditional_points
            .difference(&attempted_guard_points)
            .next()
            .is_some();
        let mut added = new_conditional_points
            && pass < options.guard_refinement.max_passes
            && remaining_domains > 0;
        if new_conditional_points && !added {
            diagnostics.guard_refinement_budget_exhausted = true;
        }
        for gap in &found.unresolved {
            if gap.reason != GuardedUnresolvedReason::UnprovedDescent {
                continue;
            }
            let Some(trigger) = residuals.iter().find(|index| gap.domain.contains(index)) else {
                continue;
            };
            if !gap.domain.is_subset_of(deformation.admitted_domain()) {
                return Err(Error::IncompleteReduction(
                    "native guard refinement escaped the admitted index domain".into(),
                ));
            }
            let Some((axis, faces)) = finite_guard_faces(
                &gap.domain,
                options.guard_refinement.max_interval_width,
                options.guard_refinement.max_added_domains,
                learned,
            )?
            else {
                continue;
            };
            let new_faces = faces.iter().filter(|face| !learned.contains(face)).count();
            // A parent subdivision is atomic: never submit a partial partition
            // just because the remaining face budget is smaller than its width.
            if pass == options.guard_refinement.max_passes
                || (residual_policy && remaining_domains == 0)
                || new_faces
                    > options
                        .guard_refinement
                        .max_added_domains
                        .saturating_sub(learned.len())
            {
                diagnostics.guard_refinement_budget_exhausted = true;
                continue;
            }
            diagnostics.guard_refinements.push(GuardRefinementRecord {
                pass: diagnostics.guard_refinement_passes + 1,
                parent_bounds: gap
                    .domain
                    .bounds()
                    .iter()
                    .map(|b| [b.lower(), b.upper()])
                    .collect(),
                native_failure: gap.detail.clone(),
                axis,
                fixed_values: faces
                    .iter()
                    .map(|face| face.bounds()[axis].lower().unwrap())
                    .collect(),
                trigger: trigger.to_vec(),
            });
            for face in faces {
                if !learned.contains(&face) {
                    learned.push(face);
                }
            }
            added = true;
        }
        if !added {
            return Ok(found);
        }
        diagnostics.guard_refinement_passes += 1;
        diagnostics.guard_refinement_added_domains = learned.len();
        if options.max_reused_rules > 0 {
            *retained = Some(found.program);
        }
        if let Some(directory) = options.checkpoints.as_deref() {
            std::fs::create_dir_all(directory)?;
            let metadata = serde_json::json!({
                "schema":1, "options":options.guard_refinement,
                "added_domains":learned.iter().map(|domain| domain.bounds().iter()
                    .map(|bound| [bound.lower(),bound.upper()]).collect::<Vec<_>>()).collect::<Vec<_>>(),
                "native_parent_gaps":diagnostics.guard_refinements,
                "parent_coverage_certified":false,
            });
            std::fs::write(
                directory.join("guard-refinements.json"),
                serde_json::to_vec_pretty(&metadata).map_err(|e| Error::Cache(e.to_string()))?,
            )?;
        }
    }
    unreachable!("last bounded pass cannot add domains")
}

/// Native reduction aggregates each requested integral independently. In
/// particular, coefficients belonging to different targets never cancel here.
#[allow(clippy::too_many_arguments)]
fn residual_search_requests<const N: usize>(
    program: &GuardedReductionProgram<N>,
    requested: &BTreeSet<[i64; N]>,
    deformation: &FixedShellDeformation<N>,
    limits: GuardedReductionLimits,
    diagnostics: &mut WeightedClosureDiagnostics,
    conditions: &mut BTreeMap<String, Atom>,
    run: &RunContext,
) -> Result<(BTreeSet<[i64; N]>, BTreeSet<[i64; N]>, bool)> {
    let mut residuals = BTreeSet::new();
    let mut guard_points = BTreeSet::new();
    let mut failed = false;
    for target in requested {
        run.cancellation.check()?;
        deformation.validate_integral(target)?;
        let reduced = program.reduce(*target, limits)?;
        diagnostics.native_rule_applications += reduced.rule_applications;
        retain_conditions(conditions, reduced.nonzero_conditions);
        for integral in reduced.terms.keys() {
            deformation.validate_integral(integral)?;
        }
        for residual in reduced.unresolved {
            deformation.validate_integral(&residual.integral)?;
            match residual.reason {
                GuardedApplicationFailure::NoApplicableRule => {
                    if !residual.coefficient.is_zero() {
                        residuals.insert(residual.integral);
                    }
                }
                GuardedApplicationFailure::ConditionVanished { .. } => {
                    // Preserve this as a failed conditional application until
                    // another original-source proof actually reduces it. Do
                    // not cancel it with any target or ordinary residual.
                    guard_points.insert(residual.integral);
                }
                _ => failed = true,
            }
        }
    }
    residuals.retain(|point| !guard_points.contains(point));
    Ok((residuals, guard_points, failed))
}

/// Pick one narrowest finite axis whose singleton faces have not all already
/// been submitted. Every other (including unbounded symbolic) axis is retained.
fn finite_guard_faces<const N: usize>(
    parent: &IndexDomain<N>,
    max_width: u64,
    max_faces: usize,
    learned: &[IndexDomain<N>],
) -> Result<Option<(usize, Vec<IndexDomain<N>>)>> {
    let mut axes = parent
        .bounds()
        .iter()
        .enumerate()
        .filter_map(|(axis, bound)| {
            let (lower, upper) = (bound.lower()?, bound.upper()?);
            let width = i128::from(upper) - i128::from(lower);
            (width > 0 && width <= i128::from(max_width) && width < max_faces as i128)
                .then_some((width, axis, lower, upper))
        })
        .collect::<Vec<_>>();
    axes.sort_unstable();
    for (_, axis, lower, upper) in axes {
        let faces = (lower..=upper)
            .map(|value| {
                let mut bounds = *parent.bounds();
                bounds[axis] = IndexBounds::fixed(value);
                IndexDomain::new(bounds).map_err(|e| Error::InvalidInput(e.to_string()))
            })
            .collect::<Result<Vec<_>>>()?;
        if faces.iter().any(|face| !learned.contains(face)) {
            return Ok(Some((axis, faces)));
        }
    }
    Ok(None)
}

fn discovery_domains<const N: usize>(
    requested: &BTreeSet<[i64; N]>,
    roles: &[IndexRole; N],
    split_ordinary_zero_faces: bool,
    requested_index_rays: bool,
) -> Result<Vec<IndexDomain<N>>> {
    let mut domains = Vec::new();
    for indices in requested {
        let mut bounds = [IndexBounds::unbounded(); N];
        for slot in 0..N {
            bounds[slot] = match (roles[slot], indices[slot]) {
                (IndexRole::Occupation, index) if index < 0 => {
                    return Err(Error::InvalidInput(
                        "negative occupation index is undefined".into(),
                    ));
                }
                (_, 0) if requested_index_rays => Ok(IndexBounds::fixed(0)),
                (_, index) if requested_index_rays && index > 0 => {
                    IndexBounds::new(Some(index), None)
                }
                (_, index) if requested_index_rays => IndexBounds::new(None, Some(index)),
                (IndexRole::Occupation, 0) => Ok(IndexBounds::fixed(0)),
                (IndexRole::Ordinary, 0) if split_ordinary_zero_faces => Ok(IndexBounds::fixed(0)),
                (IndexRole::Ordinary, index) if split_ordinary_zero_faces && index < 0 => {
                    IndexBounds::new(None, Some(-1))
                }
                (_, index) if index > 0 => IndexBounds::new(Some(1), None),
                _ => IndexBounds::new(None, Some(0)),
            }
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        }
        let domain = IndexDomain::new(bounds).map_err(|e| Error::InvalidInput(e.to_string()))?;
        if !domains.contains(&domain) {
            domains.push(domain);
        }
    }
    Ok(domains)
}

fn retain_conditions(target: &mut BTreeMap<String, Atom>, conditions: Vec<Atom>) {
    target.extend(
        conditions
            .into_iter()
            .map(|condition| (condition.to_canonical_string(), condition)),
    );
}

fn native_integral<const N: usize>(indices: &[i64; N], physical_arity: usize) -> Result<Integral> {
    validate_storage_label(indices, physical_arity)?;
    Ok(Integral(
        indices[..physical_arity]
            .iter()
            .map(|&index| {
                i16::try_from(index).map_err(|_| {
                    Error::Unsupported("weighted index exceeds RustFlow integral storage".into())
                })
            })
            .collect::<Result<Vec<_>>>()?,
    ))
}

fn linear_combination<const N: usize>(
    terms: &BTreeMap<[i64; N], Atom>,
    physical_arity: usize,
) -> Result<LinearCombination> {
    terms
        .iter()
        .filter_map(|(index, coefficient)| {
            let coefficient = coefficient.together().cancel();
            (!coefficient.is_zero())
                .then(|| Ok((native_integral(index, physical_arity)?, coefficient)))
        })
        .collect()
}

fn unclosed<const N: usize>(
    reason: &str,
    provisional_frontier: Vec<[i64; N]>,
    unresolved: Vec<GuardedAtomUnresolved<N>>,
    discovery_unresolved: Vec<GuardedUnresolved<N>>,
    conditions: BTreeMap<String, Atom>,
    diagnostics: WeightedClosureDiagnostics,
) -> WeightedClosureOutcome<N> {
    WeightedClosureOutcome::Unresolved(UnclosedWeightedSystem {
        reason: reason.into(),
        provisional_frontier,
        unresolved,
        discovery_unresolved,
        nonzero_conditions: conditions.into_values().collect(),
        diagnostics,
    })
}

fn checkpoint<const N: usize>(
    options: &WeightedClosureOptions,
    round: usize,
    closed: bool,
    program: &GuardedReductionProgram<N>,
    requested: &BTreeSet<[i64; N]>,
    frontier: &[[i64; N]],
    active_state: Option<&serde_json::Value>,
    run: &RunContext,
) -> Result<()> {
    let Some(directory) = options.checkpoints.as_deref() else {
        return Ok(());
    };
    std::fs::create_dir_all(directory)?;
    let tag = if closed { "closed" } else { "provisional" };
    let base = format!("round-{round:03}-{tag}");
    let program_path = directory.join(format!("{base}.bin"));
    let temporary = directory.join(format!("{base}.bin.part"));
    let program_bytes = program.encode(Default::default())?;
    run.cancellation.check()?;
    let program_digest = blake3::hash(&program_bytes).to_hex().to_string();
    std::fs::write(&temporary, &program_bytes)?;
    run.cancellation.check()?;
    std::fs::rename(temporary, program_path)?;
    let metadata = serde_json::json!({ "schema": 1, "round": round, "status": tag,
        "measure_id": program.native().sources().measure_id(),
        "physical_arity": program.physical_arity(), "storage_capacity": N,
        "split_ordinary_zero_faces":options.split_ordinary_zero_faces,
        "requested_index_rays":options.requested_index_rays,
        "prioritize_requested_indices":options.prioritize_requested_indices,
        "max_reused_rules":options.max_reused_rules,
        "max_domains_per_residual":options.max_domains_per_residual,
        "requested_ray_point":options.requested_ray_point,
        "discovery_schedule":if options.requested_ray_point.is_some() { "requested-ray-point-v1" } else { "legacy" },
        "program_blake3":program_digest, "program_bytes":program_bytes.len(),
        "discovery":{"max_depth":options.discovery.max_depth,
            "max_domains":options.discovery.max_domains,"sample_seed":options.discovery.sample_seed},
        "application":{"max_rule_applications":options.application.max_rule_applications,
            "max_pending_integrals":options.application.max_pending_integrals},
        "active_state":active_state,
        "max_rounds":options.max_rounds, "max_frontier":options.max_frontier,
        "max_requested":options.max_requested,
        "active_target_closure":options.active_target_closure,
        "max_direct_zero_attempts":options.max_direct_zero_attempts,
        "rule_precedence":if options.max_direct_zero_attempts > 0 {
            "new direct zeros, retained rules, ordinary fresh rules"
        } else if options.max_domains_per_residual > 0 || options.requested_ray_point.is_some() { "retained-first" } else { "fresh-first" },
        "native_rule_count":program.native().rules().len(),
        "requested": requested.iter().map(|indices| indices.to_vec()).collect::<Vec<_>>(),
        "frontier": frontier.iter().map(|indices| indices.to_vec()).collect::<Vec<_>>() });
    let path = directory.join(format!("{base}.json"));
    let temporary = directory.join(format!("{base}.json.part"));
    std::fs::write(
        &temporary,
        serde_json::to_vec_pretty(&metadata).map_err(|e| Error::Cache(e.to_string()))?,
    )?;
    run.cancellation.check()?;
    std::fs::rename(temporary, path)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::guarded::GuardedMeasureIdentity;
    use crate::finite_density::measure::WeightedMeasure;

    #[test]
    fn ordinary_zero_faces_cover_requested_points_without_broadening_admission() {
        let roles = [
            IndexRole::Ordinary,
            IndexRole::Ordinary,
            IndexRole::RequiredCut,
            IndexRole::Occupation,
        ];
        let requested = BTreeSet::from([[0, -2, 1, 0], [-3, 0, 2, 1], [2, 3, 1, 0]]);
        let broad = discovery_domains(&requested, &roles, false, false).unwrap();
        let faces = discovery_domains(&requested, &roles, true, false).unwrap();
        for point in &requested {
            assert!(faces.iter().any(|domain| domain.contains(point)));
        }
        for face in &faces {
            assert!(broad.iter().any(|domain| face.is_subset_of(domain)));
            for a in -4..=4 {
                for b in -4..=4 {
                    for cut in -1..=3 {
                        for occupation in -1..=3 {
                            let point = [a, b, cut, occupation];
                            if face.contains(&point) {
                                assert!(broad.iter().any(|domain| domain.contains(&point)));
                                assert!(occupation >= 0);
                            }
                        }
                    }
                }
            }
        }
        let first = faces
            .iter()
            .find(|domain| domain.contains(&[0, -2, 1, 0]))
            .unwrap();
        assert_eq!(first.bounds()[0], IndexBounds::fixed(0));
        assert_eq!(first.bounds()[1], IndexBounds::new(None, Some(-1)).unwrap());
        assert_eq!(first.bounds()[2], IndexBounds::new(Some(1), None).unwrap());
        assert_eq!(first.bounds()[3], IndexBounds::fixed(0));
        assert!(discovery_domains(&BTreeSet::from([[0, 0, 0, -1]]), &roles, true, false).is_err());
    }

    #[test]
    fn requested_rays_retain_symbolic_indices_and_exact_zero_faces() {
        let roles = [
            IndexRole::Ordinary,
            IndexRole::Ordinary,
            IndexRole::RequiredCut,
            IndexRole::Occupation,
            IndexRole::Occupation,
            IndexRole::Ordinary,
        ];
        let target = [-3, 4, 2, 3, 0, 0];
        let requested = BTreeSet::from([target]);
        let broad = discovery_domains(&requested, &roles, false, false).unwrap();
        let rays = discovery_domains(&requested, &roles, false, true).unwrap();
        assert_eq!(rays.len(), 1);
        assert!(rays[0].contains(&target));
        assert!(rays[0].is_subset_of(&broad[0]));
        assert!(rays[0].contains(&[-300, 400, 200, 300, 0, 0]));
        for invalid in [
            [-2, 4, 2, 3, 0, 0],
            [-3, 3, 2, 3, 0, 0],
            [-3, 4, 1, 3, 0, 0],
            [-3, 4, 2, 2, 0, 0],
            [-3, 4, 2, 3, 1, 0],
            [-3, 4, 2, 3, 0, -1],
        ] {
            assert!(!rays[0].contains(&invalid));
        }
        // Invalid occupation labels must fail before any domain narrowing.
        assert!(
            discovery_domains(&BTreeSet::from([[-3, 4, 2, -1, 0, 0]]), &roles, false, true,)
                .is_err()
        );
    }

    #[test]
    fn generated_compact_sources_close_and_replay_with_discovery_policies() {
        let (context, deformation) = compact_context();
        let target = [1, 0];
        for (
            requested_index_rays,
            prioritize_requested_indices,
            max_reused_rules,
            max_domains_per_residual,
        ) in [
            (false, false, 0, 0),
            (false, true, 0, 0),
            (true, false, 0, 0),
            (true, true, 0, 0),
            (false, false, 4096, 0),
            (false, true, 4096, 0),
            (true, false, 4096, 0),
            (true, true, 4096, 0),
            (true, false, 4096, 32),
            (true, true, 4096, 64),
        ] {
            let outcome = prepare_weighted_system(
                &context,
                &[BTreeMap::from([(target, Atom::one())])],
                &deformation,
                WeightedClosureOptions {
                    split_ordinary_zero_faces: true,
                    requested_index_rays,
                    prioritize_requested_indices,
                    max_reused_rules,
                    max_domains_per_residual,
                    ..Default::default()
                },
                &RunContext::default(),
            )
            .unwrap();
            let WeightedClosureOutcome::Closed(closed) = outcome else {
                panic!("zero-face compact closure unresolved: {outcome:?}")
            };
            closed.differential_system().unwrap().validate().unwrap();
            assert_eq!(
                closed.diagnostics.native_rule_unions > 0,
                max_reused_rules > 0
            );
            if max_reused_rules > 0 {
                assert!(closed.program.native().rules().len() <= max_reused_rules);
            }
            if max_domains_per_residual > 0 {
                assert!(closed.diagnostics.residual_search_requests > 0);
                assert!(closed.diagnostics.residual_allocated_domains > 0);
                assert_eq!(
                    closed.diagnostics.historical_discovery_domains,
                    closed.discovery_unresolved.len()
                );
                assert!(
                    closed.diagnostics.residual_allocated_domains
                        <= (closed.diagnostics.rounds + 1) * 8192
                );
            }
            let decoded = context
                .decode(
                    &closed.program.encode(Default::default()).unwrap(),
                    Default::default(),
                )
                .unwrap();
            assert!(
                decoded
                    .reduce(target, Default::default())
                    .unwrap()
                    .unresolved
                    .is_empty()
            );
            for basis in &closed.reduced.basis {
                let indices = [i64::from(basis.0[0]), i64::from(basis.0[1])];
                for derivative in deformation.derivative(indices).unwrap().keys() {
                    assert!(
                        decoded
                            .reduce(*derivative, Default::default())
                            .unwrap()
                            .unresolved
                            .is_empty()
                    );
                }
            }
        }
    }

    #[test]
    fn reused_native_rules_obey_the_explicit_union_budget() {
        let (context, deformation) = compact_context();
        let result = prepare_weighted_system(
            &context,
            &[BTreeMap::from([([1, 0], Atom::one())])],
            &deformation,
            WeightedClosureOptions {
                max_reused_rules: 1,
                ..Default::default()
            },
            &RunContext::default(),
        );
        assert!(matches!(result,
            Err(Error::Limit(ref message) | Error::Reduction(ref message))
                if message.contains("rule budget")
        ));
    }

    #[test]
    fn residual_policy_rebinds_old_terminals_before_native_discovery() {
        let (context, deformation) = compact_context();
        let target = [2, 0];
        let mut retained = Some(
            context
                .discover(Vec::new(), [target], Default::default())
                .unwrap()
                .program,
        );
        let options = WeightedClosureOptions {
            requested_index_rays: true,
            max_reused_rules: 4096,
            max_domains_per_residual: 32,
            guard_refinement: GuardRefinementOptions {
                max_passes: 0,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let found = discover_with_refinement(
            &context,
            &[],
            &[],
            &BTreeSet::from([target]),
            &deformation,
            &options,
            &mut Vec::new(),
            &mut retained,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut BTreeMap::new(),
            &RunContext::default(),
        )
        .unwrap();
        assert!(found.program.native().terminals().is_empty());
        assert_eq!(diagnostics.residual_search_requests, 1);
        assert_eq!(diagnostics.residual_allocated_domains, 32);
        assert!(
            found
                .program
                .reduce(target, options.application)
                .unwrap()
                .rule_applications
                > 0
        );
    }

    #[test]
    fn residual_policy_preserves_native_application_failure_without_searching() {
        let (context, deformation) = compact_context();
        let outcome = prepare_weighted_system(
            &context,
            &[BTreeMap::from([([1, 0], Atom::one())])],
            &deformation,
            WeightedClosureOptions {
                requested_index_rays: true,
                max_reused_rules: 4096,
                max_domains_per_residual: 32,
                application: GuardedReductionLimits {
                    max_rule_applications: 0,
                    ..Default::default()
                },
                ..Default::default()
            },
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Unresolved(failed) = outcome else {
            panic!("native WorkLimit was not propagated");
        };
        assert_eq!(failed.diagnostics.residual_allocated_domains, 0);
        assert!(
            failed
                .unresolved
                .iter()
                .any(|term| term.reason == GuardedApplicationFailure::WorkLimit)
        );
    }

    fn coupled_guard_context(
        with_point_pivot: bool,
    ) -> (GuardedContext<2>, FixedShellDeformation<2>) {
        use crate::finite_density::guarded::{GuardedIdentity, GuardedIdentityTerm};
        let a = symbol!("conditional_point_a");
        let b = symbol!("conditional_point_b");
        let eta = symbol!("conditional_point_eta");
        let parameter = symbol!("conditional_point_parameter");
        let roles = [IndexRole::Ordinary; 2];
        let domain = IndexDomain::new([
            IndexBounds::new(Some(1), None).unwrap(),
            IndexBounds::new(Some(1), None).unwrap(),
        ])
        .unwrap();
        // A synthetic native source fixture, not a physical integral claim.
        // The first original identity yields a ray rule guarded by a-b. At
        // a=b it cannot solve the target; the second original identity can.
        let identities = [
            ("coupled-pivot", Atom::var(a) - Atom::var(b)),
            ("point-pivot", Atom::var(eta)),
        ]
        .into_iter()
        .take(if with_point_pivot { 2 } else { 1 })
        .map(|(id, pivot)| GuardedIdentity {
            id: id.into(),
            terms: vec![
                GuardedIdentityTerm {
                    shift: [0, 0],
                    coefficient: pivot,
                },
                GuardedIdentityTerm {
                    shift: [-1, 0],
                    coefficient: Atom::num(-1),
                },
            ],
            domain: domain.clone(),
            nonzero_conditions: vec![Atom::var(parameter)],
        })
        .collect();
        let context = GuardedContext::new(
            GuardedMeasureIdentity {
                measure: "synthetic two-index native proof regression".into(),
                support: "original source identities only on positive indices".into(),
                orientation: "formal".into(),
                normalization: "formal".into(),
                branch: "parameter!=0; no physical period supplied".into(),
                deformation: "ordinary index0 derivative".into(),
            },
            roles,
            [a, b],
            vec![eta, parameter],
            identities,
        )
        .unwrap();
        let deformation = FixedShellDeformation::new(
            eta,
            AuxiliaryConvention::NativeMinusEta,
            roles,
            [true, false],
        )
        .unwrap();
        (context, deformation)
    }

    #[test]
    fn conditional_index_face_is_resolved_by_bounded_native_point_replay() {
        let (context, deformation) = coupled_guard_context(true);
        let target = [2, 2];
        let options = WeightedClosureOptions {
            requested_index_rays: true,
            max_reused_rules: 128,
            max_domains_per_residual: 1,
            discovery: GuardedDiscoveryOptions {
                max_domains: 2,
                ..Default::default()
            },
            guard_refinement: GuardRefinementOptions {
                max_passes: 1,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let mut conditions = BTreeMap::new();
        let found = discover_with_refinement(
            &context,
            &[],
            &[],
            &BTreeSet::from([target]),
            &deformation,
            &options,
            &mut Vec::new(),
            &mut None,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(diagnostics.residual_allocated_domains, 2);
        assert_eq!(diagnostics.conditional_point_refinements.len(), 1);
        let point = &diagnostics.conditional_point_refinements[0];
        assert_eq!(point.integral, target);
        assert_eq!(point.pass, 1);
        assert_eq!(point.allocated_domains, 1);
        assert!(found.program.native().terminals().is_empty());
        assert!(found.unresolved.iter().any(|gap| {
            gap.reason == GuardedUnresolvedReason::ExceptionalCondition
                && gap.detail.contains("coupled")
        }));
        // Round-trip replay must still validate both original-source proofs.
        let decoded = context
            .decode(
                &found.program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        let reduced = decoded.reduce(target, options.application).unwrap();
        assert!(reduced.rule_applications > 0);
        assert!(reduced.unresolved.iter().all(|leaf| {
            leaf.reason == GuardedApplicationFailure::NoApplicableRule && leaf.integral != target
        }));
        assert!(
            reduced
                .nonzero_conditions
                .iter()
                .any(|condition| { condition.contains_symbol(symbol!("conditional_point_eta")) })
        );
        assert!(reduced.nonzero_conditions.iter().any(|condition| {
            condition.contains_symbol(symbol!("conditional_point_parameter"))
        }));
    }

    #[test]
    fn conditional_index_face_remains_failure_when_refinement_budget_is_exhausted() {
        let (context, deformation) = coupled_guard_context(true);
        let target = [2, 2];
        for (max_passes, max_domains) in [(0, 2), (1, 1)] {
            let options = WeightedClosureOptions {
                requested_index_rays: true,
                max_reused_rules: 128,
                max_domains_per_residual: 1,
                discovery: GuardedDiscoveryOptions {
                    max_domains,
                    ..Default::default()
                },
                guard_refinement: GuardRefinementOptions {
                    max_passes,
                    ..Default::default()
                },
                ..Default::default()
            };
            let mut diagnostics = WeightedClosureDiagnostics::default();
            let found = discover_with_refinement(
                &context,
                &[],
                &[],
                &BTreeSet::from([target]),
                &deformation,
                &options,
                &mut Vec::new(),
                &mut None,
                &mut BTreeMap::new(),
                &mut diagnostics,
                &mut BTreeMap::new(),
                &RunContext::default(),
            )
            .unwrap();
            assert_eq!(diagnostics.residual_allocated_domains, 1);
            assert!(diagnostics.conditional_point_refinements.is_empty());
            assert!(diagnostics.guard_refinement_budget_exhausted);
            assert!(found.program.native().terminals().is_empty());
            let failed = found.program.reduce(target, options.application).unwrap();
            assert!(failed.unresolved.iter().any(|leaf| {
                leaf.integral == target
                    && matches!(
                        leaf.reason,
                        GuardedApplicationFailure::ConditionVanished { .. }
                    )
            }));
        }
    }

    #[test]
    fn unsuccessful_conditional_point_search_is_not_repeated_or_promoted() {
        let (context, deformation) = coupled_guard_context(false);
        let target = [2, 2];
        let options = WeightedClosureOptions {
            requested_index_rays: true,
            max_reused_rules: 128,
            max_domains_per_residual: 1,
            discovery: GuardedDiscoveryOptions {
                max_domains: 4,
                ..Default::default()
            },
            guard_refinement: GuardRefinementOptions {
                max_passes: 3,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let found = discover_with_refinement(
            &context,
            &[],
            &[],
            &BTreeSet::from([target]),
            &deformation,
            &options,
            &mut Vec::new(),
            &mut None,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut BTreeMap::new(),
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(diagnostics.conditional_point_refinements.len(), 1);
        assert_eq!(diagnostics.residual_allocated_domains, 2);
        assert!(found.program.native().terminals().is_empty());
        assert!(found.unresolved.iter().any(|gap| {
            gap.domain == IndexDomain::new(target.map(IndexBounds::fixed)).unwrap()
        }));
        assert!(
            found
                .program
                .reduce(target, options.application)
                .unwrap()
                .unresolved
                .iter()
                .any(|leaf| matches!(
                    leaf.reason,
                    GuardedApplicationFailure::ConditionVanished { .. }
                ))
        );
    }

    #[test]
    fn residual_policy_keeps_historical_gaps_when_final_search_is_empty() {
        let (context, deformation) = compact_context();
        let requested = BTreeSet::from([[2, 0]]);
        let options = WeightedClosureOptions {
            requested_index_rays: true,
            max_reused_rules: 4096,
            max_domains_per_residual: 1,
            discovery: GuardedDiscoveryOptions {
                max_domains: 1,
                ..Default::default()
            },
            guard_refinement: GuardRefinementOptions {
                max_passes: 0,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut retained = None;
        let mut history = BTreeMap::new();
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let mut conditions = BTreeMap::new();
        let initial = discover_with_refinement(
            &context,
            &[],
            &[],
            &requested,
            &deformation,
            &options,
            &mut Vec::new(),
            &mut retained,
            &mut history,
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert!(
            initial
                .unresolved
                .iter()
                .any(|gap| gap.reason == GuardedUnresolvedReason::DomainBudget)
        );
        let initial_gaps = initial.unresolved.len();
        let allocated = diagnostics.residual_allocated_domains;
        retained = Some(initial.program);
        // Explicit terminals are a low-level stopping set only. This check
        // makes no claim that this chosen singleton is a derivative basis.
        let final_found = discover_with_refinement(
            &context,
            &[],
            &[[2, 0]],
            &requested,
            &deformation,
            &options,
            &mut Vec::new(),
            &mut retained,
            &mut history,
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(diagnostics.residual_allocated_domains, allocated);
        assert_eq!(final_found.unresolved.len(), initial_gaps);
        assert_eq!(diagnostics.historical_discovery_domains, initial_gaps);
    }

    #[test]
    fn guard_faces_cover_only_one_finite_axis_and_preserve_symbolic_bounds() {
        let parent = IndexDomain::new([
            IndexBounds::new(Some(-3), Some(-2)).unwrap(),
            IndexBounds::new(Some(1), None).unwrap(),
            IndexBounds::new(Some(1), Some(3)).unwrap(),
        ])
        .unwrap();
        let (axis, faces) = finite_guard_faces(&parent, 2, 8, &[]).unwrap().unwrap();
        assert_eq!(axis, 0);
        assert_eq!(faces.len(), 2);
        for face in &faces {
            assert!(face.is_subset_of(&parent));
            assert_eq!(face.bounds()[1..], parent.bounds()[1..]);
        }
        for a in -4..=0 {
            for b in 0..=5 {
                for c in 0..=4 {
                    assert_eq!(
                        faces.iter().filter(|d| d.contains(&[a, b, c])).count(),
                        usize::from(parent.contains(&[a, b, c]))
                    );
                }
            }
        }
        assert!(finite_guard_faces(&parent, 0, 8, &[]).unwrap().is_none());
        assert!(finite_guard_faces(&parent, 2, 1, &[]).unwrap().is_none());
        let (axis, second) = finite_guard_faces(&parent, 2, 8, &faces).unwrap().unwrap();
        assert_eq!(axis, 2);
        assert_eq!(second.len(), 3);
        let mut learned = faces;
        learned.extend(second);
        assert!(
            finite_guard_faces(&parent, 2, 8, &learned)
                .unwrap()
                .is_none()
        );
        let extreme =
            IndexDomain::new([IndexBounds::new(Some(i64::MIN), Some(i64::MAX)).unwrap()]).unwrap();
        assert!(
            finite_guard_faces(&extreme, u64::MAX, 8, &[])
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn physical_recenter_gap_refines_through_native_discovery_and_replay() {
        use crate::finite_density::DensityInput;
        use crate::finite_density::preparation::WeightedSourcePolicy;
        let input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        let family = input
            .prepare()
            .unwrap()
            .occupied_cut(&[0, 1], 16)
            .unwrap()
            .at_physical_masses();
        let prepared = family
            .guarded_sources_with_policy::<9>(
                symbol!("guard_refinement_eps"),
                4,
                symbol!("guard_refinement_eta"),
                &[2],
                16,
                vec![],
                GuardedMeasureIdentity {
                    measure: format!("massive sunset double cut; factors={:?}", family.factors()),
                    support: "assigned positive masses; two future shells with 0<E<1".into(),
                    orientation: "future native Minkowski shells".into(),
                    normalization: "unscaled source identities".into(),
                    branch: "generic complex eta and dimension; completion powers<=0".into(),
                    deformation: "D2-eta; occupied factors fixed".into(),
                },
                WeightedSourcePolicy::TangentsThenLorentz,
            )
            .unwrap();
        let parent = IndexDomain::new([
            IndexBounds::fixed(1),
            IndexBounds::fixed(1),
            IndexBounds::fixed(1),
            IndexBounds::new(Some(-3), Some(-2)).unwrap(),
            IndexBounds::fixed(0),
            IndexBounds::fixed(0),
            IndexBounds::fixed(0),
            IndexBounds::new(Some(1), None).unwrap(),
            IndexBounds::fixed(0),
        ])
        .unwrap();
        let target = [1, 1, 1, -2, 0, 0, 0, 1, 0];
        let mut options = WeightedClosureOptions {
            guard_refinement: GuardRefinementOptions {
                max_passes: 0,
                ..Default::default()
            },
            discovery: GuardedDiscoveryOptions {
                max_depth: 3,
                max_domains: 8192,
                sample_seed: 0,
            },
            ..Default::default()
        };
        let mut learned = Vec::new();
        let mut retained = None;
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let mut conditions = BTreeMap::new();
        let requested = BTreeSet::from([target]);
        let initial = discover_with_refinement(
            &prepared.context,
            &[parent.clone()],
            &[],
            &requested,
            &prepared.deformation,
            &options,
            &mut learned,
            &mut retained,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert!(
            initial
                .program
                .reduce(target, options.application)
                .unwrap()
                .unresolved
                .iter()
                .any(|u| u.integral == target)
        );
        options.guard_refinement.max_passes = 1;
        options.guard_refinement.max_added_domains = 2;
        let refined = discover_with_refinement(
            &prepared.context,
            &[parent.clone()],
            &[],
            &requested,
            &prepared.deformation,
            &options,
            &mut learned,
            &mut retained,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(learned.len(), 2);
        assert!(
            learned
                .iter()
                .all(|face| face.bounds()[7] == parent.bounds()[7])
        );
        assert_eq!(diagnostics.guard_refinements[0].axis, 3);
        assert!(refined.unresolved.iter().any(|gap| gap.domain == parent));
        let replayed = prepared
            .context
            .decode(
                &refined.program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        let reduction = replayed.reduce(target, options.application).unwrap();
        assert!(reduction.rule_applications > 0);
        assert!(reduction.unresolved.iter().all(|u| u.integral != target));
        // The new rule is a source-certified partial reduction, not a closed
        // system: its lower-sector residuals remain explicit in this probe.
        assert!(!reduction.unresolved.is_empty());
    }

    pub(super) fn compact_context() -> (GuardedContext<2>, FixedShellDeformation<2>) {
        let z = symbol!("weighted_closure_z");
        let t = symbol!("weighted_closure_t");
        let radius = symbol!("weighted_closure_R");
        let dimension = symbol!("weighted_closure_d");
        let a = symbol!("weighted_closure_a");
        let b = symbol!("weighted_closure_b");
        let roles = [IndexRole::Ordinary, IndexRole::Occupation];
        let measure = WeightedMeasure::new(
            vec![Atom::var(z)],
            [
                Atom::var(z) + Atom::var(t),
                Atom::var(radius) - Atom::var(z),
            ],
            roles,
        )
        .unwrap();
        let mut identities = measure
            .ibp(
                "radial-dilation",
                &[a, b],
                &[Atom::num(2) * Atom::var(z)],
                &Atom::var(dimension),
                2,
            )
            .unwrap();
        identities.extend(measure.multiplication_sources().unwrap());
        let context = GuardedContext::new(
            GuardedMeasureIdentity {
                measure: "radial d-dimensional ball; z=q^2; (z+t)^(-a) H_b(R-z)".into(),
                support: "R>0 and Re(d)>0; 0<=z<=R; uncut t>0".into(),
                orientation: "positive radial Euclidean measure and delta(R-z)".into(),
                normalization: "unit angular prefactor".into(),
                branch: "positive real t, R; meromorphic dimension after exact source construction"
                    .into(),
                deformation: "only uncut z+t; R-z and its occupation derivatives fixed".into(),
            },
            roles,
            [a, b],
            vec![t, radius, dimension],
            identities,
        )
        .unwrap();
        let deformation = FixedShellDeformation::new(
            t,
            AuxiliaryConvention::EuclideanPlusT,
            roles,
            [true, false],
        )
        .unwrap();
        (context, deformation)
    }

    #[test]
    fn fixed_shell_placement_and_auxiliary_sign_are_explicit() {
        let variable = symbol!("weighted_sign_eta");
        let roles = [
            IndexRole::Ordinary,
            IndexRole::RequiredCut,
            IndexRole::Occupation,
        ];
        assert!(
            FixedShellDeformation::new(
                variable,
                AuxiliaryConvention::NativeMinusEta,
                roles,
                [false, true, false]
            )
            .is_err()
        );
        assert!(
            FixedShellDeformation::new(
                variable,
                AuxiliaryConvention::NativeMinusEta,
                roles,
                [false, false, true]
            )
            .is_err()
        );
        assert!(
            FixedShellDeformation::new(
                variable,
                AuxiliaryConvention::NativeMinusEta,
                roles,
                [false, false, false]
            )
            .is_err()
        );
        let native = FixedShellDeformation::new(
            variable,
            AuxiliaryConvention::NativeMinusEta,
            roles,
            [true, false, false],
        )
        .unwrap();
        let euclidean = FixedShellDeformation::new(
            variable,
            AuxiliaryConvention::EuclideanPlusT,
            roles,
            [true, false, false],
        )
        .unwrap();
        assert_eq!(
            native.derivative([2, 1, 0]).unwrap()[&[3, 1, 0]],
            Atom::num(2)
        );
        assert_eq!(
            euclidean.derivative([2, 1, 0]).unwrap()[&[3, 1, 0]],
            Atom::num(-2)
        );
        assert!(native.derivative([0, 1, 0]).unwrap().is_empty());
    }

    #[test]
    fn native_compact_flow_closes_without_a_caller_terminal_basis() {
        let (context, deformation) = compact_context();
        let targets = [
            BTreeMap::from([([1, 0], Atom::num(1))]),
            BTreeMap::from([([1, 0], Atom::i())]),
        ];
        let outcome = prepare_weighted_system(
            &context,
            &targets,
            &deformation,
            WeightedClosureOptions {
                search_frontier_sectors: true,
                discovery: GuardedDiscoveryOptions {
                    max_depth: 3,
                    ..Default::default()
                },
                ..Default::default()
            },
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Closed(closed) = outcome else {
            panic!("compact source closure unresolved: {outcome:?}")
        };
        assert!(!closed.reduced.basis.is_empty());
        closed.differential_system().unwrap().validate().unwrap();
        for (integral, coefficient) in &closed.reduced.targets[0] {
            assert!(
                (&closed.reduced.targets[1][integral] - Atom::i() * coefficient)
                    .together()
                    .cancel()
                    .is_zero()
            );
        }
        assert!(closed.diagnostics.rounds > 1);
        assert!(closed.diagnostics.native_frontier_requests > 0);
        // Every final stopping label was explicitly requested for native
        // discovery, even when its auxiliary derivative vanishes.
        assert!(
            closed
                .reduced
                .basis
                .iter()
                .all(|basis| closed.reduced.candidates.contains_key(basis))
        );
        assert!(!closed.reduced.nonzero_conditions.is_empty());
        // Every basis derivative must now reduce to this exact stopping set;
        // no NoApplicableRule outcome is hidden by the public Closed result.
        for basis in &closed.reduced.basis {
            let indices: [i64; 2] = basis
                .0
                .iter()
                .map(|&index| i64::from(index))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap();
            for derivative in deformation.derivative(indices).unwrap().keys() {
                let reduced = closed
                    .program
                    .reduce(*derivative, Default::default())
                    .unwrap();
                assert!(reduced.unresolved.is_empty());
            }
        }
    }

    #[test]
    fn round_budget_keeps_provisional_frontier_unresolved() {
        let (context, deformation) = compact_context();
        let targets = [BTreeMap::from([([1, 0], Atom::num(1))])];
        let outcome = prepare_weighted_system(
            &context,
            &targets,
            &deformation,
            WeightedClosureOptions {
                max_rounds: 1,
                ..Default::default()
            },
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Unresolved(unclosed) = outcome else {
            panic!("incomplete derivative search unexpectedly closed")
        };
        assert!(!unclosed.provisional_frontier.is_empty());
        assert!(unclosed.reason.contains("round budget"));
    }
}
