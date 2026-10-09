//! Bounded, native guarded discovery of a derivative-closed weighted basis.
//!
//! Unreduced integrals are provisional candidates only. A successful result
//! independently audits every target reconstruction and every basis derivative
//! with one source-replayed native program. It makes no minimality claim.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

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
            checkpoints: None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct WeightedClosureDiagnostics {
    pub rounds: usize,
    pub requested: usize,
    /// Provisional labels explicitly added for native source discovery, including
    /// auxiliary-constant lower sectors that derivatives cannot request.
    pub native_frontier_requests: usize,
    pub provisional_sizes: Vec<usize>,
    pub native_rules: usize,
    pub native_rule_applications: usize,
    /// Discovery need not solve every native index domain to close the finite
    /// requested target and derivative set. All such gaps stay inspectable.
    pub uncovered_discovery_domains: usize,
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
            &mut diagnostics,
            &mut conditions,
            run,
        )?;
        diagnostics.rounds = round + 1;
        diagnostics.requested = requested.len();
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
            options.checkpoints.as_deref(),
            options.split_ordinary_zero_faces,
            round,
            false,
            &found.program,
            &requested,
            &frontier,
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
                options.checkpoints.as_deref(),
                options.split_ordinary_zero_faces,
                round,
                true,
                &final_found.program,
                &requested,
                &frontier,
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
#[allow(clippy::too_many_arguments)]
fn discover_with_refinement<const N: usize>(
    context: &GuardedContext<N>,
    domains: &[IndexDomain<N>],
    terminals: &[[i64; N]],
    requested: &BTreeSet<[i64; N]>,
    deformation: &FixedShellDeformation<N>,
    options: &WeightedClosureOptions,
    learned: &mut Vec<IndexDomain<N>>,
    diagnostics: &mut WeightedClosureDiagnostics,
    conditions: &mut BTreeMap<String, Atom>,
    run: &RunContext,
) -> Result<GuardedDiscovery<N>> {
    for pass in 0..=options.guard_refinement.max_passes {
        run.cancellation.check()?;
        let mut submitted = domains.to_vec();
        submitted.extend(
            learned
                .iter()
                .filter(|face| !domains.contains(face))
                .cloned(),
        );
        let found = context.discover(submitted, terminals.iter().copied(), options.discovery)?;
        if options.guard_refinement.max_passes == 0 {
            return Ok(found);
        }
        let mut residuals = BTreeSet::new();
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
        let mut added = false;
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
    directory: Option<&Path>,
    split_ordinary_zero_faces: bool,
    round: usize,
    closed: bool,
    program: &GuardedReductionProgram<N>,
    requested: &BTreeSet<[i64; N]>,
    frontier: &[[i64; N]],
) -> Result<()> {
    let Some(directory) = directory else {
        return Ok(());
    };
    std::fs::create_dir_all(directory)?;
    let tag = if closed { "closed" } else { "provisional" };
    let base = format!("round-{round:03}-{tag}");
    let program_path = directory.join(format!("{base}.bin"));
    let temporary = directory.join(format!("{base}.bin.part"));
    std::fs::write(&temporary, program.encode(Default::default())?)?;
    std::fs::rename(temporary, program_path)?;
    let metadata = serde_json::json!({ "schema": 1, "round": round, "status": tag,
        "measure_id": program.native().sources().measure_id(),
        "physical_arity": program.physical_arity(), "storage_capacity": N,
        "split_ordinary_zero_faces":split_ordinary_zero_faces,
        "requested": requested.iter().map(|indices| indices.to_vec()).collect::<Vec<_>>(),
        "frontier": frontier.iter().map(|indices| indices.to_vec()).collect::<Vec<_>>() });
    let path = directory.join(format!("{base}.json"));
    let temporary = directory.join(format!("{base}.json.part"));
    std::fs::write(
        &temporary,
        serde_json::to_vec_pretty(&metadata).map_err(|e| Error::Cache(e.to_string()))?,
    )?;
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
        let broad = discovery_domains(&requested, &roles, false).unwrap();
        let faces = discovery_domains(&requested, &roles, true).unwrap();
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
        assert!(discovery_domains(&BTreeSet::from([[0, 0, 0, -1]]), &roles, true).is_err());
    }

    #[test]
    fn generated_compact_sources_close_and_replay_with_ordinary_zero_faces() {
        let (context, deformation) = compact_context();
        let target = [1, 0];
        let outcome = prepare_weighted_system(
            &context,
            &[BTreeMap::from([(target, Atom::one())])],
            &deformation,
            WeightedClosureOptions {
                split_ordinary_zero_faces: true,
                ..Default::default()
            },
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Closed(closed) = outcome else {
            panic!("zero-face compact closure unresolved: {outcome:?}")
        };
        closed.differential_system().unwrap().validate().unwrap();
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

    fn compact_context() -> (GuardedContext<2>, FixedShellDeformation<2>) {
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
