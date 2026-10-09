//! Bounded, native guarded discovery of a derivative-closed weighted basis.
//!
//! Unreduced integrals are provisional candidates only. A successful result
//! independently audits every target reconstruction and every basis derivative
//! with one source-replayed native program. It makes no minimality claim.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use symbolica::prelude::*;

use super::guarded::{
    GuardedApplicationFailure, GuardedAtomUnresolved, GuardedContext, GuardedDiscoveryOptions,
    GuardedReductionLimits, GuardedReductionProgram, GuardedUnresolved, IndexBounds, IndexDomain,
    IndexRole,
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
        })
    }

    /// Restrict valid integral labels independently of the native source guards.
    /// Invalid labels must never become provisional or constant terminal terms.
    pub fn with_admitted_domain(mut self, domain: IndexDomain<N>) -> Result<Self> {
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

#[derive(Clone, Debug)]
pub struct WeightedClosureOptions {
    pub max_rounds: usize,
    pub max_frontier: usize,
    pub max_requested: usize,
    pub discovery: GuardedDiscoveryOptions,
    pub application: GuardedReductionLimits,
    /// A distinct directory per source/deformation context. Each round stores
    /// its exact native program plus explicitly provisional/closed metadata.
    pub checkpoints: Option<PathBuf>,
}

impl Default for WeightedClosureOptions {
    fn default() -> Self {
        Self {
            max_rounds: 8,
            max_frontier: 256,
            max_requested: 4096,
            discovery: GuardedDiscoveryOptions::default(),
            application: GuardedReductionLimits::default(),
            checkpoints: None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct WeightedClosureDiagnostics {
    pub rounds: usize,
    pub requested: usize,
    pub provisional_sizes: Vec<usize>,
    pub native_rules: usize,
    pub native_rule_applications: usize,
    /// Discovery need not solve every native index domain to close the finite
    /// requested target and derivative set. All such gaps stay inspectable.
    pub uncovered_discovery_domains: usize,
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
    if context.sources().roles() != &deformation.roles {
        return Err(Error::InvalidInput(
            "weighted source/deformation index roles differ".into(),
        ));
    }
    if !context
        .sources()
        .native_sources()
        .coefficient_variables()
        .contains(&PolyVariable::Symbol(deformation.variable))
    {
        return Err(Error::InvalidInput(
            "auxiliary variable is absent from the weighted source context".into(),
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
        let domains = discovery_domains(&requested, &deformation.roles)?
            .into_iter()
            .filter_map(|domain| domain.intersection(deformation.admitted_domain()))
            .collect::<Vec<_>>();
        let found = context.discover(domains.clone(), [], options.discovery)?;
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
        if new_derivatives.is_empty() {
            // Rebuild one source-replayed program whose explicit stopping set
            // is now audited against every derivative and target. This is the
            // only point at which provisional candidates become a closed basis.
            let final_found =
                context.discover(domains, frontier.iter().copied(), options.discovery)?;
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
                    linear_combination(&weights)
                })
                .collect::<Result<Vec<_>>>()?;
            let candidates = reductions
                .iter()
                .map(|(indices, terms)| Ok((native_integral(indices)?, linear_combination(terms)?)))
                .collect::<Result<BTreeMap<_, _>>>()?;
            checkpoint(
                options.checkpoints.as_deref(),
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
                        .map(native_integral)
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

fn discovery_domains<const N: usize>(
    requested: &BTreeSet<[i64; N]>,
    roles: &[IndexRole; N],
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

fn native_integral<const N: usize>(indices: &[i64; N]) -> Result<Integral> {
    Ok(Integral(
        indices
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
) -> Result<LinearCombination> {
    terms
        .iter()
        .filter_map(|(index, coefficient)| {
            let coefficient = coefficient.together().cancel();
            (!coefficient.is_zero()).then(|| Ok((native_integral(index)?, coefficient)))
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
