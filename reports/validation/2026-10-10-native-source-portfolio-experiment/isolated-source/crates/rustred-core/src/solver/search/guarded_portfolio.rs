//! Isolated point-only source portfolio. This module is absent from production.
//! Selection changes search traversal only; original source replay is authority.
use super::*;
use crate::algebra::CoefficientPolynomial;
use crate::solver::guarded::{GuardedSearchScope, IndexRole};

#[derive(Clone, Copy, Debug, Default)]
pub struct GuardedPortfolioTrialStats {
    pub selector: &'static str,
    pub selected_source_rows: usize,
    pub attempted_rows: usize,
    pub accepted_rows: usize,
    pub seeds: usize,
    pub independent_rows: usize,
    pub exact_materialization: Duration,
    pub guard_rejected_rows: usize,
    pub empty_rows: usize,
    pub exact_trace_rows: usize,
    pub exact_trace_terms: usize,
    pub rhs_terms: Option<usize>,
    pub nonzero_conditions: usize,
    pub completion: &'static str,
    pub search_elapsed: Duration,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct GuardedPortfolioStats {
    pub baseline_rhs_terms: usize,
    pub baseline_rows: usize,
    pub baseline_seeds: usize,
    pub baseline_independent_rows: usize,
    pub baseline_exact_trace_rows: usize,
    pub baseline_direct_hit: bool,
    pub baseline_elapsed: Duration,
    pub baseline_exact_materialization: Duration,
    pub selected_rhs_terms: usize,
    /// Zero is baseline, one strict, two nonincreasing/index-dependent.
    pub selected_arm: usize,
    pub trials: [GuardedPortfolioTrialStats; 2],
    pub preseal_calls: usize,
    pub preseal_elapsed: Duration,
    /// Selection, optional searches and pre-sealing; excludes baseline search.
    pub policy_elapsed: Duration,
    /// baseline_elapsed + policy_elapsed; final outer guarded seal excluded.
    pub total_elapsed: Duration,
    pub completion: &'static str,
}

#[derive(Clone, Copy)]
pub(super) struct PortfolioLimits {
    pub attempted_rows: usize,
    pub native: RuleTrialLimits,
}
pub(super) const EXPERIMENT_LIMITS: PortfolioLimits = PortfolioLimits {
    attempted_rows: 2048,
    native: RuleTrialLimits {
        max_depth: 3,
        max_rows: 512,
        max_exact_trace_rows: 64,
        max_exact_trace_terms: 16384,
    },
};

pub(super) struct GuardedTrialTraversal<'a> {
    pub allowed: &'a [usize],
    pub max_attempted: usize,
    pub attempted: usize,
    pub guard_rejected: usize,
    pub empty_rows: usize,
    pub attempt_cap_hit: bool,
}

pub(super) fn row_selected<const N: usize>(
    row: &PolynomialRow<N>,
    roles: &[IndexRole; N],
    indices: &[usize; N],
    strict: bool,
) -> bool {
    let role_ok = !row.is_empty()
        && row.iter().all(|term| {
            roles.iter().enumerate().all(|(axis, role)| {
                let power = term.integral[axis];
                match role {
                    IndexRole::Ordinary => true,
                    IndexRole::Occupation => power.is_symbolic() && power.value() == 0,
                    IndexRole::RequiredCut => {
                        power.is_symbolic()
                            && if strict {
                                power.value() == 0
                            } else {
                                power.value() <= 0
                            }
                    }
                }
            })
        });
    let coefficient_ok = strict
        || row.iter().any(|term| {
            indices
                .iter()
                .any(|&index| term.coefficient.degree(index) != 0)
        });
    role_ok && coefficient_ok
}
pub(super) fn selected_ordinals<const N: usize>(
    sources: &crate::solver::guarded::GuardedSourceSystem<N>,
    strict: bool,
) -> Vec<usize> {
    sources
        .system
        .rows()
        .iter()
        .enumerate()
        .filter_map(|(ordinal, row)| {
            row_selected(
                row,
                &sources.roles,
                sources.system.index_variables(),
                strict,
            )
            .then_some(ordinal)
        })
        .collect()
}

/// A conservative implication: baseline nonzero conditions imply alternative
/// conditions only through nonzero constants or exact rational associates.
/// No factoring, division by a parameter or inference from sampled values.
pub(super) fn conditions_covered(
    alternative: &[CoefficientPolynomial],
    baseline: &[CoefficientPolynomial],
    variables: &[symbolica::poly::PolyVariable],
) -> bool {
    alternative.iter().all(|condition| {
        if condition.is_zero() || condition.variables().as_ref() != variables {
            return false;
        }
        if condition
            .exponents_iter()
            .all(|powers| powers.iter().all(|&p| p == 0))
        {
            return true;
        }
        baseline.iter().any(|original| {
            if original.is_zero() || original.variables() != condition.variables() {
                return false;
            }
            // Both native polynomial coefficients are exact integers. Cross
            // multiplication by any nonzero coefficient tests rational
            // proportionality without discarding or altering a condition.
            let c = condition.coefficients[0].clone();
            let b = original.coefficients[0].clone();
            condition.clone().mul_coeff(b) == original.clone().mul_coeff(c)
        })
    })
}

impl<const N: usize> SectorSolver<'_, N> {
    pub(super) fn select_guarded_portfolio(
        &self,
        case: Case<N>,
        options: SearchOptions,
        scope: &GuardedSearchScope<'_, N>,
        baseline: RuleCandidate<N>,
        limits: PortfolioLimits,
    ) -> Result<RuleCandidate<N>, SolverError> {
        // Full baseline ran first. Rays, affine cases and ordinary calls keep
        // their unchanged owner. Do not widen this isolated point pilot.
        if case.coordinate().is_none()
            || !case.is_numerical()
            || !scope
                .domain
                .bounds()
                .iter()
                .all(|b| b.lower().is_some() && b.lower() == b.upper())
            || baseline.rhs.is_empty()
        {
            return Ok(baseline);
        }
        let start = Instant::now();
        let mut stats = GuardedPortfolioStats {
            baseline_rhs_terms: baseline.rhs.len(),
            baseline_rows: baseline.stats.rows,
            baseline_seeds: baseline.stats.seeds,
            baseline_independent_rows: baseline.stats.independent_rows,
            baseline_exact_trace_rows: baseline.stats.exact_trace_rows,
            baseline_direct_hit: baseline.stats.direct_hit,
            baseline_elapsed: baseline.stats.elapsed,
            baseline_exact_materialization: baseline.stats.exact_materialization,
            selected_rhs_terms: baseline.rhs.len(),
            completion: "baseline-retained",
            ..Default::default()
        };
        let seal_start = Instant::now();
        stats.preseal_calls += 1;
        let original = scope.problem.seal_candidate(
            baseline.clone(),
            self.order.clone(),
            scope.domain.clone(),
        );
        stats.preseal_elapsed += seal_start.elapsed();
        let original = match original {
            Ok(rule) if scope.domain.is_subset_of(rule.domain()) => rule,
            Ok(_) | Err(SolverError::Certification(_)) => {
                stats.completion = "baseline-not-presealed-original-path";
                stats.policy_elapsed = start.elapsed();
                stats.total_elapsed = stats.baseline_elapsed + stats.policy_elapsed;
                let mut result = baseline;
                result.stats.guarded_portfolio = Some(stats);
                return Ok(result);
            }
            Err(error) => return Err(error),
        };
        let baseline_conditions = original.nonzero_conditions();
        let mut best = baseline;
        for arm in 0..2 {
            if best.rhs.is_empty() {
                break;
            }
            let allowed = selected_ordinals(scope.problem, arm == 0);
            let mut trial_stats = GuardedPortfolioTrialStats {
                selector: if arm == 0 {
                    "strict-protected-symbolic-zero"
                } else {
                    "nonincreasing-cut-index-dependent"
                },
                selected_source_rows: allowed.len(),
                ..Default::default()
            };
            if allowed.is_empty() {
                trial_stats.completion = "empty-selector";
                stats.trials[arm] = trial_stats;
                continue;
            }
            let mut traversal = GuardedTrialTraversal {
                allowed: &allowed,
                max_attempted: limits.attempted_rows,
                attempted: 0,
                guard_rejected: 0,
                empty_rows: 0,
                attempt_cap_hit: false,
            };
            let mut work = RuleTrialStats::default();
            let trial_start = Instant::now();
            let mut trial_options = options;
            trial_options.max_depth = Some(
                options
                    .max_depth
                    .map_or(limits.native.max_depth, |d| d.min(limits.native.max_depth)),
            );
            // account_trace=false keeps native rows as accepted instantiations;
            // traversal counts attempted rows before guard rejection separately.
            let result = self.search_attempt_inner(
                case.clone(),
                trial_options,
                None,
                Some(limits.native),
                false,
                &mut work,
                |_| {},
                Some(scope),
                Some(&mut traversal),
            );
            work.search.elapsed = trial_start.elapsed();
            trial_stats.search_elapsed = work.search.elapsed;
            trial_stats.attempted_rows = traversal.attempted;
            trial_stats.accepted_rows = work.search.rows;
            trial_stats.seeds = work.search.seeds;
            trial_stats.independent_rows = work.search.independent_rows;
            trial_stats.exact_materialization = work.search.exact_materialization;
            trial_stats.guard_rejected_rows = traversal.guard_rejected;
            trial_stats.empty_rows = traversal.empty_rows;
            trial_stats.exact_trace_rows = work.search.exact_trace_rows;
            trial_stats.exact_trace_terms = work.exact_trace_terms;
            let mut candidate = match result {
                Ok(candidate) => candidate,
                Err(TrialSearchError::Limit(kind)) => {
                    trial_stats.completion = match kind {
                        RuleTrialBudget::SourceRows if traversal.attempt_cap_hit => {
                            "attempted-row-budget"
                        }
                        RuleTrialBudget::SourceRows => "accepted-row-budget",
                        RuleTrialBudget::ExactTraceRows => "exact-trace-row-budget",
                        RuleTrialBudget::ExactTraceTerms => "exact-trace-term-budget",
                    };
                    stats.trials[arm] = trial_stats;
                    continue;
                }
                Err(TrialSearchError::Solver(SolverError::SearchExhausted { .. })) => {
                    trial_stats.completion = "search-exhausted";
                    stats.trials[arm] = trial_stats;
                    continue;
                }
                Err(TrialSearchError::Solver(error)) => return Err(error),
            };
            candidate.stats = work.search;
            trial_stats.rhs_terms = Some(candidate.rhs.len());
            // Replay even a longer candidate: diagnostics distinguish exact
            // admission/condition failure from merely a worse search result.
            stats.preseal_calls += 1;
            let seal_start = Instant::now();
            let admitted = scope.problem.seal_candidate(
                candidate.clone(),
                self.order.clone(),
                scope.domain.clone(),
            );
            stats.preseal_elapsed += seal_start.elapsed();
            let admitted = match admitted {
                Ok(rule) if scope.domain.is_subset_of(rule.domain()) => rule,
                Ok(_) => {
                    trial_stats.completion = "rejected-domain-shrink";
                    stats.trials[arm] = trial_stats;
                    continue;
                }
                Err(SolverError::Certification(_)) => {
                    trial_stats.completion = "rejected-certification";
                    stats.trials[arm] = trial_stats;
                    continue;
                }
                Err(error) => return Err(error),
            };
            trial_stats.nonzero_conditions = admitted.nonzero_conditions().len();
            if !conditions_covered(
                admitted.nonzero_conditions(),
                baseline_conditions,
                self.system.coefficient_variables(),
            ) {
                trial_stats.completion = "rejected-new-condition-locus";
            } else if candidate.rhs.len() < best.rhs.len() {
                stats.selected_arm = arm + 1;
                stats.selected_rhs_terms = candidate.rhs.len();
                stats.completion = "strictly-shorter-covered-alternative";
                trial_stats.completion = "selected-strictly-shorter";
                best = candidate;
            } else {
                trial_stats.completion = "not-strictly-shorter";
            }
            stats.trials[arm] = trial_stats;
        }
        stats.policy_elapsed = start.elapsed();
        stats.total_elapsed = stats.baseline_elapsed + stats.policy_elapsed;
        best.stats.guarded_portfolio = Some(stats);
        Ok(best)
    }
}
