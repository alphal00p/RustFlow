//! Bounded, point-specific discovery of one-original-row zero rules.
//!
//! This supplements first-valid recurrence search. It does not combine rows or
//! import an analytic zero: a source term supplies a candidate seed, and the
//! existing instantiator, canonicalizer and independent guarded replay own every
//! accepted equation.

use std::collections::BTreeSet;

use super::{GuardedSearchScope, domain_case, unresolved};
use crate::solver::guarded::{
    GuardedSolution, GuardedSourceSystem, GuardedUnresolvedReason, IndexBounds, IndexDomain,
    IndexRole,
};
use crate::solver::instantiate::canonicalize;
use crate::solver::{
    Integral, IntegralOrder, Power, RuleCandidate, SearchStats, Seed, SeedSource, SolverError,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GuardedDirectZeroSkipReason {
    UnsupportedPower,
    FixedCoordinateMismatch,
    SectorCrossing,
    SourceGuard,
    InvalidOccupationImage,
    RejectedProof,
}

/// A skipped alignment is diagnostic evidence, never a zero certificate.
#[derive(Clone, Debug)]
pub struct GuardedDirectZeroSkip<const N: usize> {
    pub point: [i64; N],
    pub source_ordinal: usize,
    pub term_ordinal: usize,
    pub seed: Option<[i64; N]>,
    pub reason: GuardedDirectZeroSkipReason,
    pub detail: String,
}

#[derive(Debug)]
pub struct GuardedDirectZeroSearch<const N: usize> {
    pub solution: GuardedSolution<N>,
    /// Every source-term alignment consumes one attempt, including rejected
    /// seeds. The budget is shared across all points in this call.
    pub attempted_rows: usize,
    pub skipped_seeds: Vec<GuardedDirectZeroSkip<N>>,
    /// A point is complete after a proved zero or exhaustive alignment search.
    /// Unsuccessful complete points still have explicit SearchExhausted gaps.
    /// Partial/unvisited budget gaps and unsupported points are not complete.
    pub completed_points: Vec<[i64; N]>,
    /// Covered by existing required-cut/explicit source zero domains. No new
    /// rule is fabricated for these points.
    pub already_zero_points: Vec<[i64; N]>,
}

impl<const N: usize> GuardedSourceSystem<N> {
    /// Try source-term-aligned seeds for exact one-row zero rules at concrete
    /// points. All points are role-validated before any work, even beyond the
    /// budget. Duplicate points are visited once, in first-occurrence order.
    ///
    /// Only a nonzero single-term instantiated row can succeed, after the
    /// native scope applies required-cut and uniformly proved measure zeros.
    /// Its canonical
    /// target must equal the requested point and its empty RHS must pass the
    /// existing original-source replay and condition reconstruction. Seeds keep
    /// ordinary/required-cut sectors, while occupation seeds may cross zero
    /// when the source guard permits it. No ordinary zero-sector heuristic is
    /// consulted. Zero attempts returns explicit budget gaps, except for points
    /// already proved zero by the supplied measure.
    pub fn direct_zero_rules_at_points(
        &self,
        points: &[[i64; N]],
        max_attempts: usize,
    ) -> Result<GuardedDirectZeroSearch<N>, SolverError> {
        if points.iter().any(|point| !self.valid_indices(point)) {
            return Err(SolverError::InvalidInput(
                "direct zero point includes a negative occupation index".into(),
            ));
        }
        let mut result = GuardedDirectZeroSearch {
            solution: GuardedSolution {
                rules: Vec::new(),
                unresolved: Vec::new(),
            },
            attempted_rows: 0,
            skipped_seeds: Vec::new(),
            completed_points: Vec::new(),
            already_zero_points: Vec::new(),
        };
        let mut visited = BTreeSet::new();
        for &point in points {
            if !visited.insert(point) {
                continue;
            }
            let domain = IndexDomain::new(point.map(IndexBounds::fixed))?;
            if self.is_zero(&point) {
                result.already_zero_points.push(point);
                result.completed_points.push(point);
                continue;
            }
            if point
                .iter()
                .any(|&n| !(i64::from(Power::MIN)..=i64::from(Power::MAX)).contains(&n))
            {
                unresolved(
                    &mut result.solution,
                    domain,
                    GuardedUnresolvedReason::UnsupportedPower,
                    "direct zero target exceeds native compact power storage",
                );
                continue;
            }
            let (sector, case) = domain_case(&domain)?;
            let target = case.integral();
            let order = IntegralOrder::new(sector, self.roles.map(|r| r == IndexRole::RequiredCut))
                .with_roles(self.roles)?;
            let scope = GuardedSearchScope::new(self, domain.clone());
            let mut proved = false;
            let mut budget_exhausted = false;
            'sources: for (source_ordinal, source) in self.system.rows().iter().enumerate() {
                for (term_ordinal, term) in source.iter().enumerate() {
                    if result.attempted_rows == max_attempts {
                        budget_exhausted = true;
                        break 'sources;
                    }
                    result.attempted_rows += 1;
                    let mut seed_point = point;
                    let mut invalid = None;
                    for axis in 0..N {
                        let power = term.integral[axis];
                        if power.is_symbolic() {
                            let Some(n) = point[axis].checked_sub(i64::from(power.value())) else {
                                invalid = Some((
                                    GuardedDirectZeroSkipReason::UnsupportedPower,
                                    "source-term alignment overflow",
                                ));
                                break;
                            };
                            seed_point[axis] = n;
                        } else if point[axis] != i64::from(power.value()) {
                            invalid = Some((
                                GuardedDirectZeroSkipReason::FixedCoordinateMismatch,
                                "fixed source coordinate does not match the target",
                            ));
                            break;
                        }
                        if !(i64::from(Power::MIN)..=i64::from(Power::MAX))
                            .contains(&seed_point[axis])
                        {
                            invalid = Some((
                                GuardedDirectZeroSkipReason::UnsupportedPower,
                                "aligned seed exceeds native compact power storage",
                            ));
                            break;
                        }
                        if self.roles[axis] != IndexRole::Occupation
                            && (seed_point[axis] > 0) != sector[axis]
                        {
                            invalid = Some((
                                GuardedDirectZeroSkipReason::SectorCrossing,
                                "aligned seed leaves the requested ordinary/cut sector",
                            ));
                            break;
                        }
                    }
                    let skip = |reason, detail: String| GuardedDirectZeroSkip {
                        point,
                        source_ordinal,
                        term_ordinal,
                        seed: Some(seed_point),
                        reason,
                        detail,
                    };
                    if let Some((reason, detail)) = invalid {
                        result.skipped_seeds.push(skip(reason, detail.into()));
                        continue;
                    }
                    let seed = Seed {
                        integral: Integral::numeric(seed_point.map(|n| n as i16))?,
                        shifts: [0; N],
                    };
                    if !scope
                        .source_domain(&seed)?
                        .is_subset_of(&self.sources[source_ordinal].domain)
                    {
                        result.skipped_seeds.push(skip(
                            GuardedDirectZeroSkipReason::SourceGuard,
                            "aligned seed violates the original source guard".into(),
                        ));
                        continue;
                    }
                    let row = match scope.instantiate(source_ordinal, source, &seed, &order) {
                        Ok(Some(row)) => row,
                        Ok(None) => {
                            result.skipped_seeds.push(skip(
                                GuardedDirectZeroSkipReason::InvalidOccupationImage,
                                "instantiated source has an invalid occupation image".into(),
                            ));
                            continue;
                        }
                        Err(SolverError::Power(error)) => {
                            result.skipped_seeds.push(skip(
                                GuardedDirectZeroSkipReason::UnsupportedPower,
                                error.to_string(),
                            ));
                            continue;
                        }
                        Err(error) => return Err(error),
                    };
                    if row.len() != 1 || row[0].integral != target {
                        continue;
                    }
                    let (canonical_target, rhs) = canonicalize(row, self.system.index_variables())?;
                    if canonical_target != target || !rhs.is_empty() {
                        continue;
                    }
                    let candidate = RuleCandidate {
                        case: case.clone().into(),
                        target,
                        rhs,
                        sources: vec![SeedSource {
                            basis_row: source_ordinal,
                            seed,
                        }],
                        stats: SearchStats {
                            seeds: 1,
                            rows: 1,
                            exact_trace_rows: 1,
                            direct_hit: true,
                            ..Default::default()
                        },
                    };
                    match self.seal_candidate(candidate, order.clone(), domain.clone()) {
                        Ok(rule) => {
                            result.solution.rules.push(rule);
                            proved = true;
                            break 'sources;
                        }
                        Err(SolverError::Certification(detail)) => result
                            .skipped_seeds
                            .push(skip(GuardedDirectZeroSkipReason::RejectedProof, detail)),
                        Err(error) => return Err(error),
                    }
                }
            }
            if proved {
                result.completed_points.push(point);
            } else if budget_exhausted {
                unresolved(
                    &mut result.solution,
                    domain,
                    GuardedUnresolvedReason::DomainBudget,
                    "direct zero global source-term attempt budget exhausted",
                );
            } else {
                result.completed_points.push(point);
                unresolved(
                    &mut result.solution,
                    domain,
                    GuardedUnresolvedReason::SearchExhausted,
                    "no one-original-row zero rule among all admissible aligned seeds",
                );
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
