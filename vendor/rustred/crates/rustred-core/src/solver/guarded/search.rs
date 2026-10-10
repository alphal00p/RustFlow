use super::{
    GuardedRule, GuardedSolution, GuardedSourceSystem, GuardedUnresolved, GuardedUnresolvedReason,
    IndexBounds, IndexDomain, IndexRole,
};
use crate::algebra::Coefficient;
use crate::solver::{
    Case, CoordinateCase, ExactRow, IntegralOrder, PolynomialRow, Power, RuleCandidate,
    SearchOptions, SectorSolver, Seed, SolverError, Term, extract_exceptions,
};
use std::collections::VecDeque;

pub(in crate::solver) struct GuardedSearchScope<'a, const N: usize> {
    pub problem: &'a GuardedSourceSystem<N>,
    pub domain: IndexDomain<N>,
}

impl<'a, const N: usize> GuardedSearchScope<'a, N> {
    pub fn new(problem: &'a GuardedSourceSystem<N>, domain: IndexDomain<N>) -> Self {
        Self { problem, domain }
    }

    pub fn source_domain(&self, seed: &Seed<N>) -> Result<IndexDomain<N>, SolverError> {
        image_domain(&self.domain, &seed.integral)
    }

    pub fn instantiate(
        &self,
        ordinal: usize,
        source: &PolynomialRow<N>,
        seed: &Seed<N>,
        order: &IntegralOrder<N>,
    ) -> Result<Option<ExactRow<N>>, SolverError> {
        let info = self.problem.sources.get(ordinal).ok_or_else(|| {
            SolverError::InvalidInput("guarded source ordinal is out of range".into())
        })?;
        if !self.source_domain(seed)?.is_subset_of(&info.domain) {
            return Ok(None);
        }
        let row = crate::solver::instantiate::instantiate(
            source,
            seed,
            self.problem.system.index_variables(),
            self.problem.system.fixed(),
            order,
            &[],
            None,
        )?;
        let mut result: ExactRow<N> = Vec::with_capacity(row.len());
        for term in row {
            if let Some(last) = result
                .last_mut()
                .filter(|last| last.integral == term.integral)
            {
                last.coefficient = &last.coefficient + &term.coefficient;
            } else {
                result.push(term);
            }
        }
        result.retain(|term| !term.coefficient.is_zero());
        let valid = IndexDomain::for_roles(&self.problem.roles);
        for term in &result {
            if !image_domain(&self.domain, &term.integral)?.is_subset_of(&valid) {
                return Ok(None);
            }
        }
        Ok(Some(result))
    }
}

pub(super) fn image_domain<const N: usize>(
    domain: &IndexDomain<N>,
    integral: &crate::solver::Integral<N>,
) -> Result<IndexDomain<N>, SolverError> {
    let mut bounds = *domain.bounds();
    for (axis, bound) in bounds.iter_mut().enumerate() {
        let power = integral[axis];
        if power.is_symbolic() {
            let shift = i64::from(power.value());
            let add = |value: i64| {
                value.checked_add(shift).ok_or_else(|| {
                    SolverError::InvalidInput("source-domain translation overflow".into())
                })
            };
            *bound = IndexBounds::new(
                bound.lower().map(add).transpose()?,
                bound.upper().map(add).transpose()?,
            )?;
        } else {
            *bound = IndexBounds::fixed(i64::from(power.value()));
        }
    }
    IndexDomain::new(bounds)
}

/// Search hints only: neither queue grants coverage or changes an identity.
struct PendingDomains<'a, const N: usize> {
    points: &'a [[i64; N]],
    prioritized: VecDeque<(IndexDomain<N>, Vec<usize>)>,
    ordinary: VecDeque<IndexDomain<N>>,
}

impl<'a, const N: usize> PendingDomains<'a, N> {
    fn new(domains: Vec<IndexDomain<N>>, points: &'a [[i64; N]]) -> Self {
        let mut pending = Self {
            points,
            prioritized: VecDeque::new(),
            ordinary: VecDeque::new(),
        };
        let all = (0..points.len()).collect::<Vec<_>>();
        pending.extend(domains, &all);
        pending
    }

    // Every child supplied by discovery is a subset of its parent. Filtering
    // the parent's matching point indices therefore retains precisely all
    // global priority points in that child without rescanning unrelated ones.
    fn extend(&mut self, domains: impl IntoIterator<Item = IndexDomain<N>>, parent: &[usize]) {
        for domain in domains {
            let matching = parent
                .iter()
                .copied()
                .filter(|&index| domain.contains(&self.points[index]))
                .collect::<Vec<_>>();
            if matching.is_empty() {
                self.ordinary.push_back(domain);
            } else {
                self.prioritized.push_back((domain, matching));
            }
        }
    }

    fn pop_front(&mut self) -> Option<(IndexDomain<N>, Vec<usize>)> {
        self.prioritized
            .pop_front()
            .or_else(|| self.ordinary.pop_front().map(|domain| (domain, Vec::new())))
    }

    fn into_domains(self) -> impl Iterator<Item = IndexDomain<N>> {
        self.prioritized
            .into_iter()
            .map(|(domain, _)| domain)
            .chain(self.ordinary)
    }
}

impl<const N: usize> GuardedSourceSystem<N> {
    /// Discover reusable rules on one explicit integer box. Search and domain
    /// traversal are bounded; uncovered pieces remain in `unresolved`.
    pub fn solve_domain(
        &self,
        domain: IndexDomain<N>,
        options: SearchOptions,
    ) -> Result<GuardedSolution<N>, SolverError> {
        self.solve_domains(vec![domain], options, 256)
    }

    pub fn solve_domains(
        &self,
        domains: Vec<IndexDomain<N>>,
        options: SearchOptions,
        max_domains: usize,
    ) -> Result<GuardedSolution<N>, SolverError> {
        self.solve_domains_with_priority_points(domains, &[], options, max_domains)
    }

    /// Allocate a separate bounded search to each input box, in input order.
    /// Each allocation is charged in full even if that search visits fewer
    /// domains. Thus the sum of actual visits cannot exceed `max_domains`.
    /// No rule, guard, ordering or proof is changed by this scheduling policy.
    /// All local gaps and every unsubmitted box remain explicit. A zero total
    /// budget is allowed and returns all valid input boxes as budget gaps.
    /// Callers must replay the combined rules with `GuardedProgram::new`.
    pub fn solve_domains_partitioned(
        &self,
        domains: Vec<IndexDomain<N>>,
        priority_points: &[[i64; N]],
        options: SearchOptions,
        max_domains_per_domain: usize,
        max_domains: usize,
    ) -> Result<GuardedSolution<N>, SolverError> {
        if options.max_depth.is_none() || max_domains_per_domain == 0 {
            return Err(SolverError::InvalidInput(
                "partitioned guarded discovery requires finite depth and a positive per-domain allocation".into(),
            ));
        }
        // Validate the complete request, including work beyond the budget.
        let valid = IndexDomain::for_roles(&self.roles);
        if domains.iter().any(|domain| !domain.is_subset_of(&valid)) {
            return Err(SolverError::InvalidInput(
                "requested domain includes a negative occupation index".into(),
            ));
        }
        if priority_points.iter().any(|point| !valid.contains(point)) {
            return Err(SolverError::InvalidInput(
                "priority point includes a negative occupation index".into(),
            ));
        }
        let mut remaining = max_domains;
        let mut solution = GuardedSolution {
            rules: Vec::new(),
            unresolved: Vec::new(),
        };
        for domain in domains {
            let allocation = max_domains_per_domain.min(remaining);
            if allocation == 0 {
                unresolved(
                    &mut solution,
                    domain,
                    GuardedUnresolvedReason::DomainBudget,
                    "partitioned domain allocation budget exhausted",
                );
                continue;
            }
            remaining -= allocation;
            let local = self.solve_domains_with_priority_points(
                vec![domain],
                priority_points,
                options,
                allocation,
            )?;
            solution.rules.extend(local.rules);
            solution.unresolved.extend(local.unresolved);
        }
        Ok(solution)
    }

    /// Discover exactly the requested boxes, visiting pieces containing a
    /// priority point before other pieces. Points are scheduling hints, not
    /// terminal labels or coverage certificates; admissible points outside all
    /// requested boxes have no effect. Negative occupation indices are invalid
    /// hints. Unvisited pieces remain explicit at the budget.
    /// With no points, traversal is identical to `solve_domains`.
    pub fn solve_domains_with_priority_points(
        &self,
        domains: Vec<IndexDomain<N>>,
        priority_points: &[[i64; N]],
        options: SearchOptions,
        max_domains: usize,
    ) -> Result<GuardedSolution<N>, SolverError> {
        if options.max_depth.is_none() || max_domains == 0 {
            return Err(SolverError::InvalidInput(
                "guarded discovery requires a finite search depth and positive domain budget"
                    .into(),
            ));
        }
        let valid = IndexDomain::for_roles(&self.roles);
        if domains.iter().any(|domain| !domain.is_subset_of(&valid)) {
            return Err(SolverError::InvalidInput(
                "requested domain includes a negative occupation index".into(),
            ));
        }
        if priority_points.iter().any(|point| !valid.contains(point)) {
            return Err(SolverError::InvalidInput(
                "priority point includes a negative occupation index".into(),
            ));
        }
        let mut pending = PendingDomains::new(domains, priority_points);
        let mut visited = Vec::new();
        let mut solution = GuardedSolution {
            rules: Vec::new(),
            unresolved: Vec::new(),
        };
        let mut work = 0;
        while let Some((domain, priority_indices)) = pending.pop_front() {
            if visited.contains(&domain) {
                continue;
            }
            if work >= max_domains {
                unresolved(
                    &mut solution,
                    domain,
                    GuardedUnresolvedReason::DomainBudget,
                    "domain traversal budget exhausted",
                );
                for domain in pending.into_domains() {
                    unresolved(
                        &mut solution,
                        domain,
                        GuardedUnresolvedReason::DomainBudget,
                        "domain traversal budget exhausted",
                    );
                }
                break;
            }
            work += 1;
            visited.push(domain.clone());
            // Ordinary signs remain sector coordinates; occupation zero is a
            // fixed bulk face. Neither partition asserts a zero integral.
            if let Some(axis) = domain
                .bounds()
                .iter()
                .position(|b| b.lower().is_none_or(|n| n <= 0) && b.upper().is_none_or(|n| n > 0))
            {
                let (left, right) = domain.split(axis, 0)?;
                pending.extend(left, &priority_indices);
                pending.extend(right, &priority_indices);
                continue;
            }
            if self.roles.iter().enumerate().any(|(axis, role)| {
                *role == IndexRole::RequiredCut
                    && domain.bounds()[axis].upper().is_some_and(|n| n <= 0)
            }) || self
                .zero_domains
                .iter()
                .any(|zero| domain.is_subset_of(zero))
            {
                continue;
            }
            if domain.bounds().iter().any(|bound| {
                bound
                    .lower()
                    .zip(bound.upper())
                    .is_some_and(|(lower, upper)| {
                        lower == upper
                            && !(i64::from(Power::MIN)..=i64::from(Power::MAX)).contains(&lower)
                    })
            }) {
                unresolved(
                    &mut solution,
                    domain,
                    GuardedUnresolvedReason::UnsupportedPower,
                    "fixed target exceeds native compact power storage",
                );
                continue;
            }
            let (sector, case) = domain_case(&domain)?;
            let scope = GuardedSearchScope::new(self, domain.clone());
            let solver = SectorSolver::new_guarded(&self.system, sector, self.roles)?;
            let candidate = match solver.solve_guarded_case(case.clone().into(), options, &scope) {
                Ok(candidate) => candidate,
                Err(SolverError::SearchExhausted { .. }) => {
                    // A boundary-only identity cannot be admitted uniformly on
                    // a ray. Refine at guard endpoints reachable by this seed
                    // budget, then let the native search handle each face.
                    if let Some((axis, at)) =
                        self.source_boundary(&domain, options.max_depth.unwrap())
                    {
                        let (left, right) = domain.split(axis, at)?;
                        pending.extend(left, &priority_indices);
                        pending.extend(right, &priority_indices);
                        continue;
                    }
                    unresolved(
                        &mut solution,
                        domain,
                        GuardedUnresolvedReason::SearchExhausted,
                        "no rule found within the source-seed budget",
                    );
                    continue;
                }
                Err(SolverError::Power(error)) => {
                    unresolved(
                        &mut solution,
                        domain,
                        GuardedUnresolvedReason::UnsupportedPower,
                        &error.to_string(),
                    );
                    continue;
                }
                Err(error) => return Err(error),
            };
            // Make physical sign changes exact before certifying descent.
            // The refined boundary is searched again with its own fixed face.
            if let Some((axis, at)) = candidate.rhs.iter().find_map(|term| {
                (0..N).find_map(|axis| {
                    if self.roles[axis] == IndexRole::Occupation
                        || !term.integral[axis].is_symbolic()
                    {
                        return None;
                    }
                    let at = -i64::from(term.integral[axis].value());
                    let bound = domain.bounds()[axis];
                    (bound.lower().is_none_or(|n| n <= at) && bound.upper().is_none_or(|n| n > at))
                        .then_some((axis, at))
                })
            }) {
                let (left, right) = domain.split(axis, at)?;
                pending.extend(left, &priority_indices);
                pending.extend(right, &priority_indices);
                continue;
            }
            let rule =
                match self.seal_candidate(candidate, solver.ordering().clone(), domain.clone()) {
                    Ok(rule) => rule,
                    Err(SolverError::Certification(detail)) => {
                        unresolved(
                            &mut solution,
                            domain,
                            GuardedUnresolvedReason::UnprovedDescent,
                            &detail,
                        );
                        continue;
                    }
                    Err(error) => return Err(error),
                };
            // A box difference has at most 2N pieces. Preserve them even when
            // the traversal budget is smaller; the queue reports each unvisited
            // piece as unresolved instead of discarding the discovered rule.
            pending.extend(
                domain.difference(&rule.domain, N.saturating_mul(2))?,
                &priority_indices,
            );
            // Reuse native polynomial exception extraction and exact case
            // intersections for index-only poles, including source weights.
            let exceptions = condition_exceptions(&rule, self.system.index_variables(), &sector)?;
            for branch in exceptions.branches {
                match Case::from(case.clone()).intersect_many(
                    &branch,
                    self.system.index_variables(),
                    &sector,
                    Default::default(),
                ) {
                    Ok(intersection) => {
                        for child in intersection.cases {
                            if let Case::Coordinate(child) = child {
                                if let Some(child_domain) =
                                    IndexDomain::from_sector_case(&sector, &child)
                                        .and_then(|d| d.intersection(&rule.domain))
                                {
                                    if child_domain == domain {
                                        unresolved(
                                            &mut solution,
                                            child_domain,
                                            GuardedUnresolvedReason::ExceptionalCondition,
                                            "condition vanishes throughout the requested case",
                                        );
                                    } else {
                                        pending.extend([child_domain], &priority_indices);
                                    }
                                }
                            } else {
                                unresolved(
                                    &mut solution,
                                    rule.domain.clone(),
                                    GuardedUnresolvedReason::ExceptionalCondition,
                                    "coupled exceptional locus remains guarded; coordinate-box refinement is unavailable",
                                );
                            }
                        }
                    }
                    Err(error) => unresolved(
                        &mut solution,
                        rule.domain.clone(),
                        GuardedUnresolvedReason::ExceptionalCondition,
                        &error.to_string(),
                    ),
                }
            }
            solution.rules.push(rule);
        }
        Ok(solution)
    }

    fn source_boundary(&self, domain: &IndexDomain<N>, depth: u32) -> Option<(usize, i64)> {
        // Native compact symbolic displacements cannot exceed this range.
        let depth = i64::from(depth.min(127));
        for source in &self.sources {
            for axis in 0..N {
                let guard = source.domain.bounds()[axis];
                for shift in -depth..=depth {
                    let endpoints = [
                        guard
                            .lower()
                            .and_then(|n| n.checked_sub(shift)?.checked_sub(1)),
                        guard.upper().and_then(|n| n.checked_sub(shift)),
                    ];
                    for at in endpoints.into_iter().flatten() {
                        let bound = domain.bounds()[axis];
                        if bound.lower().is_none_or(|n| n <= at)
                            && bound.upper().is_none_or(|n| n > at)
                        {
                            return Some((axis, at));
                        }
                    }
                }
            }
        }
        None
    }
}

fn domain_case<const N: usize>(
    domain: &IndexDomain<N>,
) -> Result<([bool; N], CoordinateCase<N>), SolverError> {
    let sector = std::array::from_fn(|axis| domain.bounds()[axis].lower().is_some_and(|n| n >= 1));
    let mut fixed = [None; N];
    for (axis, bound) in domain.bounds().iter().enumerate() {
        if let (Some(lower), Some(upper)) = (bound.lower(), bound.upper()) {
            if lower == upper {
                let value = i16::try_from(lower).map_err(|_| {
                    SolverError::InvalidInput("fixed target exceeds compact index storage".into())
                })?;
                Power::new(false, value)?;
                fixed[axis] = Some(value);
            }
        }
    }
    Ok((sector, CoordinateCase::new(fixed)?))
}

fn unresolved<const N: usize>(
    solution: &mut GuardedSolution<N>,
    domain: IndexDomain<N>,
    reason: GuardedUnresolvedReason,
    detail: &str,
) {
    solution.unresolved.push(GuardedUnresolved {
        domain,
        reason,
        detail: detail.into(),
    });
}

fn condition_exceptions<const N: usize>(
    rule: &GuardedRule<N>,
    indices: &[usize; N],
    sector: &[bool; N],
) -> Result<crate::solver::ExceptionalConditions, SolverError> {
    use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;
    let rhs = rule
        .nonzero_conditions
        .iter()
        .map(|p| Term {
            integral: rule.candidate.target,
            coefficient: Coefficient::from_num_den(
                p.one(),
                p.clone(),
                &symbolica::prelude::Z,
                false,
            ),
        })
        .collect();
    let candidate = RuleCandidate {
        case: rule.candidate.case.clone(),
        target: rule.candidate.target,
        rhs,
        sources: Vec::new(),
        stats: Default::default(),
    };
    extract_exceptions(&candidate, indices, sector)
        .map_err(|e| SolverError::ExactReplay(e.to_string()))
}

#[cfg(test)]
mod priority_tests {
    use super::super::{GuardedApplicationStatus, GuardedProgram, GuardedSource};
    use super::*;
    use crate::algebra::CoefficientContext;
    use crate::solver::{Integral, Term};
    use std::sync::Arc;

    fn source() -> (CoefficientContext, Arc<GuardedSourceSystem<1>>) {
        let context = CoefficientContext::new(["priority_n", "priority_x"]);
        let x = context.parameter("priority_x").unwrap();
        let row = GuardedSource::new(
            "positive-ray-recurrence",
            vec![
                Term {
                    integral: Integral::symbolic([0]).unwrap(),
                    coefficient: context.one().numerator,
                },
                Term {
                    integral: Integral::symbolic([-1]).unwrap(),
                    coefficient: (-x.clone()).numerator,
                },
            ],
            IndexDomain::new([IndexBounds::new(Some(2), None).unwrap()]).unwrap(),
        )
        .with_nonzero_conditions(vec![x.numerator]);
        let source = GuardedSourceSystem::new(
            "priority-domain-regression",
            [IndexRole::Ordinary],
            [0],
            vec![row],
        )
        .unwrap();
        (context, Arc::new(source))
    }

    fn options() -> SearchOptions {
        SearchOptions {
            max_depth: Some(2),
            sample_seed: 0,
            ..Default::default()
        }
    }

    #[test]
    fn requested_sign_child_is_searched_before_unrelated_siblings() {
        let (context, source) = source();
        let domain = IndexDomain::unrestricted();
        let old = source
            .solve_domains(vec![domain.clone()], options(), 3)
            .unwrap();
        let old_program = GuardedProgram::new(source.clone(), old.rules, []).unwrap();
        assert!(!matches!(
            old_program.apply(&[2]).unwrap().status,
            GuardedApplicationStatus::Applied { .. }
        ));

        let found = source
            .solve_domains_with_priority_points(vec![domain.clone()], &[[2]], options(), 3)
            .unwrap();
        // Scheduling may postpone regions, but never certifies or discards
        // them. Every sampled point outside the accepted rule is still in an
        // explicitly unresolved box when the three-domain budget stops.
        for point in -8..=8 {
            assert!(
                found
                    .rules
                    .iter()
                    .any(|rule| rule.domain().contains(&[point]))
                    || found
                        .unresolved
                        .iter()
                        .any(|gap| gap.domain.contains(&[point]))
            );
        }
        assert!(
            found
                .unresolved
                .iter()
                .all(|gap| gap.reason == GuardedUnresolvedReason::DomainBudget)
        );
        let program = GuardedProgram::new(source.clone(), found.rules, []).unwrap();
        let applied = program.apply(&[2]).unwrap();
        assert!(matches!(
            applied.status,
            GuardedApplicationStatus::Applied { .. }
        ));
        let x = context.parameter("priority_x").unwrap();
        assert_eq!(applied.terms.get(&[1]), Some(&x));
        assert!(applied.nonzero_conditions.contains(&x.numerator));
        let replay = GuardedProgram::decode_generated(
            &program.encode_native(Default::default()).unwrap(),
            source,
            Default::default(),
        )
        .unwrap();
        assert_eq!(replay.apply(&[2]).unwrap().terms, applied.terms);
        assert_eq!(
            replay.apply(&[2]).unwrap().nonzero_conditions,
            applied.nonzero_conditions
        );
    }

    #[test]
    fn empty_priority_hints_preserve_rule_transport_and_unresolved_order() {
        let (_, source) = source();
        let domains = vec![IndexDomain::unrestricted()];
        let plain = source
            .solve_domains(domains.clone(), options(), 32)
            .unwrap();
        let hinted = source
            .solve_domains_with_priority_points(domains, &[], options(), 32)
            .unwrap();
        assert_eq!(plain.unresolved.len(), hinted.unresolved.len());
        for (a, b) in plain.unresolved.iter().zip(&hinted.unresolved) {
            assert_eq!(a.domain, b.domain);
            assert_eq!(a.reason, b.reason);
            assert_eq!(a.detail, b.detail);
        }
        let plain = GuardedProgram::new(source.clone(), plain.rules, []).unwrap();
        let hinted = GuardedProgram::new(source, hinted.rules, []).unwrap();
        assert!(
            crate::persistence::equivalent_generated_programs(
                &plain.encode_native(Default::default()).unwrap(),
                &hinted.encode_native(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap()
        );
    }

    #[test]
    fn occupation_priority_points_cannot_be_undefined() {
        let program = super::super::lifecycle::tests::sample("priority-role-validation");
        let result = program.sources().solve_domains_with_priority_points(
            vec![IndexDomain::for_roles(&[IndexRole::Occupation])],
            &[[-1]],
            options(),
            8,
        );
        assert!(
            matches!(result, Err(SolverError::InvalidInput(message)) if message.contains("priority point"))
        );
    }

    #[test]
    fn partitioned_search_preserves_rules_conditions_and_all_unvisited_boxes() {
        let (context, source) = source();
        let positive = IndexDomain::new([IndexBounds::new(Some(2), None).unwrap()]).unwrap();
        let unsubmitted = IndexDomain::new([IndexBounds::fixed(9)]).unwrap();
        let found = source
            .solve_domains_partitioned(
                vec![
                    IndexDomain::unrestricted(),
                    positive.clone(),
                    unsubmitted.clone(),
                ],
                &[],
                options(),
                1,
                2,
            )
            .unwrap();
        assert!(
            found.unresolved.iter().any(|gap| gap.domain == unsubmitted
                && gap.reason == GuardedUnresolvedReason::DomainBudget)
        );
        for point in -8..=12 {
            assert!(
                found
                    .rules
                    .iter()
                    .any(|rule| rule.domain().contains(&[point]))
                    || found
                        .unresolved
                        .iter()
                        .any(|gap| gap.domain.contains(&[point]))
            );
        }
        let local = source.solve_domains(vec![positive], options(), 1).unwrap();
        let program = GuardedProgram::new(source.clone(), found.rules, []).unwrap();
        let local_program = GuardedProgram::new(source.clone(), local.rules, []).unwrap();
        let applied = program.apply(&[2]).unwrap();
        assert!(matches!(
            applied.status,
            GuardedApplicationStatus::Applied { .. }
        ));
        assert_eq!(applied.terms, local_program.apply(&[2]).unwrap().terms);
        let x = context.parameter("priority_x").unwrap();
        assert_eq!(applied.terms.get(&[1]), Some(&x));
        assert!(applied.nonzero_conditions.contains(&x.numerator));
        let decoded = GuardedProgram::decode_generated(
            &program.encode_native(Default::default()).unwrap(),
            source,
            Default::default(),
        )
        .unwrap();
        assert_eq!(decoded.apply(&[2]).unwrap().terms, applied.terms);
    }

    #[test]
    fn partitioned_zero_budget_reports_gaps_but_still_validates_all_requests() {
        let program = super::super::lifecycle::tests::sample("partitioned-role-validation");
        let source = program.sources();
        let valid = IndexDomain::for_roles(&[IndexRole::Occupation]);
        let found = source
            .solve_domains_partitioned(vec![valid.clone()], &[], options(), 1, 0)
            .unwrap();
        assert!(found.rules.is_empty());
        assert_eq!(found.unresolved.len(), 1);
        assert_eq!(found.unresolved[0].domain, valid);
        assert_eq!(
            found.unresolved[0].reason,
            GuardedUnresolvedReason::DomainBudget
        );
        let invalid = IndexDomain::new([IndexBounds::fixed(-1)]).unwrap();
        assert!(
            source
                .solve_domains_partitioned(vec![valid.clone(), invalid], &[], options(), 1, 0,)
                .is_err()
        );
        assert!(
            source
                .solve_domains_partitioned(vec![valid.clone()], &[[-1]], options(), 1, 0,)
                .is_err()
        );
        assert!(
            source
                .solve_domains_partitioned(vec![valid], &[], options(), 0, 1)
                .is_err()
        );
    }
}
