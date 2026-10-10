//! Compiled only under explicit isolated cfg rustflow_portfolio_controls.
//! External runner calls these actual native branch controls without compiling
//! the unrelated monolithic native test suite. No physical values are supplied.
use super::guarded_portfolio::*;
use super::*;
use crate::algebra::{CoefficientContext, CoefficientPolynomial};
use crate::solver::guarded::{
    GuardedApplicationStatus, GuardedProgram, GuardedSearchScope, GuardedSource,
    GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole,
};
use crate::solver::{CoordinateCase, Power};

fn context() -> CoefficientContext {
    CoefficientContext::try_new(["portfolio_n", "portfolio_x", "portfolio_y"]).unwrap()
}
fn point(n: i64) -> IndexDomain<1> {
    IndexDomain::new([IndexBounds::fixed(n)]).unwrap()
}
fn row(c: &CoefficientPolynomial, shifts: &[i16]) -> PolynomialRow<1> {
    shifts
        .iter()
        .map(|&shift| Term {
            integral: Integral::symbolic([shift]).unwrap(),
            coefficient: c.clone(),
        })
        .collect()
}
fn make(
    rows: Vec<(PolynomialRow<1>, IndexDomain<1>, Vec<CoefficientPolynomial>)>,
    role: IndexRole,
) -> Arc<GuardedSourceSystem<1>> {
    let sources = rows
        .into_iter()
        .enumerate()
        .map(|(i, (r, d, c))| {
            GuardedSource::new(format!("original-{i}"), r, d).with_nonzero_conditions(c)
        })
        .collect();
    Arc::new(
        GuardedSourceSystem::new("isolated-point-portfolio-controls", [role], [0], sources)
            .unwrap(),
    )
}
fn fixture(
    first_indirect: bool,
    trial_condition: Option<CoefficientPolynomial>,
) -> Arc<GuardedSourceSystem<1>> {
    let c = context();
    let one = c.one().numerator;
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let x = c.parameter("portfolio_x").unwrap().numerator;
    let mut rows = if first_indirect {
        vec![
            (row(&one, &[1, 0]), point(1), vec![x.clone()]),
            (row(&one, &[1, -1]), point(1), vec![x.clone()]),
        ]
    } else {
        vec![(row(&one, &[0, -1]), point(1), vec![x.clone()])]
    };
    let conditions = trial_condition.map(|p| vec![p]).unwrap_or_else(|| vec![x]);
    rows.push((row(&n, &[2, 0]), point(1), conditions.clone()));
    rows.push((row(&n, &[2]), point(1), conditions));
    make(rows, IndexRole::Ordinary)
}
fn raw<const N: usize>(
    solver: &SectorSolver<'_, N>,
    case: Case<N>,
    options: SearchOptions,
    scope: &GuardedSearchScope<'_, N>,
) -> Result<RuleCandidate<N>, SolverError> {
    solver.validate_search(&case, options)?;
    let mut work = RuleTrialStats::default();
    match solver.search_attempt_inner(
        case,
        options,
        None,
        None,
        false,
        &mut work,
        |_| {},
        Some(scope),
        None,
    ) {
        Ok(mut c) => {
            c.stats = work.search;
            Ok(c)
        }
        Err(TrialSearchError::Solver(e)) => Err(e),
        Err(TrialSearchError::Limit(_)) => unreachable!(),
    }
}
fn run(
    system: &GuardedSourceSystem<1>,
    target: i16,
    limits: PortfolioLimits,
    prime: u64,
) -> Result<RuleCandidate<1>, SolverError> {
    let solver = SectorSolver::new_guarded(&system.system, [target > 0], system.roles).unwrap();
    let scope = GuardedSearchScope::new(system, point(i64::from(target)));
    let case: Case<1> = CoordinateCase::new([Some(target)]).unwrap().into();
    let options = SearchOptions {
        max_depth: Some(0),
        prime,
        ..Default::default()
    };
    let baseline = raw(&solver, case.clone(), options, &scope)?;
    solver.select_guarded_portfolio(case, options, &scope, baseline, limits)
}
fn standard(system: &GuardedSourceSystem<1>) -> RuleCandidate<1> {
    run(system, 1, EXPERIMENT_LIMITS, SearchOptions::default().prime).unwrap()
}
fn publish(system: Arc<GuardedSourceSystem<1>>, candidate: RuleCandidate<1>) {
    let solver = SectorSolver::new_guarded(&system.system, [true], system.roles).unwrap();
    let rule = system
        .seal_candidate(candidate, solver.ordering().clone(), point(1))
        .unwrap();
    let p = GuardedProgram::new(system.clone(), vec![rule], []).unwrap();
    let a = p.apply(&[1]).unwrap();
    assert!(matches!(a.status, GuardedApplicationStatus::Applied { .. }));
    let encoded = p.encode_native(Default::default()).unwrap();
    let q = GuardedProgram::decode_generated(&encoded, system, Default::default()).unwrap();
    let b = q.apply(&[1]).unwrap();
    assert_eq!(a.terms, b.terms);
    assert_eq!(a.nonzero_conditions, b.nonzero_conditions);
    q.sources().replay_rule(&q.rules()[0]).unwrap();
}
fn condition_gate_exact_and_map_safe() {
    let c = context();
    let x = c.parameter("portfolio_x").unwrap().numerator;
    let y = c.parameter("portfolio_y").unwrap().numerator;
    let vars = x.variables().clone();
    let six = x.clone().mul_coeff(6.into());
    let minus = x.clone().mul_coeff((-3).into());
    assert!(conditions_covered(&[six], &[minus], vars.as_ref()));
    assert!(conditions_covered(&[], &[x.clone()], vars.as_ref()));
    assert!(conditions_covered(
        &[c.integer(-7).numerator],
        &[],
        vars.as_ref()
    ));
    assert!(!conditions_covered(
        &[c.zero().numerator],
        &[x.clone()],
        vars.as_ref()
    ));
    assert!(!conditions_covered(&[&x * &x], &[x.clone()], vars.as_ref()));
    assert!(!conditions_covered(
        &[&x * &y],
        &[x.clone(), y.clone()],
        vars.as_ref()
    ));
    assert!(!conditions_covered(&[y], &[x.clone()], vars.as_ref()));
    let foreign = CoefficientContext::try_new(["portfolio_foreign"])
        .unwrap()
        .one()
        .numerator;
    assert!(!conditions_covered(&[foreign], &[x], vars.as_ref()));
}
fn selectors_protect_all_roles_and_literal_indices() {
    let c = CoefficientContext::try_new(["ps0", "ps1", "ps2", "ps3", "ps4"]).unwrap();
    let n = c.parameter("ps4").unwrap().numerator;
    let roles = [
        IndexRole::RequiredCut,
        IndexRole::Occupation,
        IndexRole::RequiredCut,
        IndexRole::Occupation,
        IndexRole::Ordinary,
    ];
    let term = Term {
        integral: Integral::symbolic([0, 0, 0, 0, 1]).unwrap(),
        coefficient: n.clone(),
    };
    assert!(row_selected(
        &vec![term.clone()],
        &roles,
        &[0, 1, 2, 3, 4],
        true
    ));
    assert!(row_selected(
        &vec![term.clone()],
        &roles,
        &[0, 1, 2, 3, 4],
        false
    ));
    for axis in 0..4 {
        let mut t = term.clone();
        let mut powers = *t.integral.powers();
        powers[axis] = Power::new(false, 0).unwrap();
        t.integral = Integral::new(powers);
        assert!(!row_selected(&vec![t], &roles, &[0, 1, 2, 3, 4], true));
        let mut t = term.clone();
        let mut powers = *t.integral.powers();
        powers[axis] = Power::new(true, 1).unwrap();
        t.integral = Integral::new(powers);
        assert!(!row_selected(&vec![t], &roles, &[0, 1, 2, 3, 4], false));
    }
    let mut lower = term.clone();
    let mut powers = *lower.integral.powers();
    powers[2] = Power::new(true, -1).unwrap();
    lower.integral = Integral::new(powers);
    assert!(!row_selected(
        &vec![lower.clone()],
        &roles,
        &[0, 1, 2, 3, 4],
        true
    ));
    assert!(row_selected(&vec![lower], &roles, &[0, 1, 2, 3, 4], false));
    let constant = Term {
        coefficient: c.one().numerator,
        ..term
    };
    assert!(!row_selected(
        &vec![constant],
        &roles,
        &[0, 1, 2, 3, 4],
        false
    ));
    assert!(!row_selected::<5>(&vec![], &roles, &[0, 1, 2, 3, 4], true));
}
fn first_direct_shorter_zero_keeps_full_ordinals_and_roundtrips() {
    let s = fixture(false, None);
    let c = standard(&s);
    assert!(c.rhs.is_empty());
    let mut ordinals = c.sources.iter().map(|s| s.basis_row).collect::<Vec<_>>();
    ordinals.sort_unstable();
    assert_eq!(ordinals, vec![1, 2]);
    let d = c.stats.guarded_portfolio.unwrap();
    assert_eq!(d.baseline_rhs_terms, 1);
    assert_eq!(d.selected_arm, 2);
    assert_eq!(d.preseal_calls, 3);
    assert_eq!(s.sources().len(), 3);
    publish(s, c);
}
fn first_indirect_is_eligible_and_keeps_full_ordinals() {
    let s = fixture(true, None);
    let c = standard(&s);
    assert!(c.rhs.is_empty());
    let mut ordinals = c.sources.iter().map(|s| s.basis_row).collect::<Vec<_>>();
    ordinals.sort_unstable();
    assert_eq!(ordinals, vec![2, 3]);
    assert_eq!(c.stats.guarded_portfolio.unwrap().selected_arm, 2);
    publish(s, c);
}
fn rational_associate_conditions_accepted_unchanged() {
    let c = context();
    let twice = c
        .parameter("portfolio_x")
        .unwrap()
        .numerator
        .mul_coeff(2.into());
    let s = fixture(false, Some(twice.clone()));
    let c = standard(&s);
    assert!(c.rhs.is_empty());
    let solver = SectorSolver::new_guarded(&s.system, [true], s.roles).unwrap();
    let rule = s
        .seal_candidate(c, solver.ordering().clone(), point(1))
        .unwrap();
    assert!(rule.nonzero_conditions().contains(&twice));
}
fn new_parameter_condition_retains_baseline() {
    let c = context();
    let s = fixture(false, Some(c.parameter("portfolio_y").unwrap().numerator));
    let c = standard(&s);
    assert_eq!(c.rhs.len(), 1);
    let d = c.stats.guarded_portfolio.unwrap();
    assert_eq!(d.selected_arm, 0);
    assert_eq!(d.trials[1].completion, "rejected-new-condition-locus");
    publish(s, c);
}
fn certification_failure_retains_baseline() {
    let c = context();
    let bad = &c.parameter("portfolio_n").unwrap().numerator - &c.one().numerator;
    let s = fixture(false, Some(bad));
    let c = standard(&s);
    assert_eq!(c.rhs.len(), 1);
    assert_eq!(
        c.stats.guarded_portfolio.unwrap().trials[1].completion,
        "rejected-certification"
    );
    publish(s, c);
}
fn baseline_certification_failure_returns_unchanged_candidate() {
    let c = context();
    let bad = &c.parameter("portfolio_n").unwrap().numerator - &c.one().numerator;
    let s = make(
        vec![(row(&c.one().numerator, &[0, -1]), point(1), vec![bad])],
        IndexRole::Ordinary,
    );
    let c = standard(&s);
    assert_eq!(c.rhs.len(), 1);
    assert_eq!(
        c.stats.guarded_portfolio.unwrap().completion,
        "baseline-not-presealed-original-path"
    );
    let solver = SectorSolver::new_guarded(&s.system, [true], s.roles).unwrap();
    assert!(matches!(
        s.seal_candidate(c, solver.ordering().clone(), point(1)),
        Err(SolverError::Certification(_))
    ));
}
fn attempted_guard_rejections_and_accepted_budget_are_distinct() {
    let c = context();
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let s = make(
        vec![
            (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
            (row(&n, &[2, 0]), point(2), vec![]),
            (row(&n, &[2]), point(2), vec![]),
        ],
        IndexRole::Ordinary,
    );
    let mut lim = EXPERIMENT_LIMITS;
    lim.attempted_rows = 2;
    let c = run(&s, 1, lim, SearchOptions::default().prime).unwrap();
    let d = c.stats.guarded_portfolio.unwrap();
    // Depth zero exhausts immediately after this seed, so use depth-independent
    // evidence: both rejected attempts are charged and no row is accepted.
    assert_eq!(d.trials[1].attempted_rows, 2);
    assert_eq!(d.trials[1].guard_rejected_rows, 2);
    assert_eq!(d.trials[1].accepted_rows, 0);
    let s = fixture(false, None);
    let mut lim = EXPERIMENT_LIMITS;
    lim.native.max_rows = 1;
    let c = run(&s, 1, lim, SearchOptions::default().prime).unwrap();
    assert_eq!(c.rhs.len(), 1);
    let t = c.stats.guarded_portfolio.unwrap().trials[1];
    assert_eq!(t.completion, "accepted-row-budget");
    assert_eq!(t.accepted_rows, 1);
}
fn attempt_cap_stops_before_next_selected_instantiation() {
    let c = context();
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let s = make(
        vec![
            (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
            (row(&n, &[2, 0]), point(2), vec![]),
            (row(&n, &[2]), point(2), vec![]),
            (row(&n, &[0]), point(1), vec![]),
        ],
        IndexRole::Ordinary,
    );
    let mut lim = EXPERIMENT_LIMITS;
    lim.attempted_rows = 2;
    let c = run(&s, 1, lim, SearchOptions::default().prime).unwrap();
    assert_eq!(c.rhs.len(), 1);
    let t = c.stats.guarded_portfolio.unwrap().trials[1];
    assert_eq!(t.completion, "attempted-row-budget");
    assert_eq!(t.attempted_rows, 2);
    assert_eq!(t.guard_rejected_rows, 2);
}
fn trace_caps_return_baseline_without_truncated_proof() {
    for rows in [true, false] {
        let s = fixture(false, None);
        let mut lim = EXPERIMENT_LIMITS;
        if rows {
            lim.native.max_exact_trace_rows = 1;
        } else {
            lim.native.max_exact_trace_terms = 2;
        }
        let c = run(&s, 1, lim, SearchOptions::default().prime).unwrap();
        assert_eq!(c.rhs.len(), 1);
        assert_eq!(
            c.stats.guarded_portfolio.unwrap().trials[1].completion,
            if rows {
                "exact-trace-row-budget"
            } else {
                "exact-trace-term-budget"
            }
        );
        publish(s, c);
    }
}
fn optional_power_failure_propagates() {
    let c = context();
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let s = make(
        vec![
            (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
            (row(&n, &[Power::MAX]), point(1), vec![]),
        ],
        IndexRole::Ordinary,
    );
    assert!(matches!(
        run(&s, 1, EXPERIMENT_LIMITS, SearchOptions::default().prime),
        Err(SolverError::Power(_))
    ));
}
fn optional_exact_lift_failure_propagates() {
    let c = context();
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let x = c.parameter("portfolio_x").unwrap().numerator;
    let q = &(&x - &c.one().numerator) * &(&x - &c.integer(2).numerator);
    let r = vec![
        Term {
            integral: Integral::symbolic([2]).unwrap(),
            coefficient: &q * &n,
        },
        Term {
            integral: Integral::symbolic([1]).unwrap(),
            coefficient: n.clone(),
        },
        Term {
            integral: Integral::symbolic([0]).unwrap(),
            coefficient: n.clone(),
        },
    ];
    let s = make(
        vec![
            (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
            (r, point(1), vec![]),
            (row(&n, &[1]), point(1), vec![]),
        ],
        IndexRole::Ordinary,
    );
    assert!(matches!(
        run(&s, 1, EXPERIMENT_LIMITS, 3),
        Err(SolverError::ExactReplay(_))
    ));
}
fn empty_subset_and_initial_zero_skip_safely() {
    let c = context();
    let s = make(
        vec![(row(&c.one().numerator, &[0, -1]), point(2), vec![])],
        IndexRole::RequiredCut,
    );
    let c = run(&s, 2, EXPERIMENT_LIMITS, SearchOptions::default().prime).unwrap();
    assert_eq!(c.rhs.len(), 1);
    let d = c.stats.guarded_portfolio.unwrap();
    assert!(d.trials.iter().all(|t| t.completion == "empty-selector"));
    let c = context();
    let s = make(
        vec![(row(&c.one().numerator, &[0]), point(1), vec![])],
        IndexRole::Ordinary,
    );
    let c = standard(&s);
    assert!(c.rhs.is_empty());
    assert!(c.stats.guarded_portfolio.is_none());
}
fn ordinary_and_ray_paths_have_no_portfolio() {
    let c = context();
    let s = make(
        vec![(
            row(&c.one().numerator, &[0, -1]),
            IndexDomain::new([IndexBounds::new(Some(1), None).unwrap()]).unwrap(),
            vec![],
        )],
        IndexRole::Ordinary,
    );
    let solver = SectorSolver::new(&s.system, [true], Default::default()).unwrap();
    let c = solver
        .solve_case(
            CoordinateCase::new([Some(2)]).unwrap(),
            SearchOptions {
                max_depth: Some(0),
                ..Default::default()
            },
        )
        .unwrap();
    assert!(c.stats.guarded_portfolio.is_none());
    let solver = SectorSolver::new_guarded(&s.system, [true], s.roles).unwrap();
    let scope = GuardedSearchScope::new(
        &s,
        IndexDomain::new([IndexBounds::new(Some(2), None).unwrap()]).unwrap(),
    );
    let c = solver
        .solve_guarded_case(
            Case::generic(),
            SearchOptions {
                max_depth: Some(0),
                ..Default::default()
            },
            &scope,
        )
        .unwrap();
    assert!(c.stats.guarded_portfolio.is_none());
}

fn coverage_compares_original_baseline_after_conditionless_A() {
    let c = CoefficientContext::try_new(["pcut", "pn", "px"]).unwrap();
    let n = c.parameter("pn").unwrap().numerator;
    let x = c.parameter("px").unwrap().numerator;
    let one = c.one().numerator;
    let row2 = |shifts: &[[i16; 2]], coefficient: &CoefficientPolynomial| {
        shifts
            .iter()
            .map(|&a| Term {
                integral: Integral::symbolic(a).unwrap(),
                coefficient: coefficient.clone(),
            })
            .collect()
    };
    let domain = IndexDomain::new([IndexBounds::fixed(2), IndexBounds::fixed(1)]).unwrap();
    let system = Arc::new(
        GuardedSourceSystem::new(
            "original-baseline-condition-control",
            [IndexRole::RequiredCut, IndexRole::Ordinary],
            [0, 1],
            vec![
                GuardedSource::new(
                    "baseline",
                    row2(&[[0, 0], [-1, 0], [0, -1]], &one),
                    domain.clone(),
                )
                .with_nonzero_conditions(vec![x.clone()]),
                GuardedSource::new("A", row2(&[[0, 0], [0, -1]], &one), domain.clone()),
                GuardedSource::new("C", row2(&[[0, 0]], &n), domain.clone())
                    .with_nonzero_conditions(vec![x]),
            ],
        )
        .unwrap(),
    );
    let solver = SectorSolver::new_guarded(&system.system, [true, true], system.roles).unwrap();
    let scope = GuardedSearchScope::new(&system, domain.clone());
    let case: Case<2> = CoordinateCase::new([Some(2), Some(1)]).unwrap().into();
    let options = SearchOptions {
        max_depth: Some(0),
        ..Default::default()
    };
    let baseline = raw(&solver, case.clone(), options, &scope).unwrap();
    assert_eq!(baseline.rhs.len(), 2);
    let candidate = solver
        .select_guarded_portfolio(case, options, &scope, baseline, EXPERIMENT_LIMITS)
        .unwrap();
    assert!(candidate.rhs.is_empty());
    let stats = candidate.stats.guarded_portfolio.unwrap();
    assert_eq!(stats.trials[0].rhs_terms, Some(1));
    assert_eq!(stats.trials[0].nonzero_conditions, 0);
    assert_eq!(stats.trials[0].completion, "selected-strictly-shorter");
    assert_eq!(stats.trials[1].completion, "selected-strictly-shorter");
    assert_eq!(stats.selected_arm, 2);
    let rule = system
        .seal_candidate(candidate, solver.ordering().clone(), domain)
        .unwrap();
    let p = GuardedProgram::new(system.clone(), vec![rule], []).unwrap();
    let encoded = p.encode_native(Default::default()).unwrap();
    let q = GuardedProgram::decode_generated(&encoded, system, Default::default()).unwrap();
    assert_eq!(
        p.apply(&[2, 1]).unwrap().terms,
        q.apply(&[2, 1]).unwrap().terms
    );
}

fn equal_and_larger_candidates_keep_first_baseline() {
    for shifts in [&[0, -1][..], &[0, -1, -2][..]] {
        let c = context();
        let n = c.parameter("portfolio_n").unwrap().numerator;
        let system = make(
            vec![
                (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
                (row(&n, shifts), point(1), vec![]),
            ],
            IndexRole::Ordinary,
        );
        let candidate = standard(&system);
        assert_eq!(candidate.rhs.len(), 1);
        let stats = candidate.stats.guarded_portfolio.unwrap();
        assert_eq!(stats.selected_arm, 0);
        assert_eq!(stats.trials[1].completion, "not-strictly-shorter");
        assert_eq!(candidate.sources[0].basis_row, 0);
        publish(system, candidate);
    }
}
fn bounded_search_exhaustion_retains_baseline() {
    let c = context();
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let system = make(
        vec![
            (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
            (row(&n, &[2]), point(1), vec![]),
        ],
        IndexRole::Ordinary,
    );
    let candidate = standard(&system);
    assert_eq!(candidate.rhs.len(), 1);
    assert_eq!(
        candidate.stats.guarded_portfolio.unwrap().trials[1].completion,
        "search-exhausted"
    );
    publish(system, candidate);
}
fn source_zero_premise_survives_full_context_replay() {
    let c = context();
    let n = c.parameter("portfolio_n").unwrap().numerator;
    let base = make(
        vec![
            (row(&c.one().numerator, &[0, -1]), point(1), vec![]),
            (row(&n, &[0, 2]), point(1), vec![]),
        ],
        IndexRole::Ordinary,
    );
    let system = Arc::new(
        Arc::try_unwrap(base)
            .unwrap()
            .with_zero_domains(vec![point(3)])
            .unwrap(),
    );
    let candidate = standard(&system);
    assert!(candidate.rhs.is_empty());
    assert_eq!(candidate.sources[0].basis_row, 1);
    assert_eq!(system.zero_domains(), &[point(3)]);
    publish(system, candidate);
}

pub fn run_guarded_portfolio_controls() -> Vec<&'static str> {
    let controls: [(&str, fn()); 19] = [
        (
            "condition_gate_exact_and_map_safe",
            condition_gate_exact_and_map_safe,
        ),
        (
            "coverage_compares_original_baseline_after_conditionless_A",
            coverage_compares_original_baseline_after_conditionless_A,
        ),
        (
            "selectors_protect_all_roles_and_literal_indices",
            selectors_protect_all_roles_and_literal_indices,
        ),
        (
            "first_direct_shorter_zero_keeps_full_ordinals_and_roundtrips",
            first_direct_shorter_zero_keeps_full_ordinals_and_roundtrips,
        ),
        (
            "first_indirect_is_eligible_and_keeps_full_ordinals",
            first_indirect_is_eligible_and_keeps_full_ordinals,
        ),
        (
            "rational_associate_conditions_accepted_unchanged",
            rational_associate_conditions_accepted_unchanged,
        ),
        (
            "new_parameter_condition_retains_baseline",
            new_parameter_condition_retains_baseline,
        ),
        (
            "certification_failure_retains_baseline",
            certification_failure_retains_baseline,
        ),
        (
            "baseline_certification_failure_returns_unchanged_candidate",
            baseline_certification_failure_returns_unchanged_candidate,
        ),
        (
            "attempted_guard_rejections_and_accepted_budget_are_distinct",
            attempted_guard_rejections_and_accepted_budget_are_distinct,
        ),
        (
            "attempt_cap_stops_before_next_selected_instantiation",
            attempt_cap_stops_before_next_selected_instantiation,
        ),
        (
            "trace_caps_return_baseline_without_truncated_proof",
            trace_caps_return_baseline_without_truncated_proof,
        ),
        (
            "optional_power_failure_propagates",
            optional_power_failure_propagates,
        ),
        (
            "optional_exact_lift_failure_propagates",
            optional_exact_lift_failure_propagates,
        ),
        (
            "empty_subset_and_initial_zero_skip_safely",
            empty_subset_and_initial_zero_skip_safely,
        ),
        (
            "ordinary_and_ray_paths_have_no_portfolio",
            ordinary_and_ray_paths_have_no_portfolio,
        ),
        (
            "equal_and_larger_candidates_keep_first_baseline",
            equal_and_larger_candidates_keep_first_baseline,
        ),
        (
            "bounded_search_exhaustion_retains_baseline",
            bounded_search_exhaustion_retains_baseline,
        ),
        (
            "source_zero_premise_survives_full_context_replay",
            source_zero_premise_survives_full_context_replay,
        ),
    ];
    let mut passed = Vec::new();
    for (name, control) in controls {
        eprintln!("running {name}");
        control();
        eprintln!("passed {name}");
        passed.push(name);
    }
    passed
}
