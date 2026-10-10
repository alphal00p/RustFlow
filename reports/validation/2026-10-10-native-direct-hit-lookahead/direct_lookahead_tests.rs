//! Isolated experiment controls. No physical periods or analytic zeros supplied.
use super::*;
use crate::algebra::CoefficientContext;
use crate::solver::{CoordinateCase, Power};
use crate::solver::guarded::{
    GuardedApplicationStatus, GuardedProgram, GuardedSearchScope, GuardedSource,
    GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole,
};

fn point(n: i64) -> IndexDomain<1> {
    IndexDomain::new([IndexBounds::fixed(n)]).unwrap()
}
fn fixture(rows: &[(&[i16], i64)]) -> Arc<GuardedSourceSystem<1>> {
    let c = CoefficientContext::new(["lookahead_n", "lookahead_x"]);
    let sources = rows
        .iter()
        .enumerate()
        .map(|(i, (shifts, at))| {
            GuardedSource::new(
                format!("source-{i}"),
                shifts
                    .iter()
                    .map(|&s| Term {
                        integral: Integral::symbolic([s]).unwrap(),
                        coefficient: c.one().numerator,
                    })
                    .collect(),
                point(*at),
            )
            .with_nonzero_conditions(vec![c.parameter("lookahead_x").unwrap().numerator])
        })
        .collect();
    Arc::new(
        GuardedSourceSystem::new(
            "isolated-direct-lookahead",
            [IndexRole::Ordinary],
            [0],
            sources,
        )
        .unwrap(),
    )
}
fn run(
    system: &GuardedSourceSystem<1>,
    cap: Option<DirectLookaheadLimits>,
    depth: u32,
) -> Result<RuleCandidate<1>, TrialSearchError> {
    let solver = SectorSolver::new_guarded(&system.system, [true], system.roles).unwrap();
    let scope = GuardedSearchScope::new(system, point(1));
    let mut work = RuleTrialStats::default();
    solver.search_attempt_inner(
        CoordinateCase::new([Some(1)]).unwrap().into(),
        SearchOptions {
            max_depth: Some(depth),
            ..Default::default()
        },
        None,
        None,
        false,
        &mut work,
        |_| {},
        Some(&scope),
        cap,
    )
}
fn success(result: Result<RuleCandidate<1>, TrialSearchError>) -> RuleCandidate<1> {
    match result {
        Ok(c) => c,
        Err(TrialSearchError::Solver(e)) => panic!("{e}"),
        Err(TrialSearchError::Limit(e)) => panic!("{e:?}"),
    }
}
fn caps(rows: usize, trace: usize, terms: usize) -> Option<DirectLookaheadLimits> {
    Some(DirectLookaheadLimits {
        extra_rows: rows,
        trace_rows: trace,
        trace_terms: terms,
    })
}

#[test]
fn later_native_indirect_zero_beats_direct_and_roundtrips_with_conditions() {
    // Baseline I(1)=-I(0); later I(2)+I(1)=0 and I(2)=0 imply I(1)=0.
    // Neither later row is a target-leading source at the initial point.
    let system = fixture(&[(&[0, -1], 1), (&[1, 0], 1), (&[1], 1)]);
    let baseline = success(run(&system, None, 0));
    assert_eq!(baseline.rhs.len(), 1);
    let candidate = success(run(&system, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert!(candidate.rhs.is_empty());
    assert!(!candidate.stats.direct_hit);
    assert_eq!(candidate.sources.len(), 2);
    let stats = candidate.stats.direct_lookahead.unwrap();
    assert_eq!(stats.extra_attempted_rows, 2);
    assert_eq!(stats.indirect_trace_rows, 2);
    assert_eq!(stats.completion, "strictly-smaller-indirect");
    let found = system
        .solve_domain(
            point(1),
            SearchOptions {
                max_depth: Some(0),
                ..Default::default()
            },
        )
        .unwrap();
    let program = GuardedProgram::new(system.clone(), found.rules, []).unwrap();
    let bytes = program.encode_native(Default::default()).unwrap();
    let decoded = GuardedProgram::decode_generated(&bytes, system, Default::default()).unwrap();
    let applied = decoded.apply(&[1]).unwrap();
    assert!(matches!(
        applied.status,
        GuardedApplicationStatus::Applied { .. }
    ));
    assert!(applied.terms.is_empty());
    assert_eq!(applied.nonzero_conditions.len(), 1);
}

#[test]
fn extra_row_limit_preserves_identical_baseline_proof() {
    let system = fixture(&[(&[0, -1], 1), (&[1, 0], 1), (&[1], 1)]);
    let baseline = success(run(&system, None, 0));
    for n in [0, 1] {
        let got = success(run(&system, caps(n, 64, 16384), 0));
        assert_eq!(got.rhs, baseline.rhs);
        assert_eq!(got.sources, baseline.sources);
        assert!(got.stats.direct_hit);
        assert_eq!(got.stats.exact_trace_rows, 1);
        assert_eq!(got.stats.direct_lookahead.unwrap().extra_attempted_rows, n);
        assert_eq!(
            got.stats.direct_lookahead.unwrap().completion,
            "source-row-budget-fallback"
        );
    }
}

#[test]
fn guard_rejected_attempts_consume_extra_budget() {
    let system = fixture(&[(&[0, -1], 1), (&[1, 0], 2), (&[1], 1)]);
    let got = success(run(&system, caps(1, 64, 16384), 0));
    let extra = got.stats.direct_lookahead.unwrap();
    assert_eq!(extra.extra_attempted_rows, 1);
    assert_eq!(extra.guard_rejected_rows, 1);
    assert_eq!(extra.completion, "source-row-budget-fallback");
    assert_eq!(got.rhs.len(), 1);
}

#[test]
fn trace_caps_return_complete_direct_proof_never_truncated_lift() {
    let system = fixture(&[(&[0, -1], 1), (&[1, 0], 1), (&[1], 1)]);
    for cap in [caps(512, 1, 16384), caps(512, 64, 1)] {
        let got = success(run(&system, cap, 0));
        assert_eq!(got.rhs.len(), 1);
        assert_eq!(got.sources.len(), 1);
        assert!(got.stats.direct_hit);
        assert_eq!(
            got.stats.direct_lookahead.unwrap().completion,
            "indirect-trace-budget-fallback"
        );
    }
}

#[test]
fn first_indirect_and_initial_zero_do_not_launch_lookahead() {
    let indirect = fixture(&[(&[1, 0], 1), (&[1], 1)]);
    let got = success(run(&indirect, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert!(got.rhs.is_empty());
    assert!(!got.stats.direct_hit);
    assert!(got.stats.direct_lookahead.is_none());
    let zero = fixture(&[(&[0], 1), (&[1, 0], 1), (&[1], 1)]);
    let got = success(run(&zero, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert!(got.rhs.is_empty());
    assert!(got.stats.direct_hit);
    assert!(got.stats.direct_lookahead.is_none());
}

#[test]
fn larger_indirect_rhs_and_depth_exhaustion_keep_fallback() {
    let larger = fixture(&[(&[0, -1], 1), (&[1, 0, -1, -2], 1), (&[1], 1)]);
    let got = success(run(&larger, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert_eq!(got.rhs.len(), 1);
    assert!(got.stats.direct_hit);
    assert_eq!(
        got.stats.direct_lookahead.unwrap().indirect_rhs_terms,
        Some(2)
    );
    assert_eq!(
        got.stats.direct_lookahead.unwrap().completion,
        "indirect-not-smaller-fallback"
    );
    let exhausted = fixture(&[(&[0, -1], 1)]);
    let got = success(run(&exhausted, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert_eq!(got.rhs.len(), 1);
    assert_eq!(
        got.stats.direct_lookahead.unwrap().completion,
        "depth-exhaustion-fallback"
    );
}

#[test]
fn later_direct_zero_is_immediately_optimal() {
    let system = fixture(&[(&[0, -1], 1), (&[0], 1), (&[1, 0], 1), (&[1], 1)]);
    let got = success(run(&system, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert!(got.rhs.is_empty());
    assert!(got.stats.direct_hit);
    assert_eq!(got.sources.len(), 1);
    let extra = got.stats.direct_lookahead.unwrap();
    assert_eq!(extra.extra_attempted_rows, 1);
    assert_eq!(extra.completion, "exact-direct-zero");
}

#[test]
fn ordinary_first_direct_and_no_direct_search_do_not_use_experiment() {
    for (rows, want_direct, want_terms) in [
        (
            vec![(&[0, -1][..], 1), (&[1, 0][..], 1), (&[1][..], 1)],
            true,
            1,
        ),
        (vec![(&[1, 0][..], 1), (&[1][..], 1)], false, 0),
    ] {
        let system = fixture(&rows);
        let solver = SectorSolver {
            system: &system.system,
            basis: system.system.rows().to_vec(),
            order: IntegralOrder::new([true], [false]),
            config: SectorConfig::default(),
        };
        let got = solver
            .solve_case(
                CoordinateCase::new([Some(1)]).unwrap(),
                SearchOptions {
                    max_depth: Some(0),
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(got.stats.direct_hit, want_direct);
        assert_eq!(got.rhs.len(), want_terms);
        assert!(got.stats.direct_lookahead.is_none());
    }
}

#[test]
fn optional_fatal_power_error_is_not_hidden_by_direct_fallback() {
    let system = fixture(&[(&[0, -1], 1), (&[Power::MAX], 1)]);
    let baseline = success(run(&system, None, 0));
    assert_eq!(baseline.rhs.len(), 1);
    assert!(matches!(
        run(&system, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0),
        Err(TrialSearchError::Solver(SolverError::Power(
            crate::solver::PowerError::OutOfRange { value: 64 }
        )))
    ));
}

#[test]
fn equal_size_indirect_rhs_keeps_original_direct_proof() {
    let system = fixture(&[(&[0, -1], 1), (&[1, 0], 1), (&[1, -2], 1)]);
    let baseline = success(run(&system, None, 0));
    let got = success(run(&system, Some(EXPERIMENT_DIRECT_LOOKAHEAD), 0));
    assert_eq!(got.rhs, baseline.rhs);
    assert_eq!(got.sources, baseline.sources);
    assert_eq!(
        got.stats.direct_lookahead.unwrap().indirect_rhs_terms,
        Some(1)
    );
    assert_eq!(
        got.stats.direct_lookahead.unwrap().completion,
        "indirect-not-smaller-fallback"
    );
}

#[test]
fn smaller_indirect_retains_its_own_nonzero_conditions_after_sealing() {
    let c = CoefficientContext::new([
        "lookahead_condition_n",
        "lookahead_condition_x",
        "lookahead_condition_y",
    ]);
    let specs: Vec<(&[i16], &str)> = vec![
        (&[0, -1], "lookahead_condition_x"),
        (&[1, 0], "lookahead_condition_y"),
        (&[1], "lookahead_condition_y"),
    ];
    let rows = specs
        .iter()
        .enumerate()
        .map(|(i, (shifts, condition))| {
            GuardedSource::new(
                format!("conditional-{i}"),
                shifts
                    .iter()
                    .map(|&s| Term {
                        integral: Integral::symbolic([s]).unwrap(),
                        coefficient: c.one().numerator,
                    })
                    .collect(),
                point(1),
            )
            .with_nonzero_conditions(vec![c.parameter(condition).unwrap().numerator])
        })
        .collect();
    let system = Arc::new(
        GuardedSourceSystem::new(
            "isolated-lookahead-conditional-choice",
            [IndexRole::Ordinary],
            [0],
            rows,
        )
        .unwrap(),
    );
    let found = system
        .solve_domain(
            point(1),
            SearchOptions {
                max_depth: Some(0),
                ..Default::default()
            },
        )
        .unwrap();
    let program = GuardedProgram::new(system.clone(), found.rules, []).unwrap();
    let bytes = program.encode_native(Default::default()).unwrap();
    let decoded = GuardedProgram::decode_generated(&bytes, system, Default::default()).unwrap();
    let applied = decoded.apply(&[1]).unwrap();
    assert!(applied.terms.is_empty());
    assert!(
        applied
            .nonzero_conditions
            .contains(&c.parameter("lookahead_condition_y").unwrap().numerator)
    );
}
