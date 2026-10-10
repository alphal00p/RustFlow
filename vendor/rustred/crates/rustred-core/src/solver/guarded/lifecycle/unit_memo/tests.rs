use super::*;
use crate::algebra::CoefficientContext;
use crate::solver::guarded::{GuardedSource, IndexBounds, IndexDomain};
use crate::solver::{Integral, SearchOptions, Term};
fn discover(
    sources: Arc<GuardedSourceSystem<1>>,
    points: impl IntoIterator<Item = i64>,
    terminals: impl IntoIterator<Item = [i64; 1]>,
) -> GuardedProgram<1> {
    let mut rules = Vec::new();
    for point in points {
        let found = sources
            .solve_domains(
                vec![IndexDomain::new([IndexBounds::fixed(point)]).unwrap()],
                SearchOptions {
                    max_depth: Some(0),
                    sample_seed: 0,
                    ..Default::default()
                },
                1,
            )
            .unwrap();
        assert!(
            !found.rules.is_empty(),
            "missing direct row at {point}: {:?}",
            found.unresolved
        );
        rules.extend(found.rules);
    }
    GuardedProgram::new(sources, rules, terminals).unwrap()
}

fn diamond_program(cancel: bool) -> GuardedProgram<1> {
    let context = CoefficientContext::try_new(["guarded_diamond_n", "guarded_diamond_x"]).unwrap();
    let x = context.parameter("guarded_diamond_x").unwrap();
    let target = Integral::symbolic([0]).unwrap();
    let branches = [
        (4, vec![(-1, context.one()), (-2, context.one())]),
        (2, vec![(-1, x.clone())]),
        (3, vec![(-2, if cancel { -x.clone() } else { x.clone() })]),
    ];
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "replayed-uncovered-diamond",
            [IndexRole::Occupation],
            [0],
            branches
                .iter()
                .enumerate()
                .map(|(ordinal, (point, rhs))| {
                    let mut row = vec![Term {
                        integral: target,
                        coefficient: context.one().numerator,
                    }];
                    row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                        integral: Integral::symbolic([*shift]).unwrap(),
                        coefficient: (-coefficient.clone()).numerator,
                    }));
                    GuardedSource::new(
                        format!("diamond-{ordinal}"),
                        row,
                        IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                    )
                    .with_nonzero_conditions(vec![x.numerator.clone()])
                })
                .collect(),
        )
        .unwrap(),
    );
    discover(sources, branches.iter().map(|(point, ..)| *point), [])
}

fn guarded_pending_diamond(cancel: bool) -> GuardedProgram<1> {
    let context = CoefficientContext::try_new(["pending_n", "pending_x", "pending_y"]).unwrap();
    let x = context.parameter("pending_x").unwrap();
    let y = context.parameter("pending_y").unwrap();
    let target = Integral::symbolic([0]).unwrap();
    let branches = [
        (
            6,
            context.one(),
            vec![(-1, context.one()), (-4, context.one())],
        ),
        (5, context.one(), vec![(-1, context.one())]),
        (
            4,
            context.one(),
            vec![(
                -2,
                if cancel {
                    -context.one()
                } else {
                    context.one()
                },
            )],
        ),
        (2, y.clone(), vec![(-1, &context.one() / &y)]),
    ];
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "pending-native-order-guarded-diamond",
            [IndexRole::Occupation],
            [0],
            branches
                .iter()
                .enumerate()
                .map(|(ordinal, (point, pivot, rhs))| {
                    let mut row = vec![Term {
                        integral: target,
                        coefficient: pivot.numerator.clone(),
                    }];
                    row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                        integral: Integral::symbolic([*shift]).unwrap(),
                        coefficient: (-(coefficient * pivot)).numerator,
                    }));
                    let source = GuardedSource::new(
                        format!("pending-{ordinal}"),
                        row,
                        IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                    );
                    if *point == 5 {
                        source.with_nonzero_conditions(vec![x.numerator.clone()])
                    } else {
                        source
                    }
                })
                .collect(),
        )
        .unwrap(),
    );
    discover(sources, branches.iter().map(|(point, ..)| *point), [])
}

fn capacity_fixture() -> GuardedProgram<1> {
    let context =
        CoefficientContext::try_new(["balanced_capacity_n", "balanced_capacity_x"]).unwrap();
    let target = Integral::symbolic([0]).unwrap();
    let branches = [
        (8, vec![(-1, 1), (-2, 1)]),
        (7, vec![(-6, 1)]),
        (6, vec![(-5, -1), (-4, 1), (-3, 1)]),
    ];
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "balanced-budget-cancellation",
            [IndexRole::Occupation],
            [0],
            branches
                .iter()
                .enumerate()
                .map(|(ordinal, (point, rhs))| {
                    let mut row = vec![Term {
                        integral: target,
                        coefficient: context.one().numerator,
                    }];
                    row.extend(rhs.iter().map(|(shift, c)| Term {
                        integral: Integral::symbolic([*shift]).unwrap(),
                        coefficient: context.integer(-i64::from(*c)).numerator,
                    }));
                    GuardedSource::new(
                        format!("budget-{ordinal}"),
                        row,
                        IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                    )
                })
                .collect(),
        )
        .unwrap(),
    );
    discover(sources, branches.iter().map(|(point, _)| *point), [])
}
fn assert_equal<const N: usize>(a: &GuardedReduction<N>, b: &GuardedReduction<N>) {
    assert_eq!(a.terms, b.terms);
    assert_eq!(a.rule_applications, b.rule_applications);
    assert_eq!(a.nonzero_conditions, b.nonzero_conditions);
    assert_eq!(a.unresolved.len(), b.unresolved.len());
    for (x, y) in a.unresolved.iter().zip(&b.unresolved) {
        assert_eq!(x.integral, y.integral);
        assert_eq!(x.reason, y.reason);
        assert_eq!(x.coefficient, y.coefficient);
    }
}
fn enabled<const N: usize>(p: GuardedProgram<N>) -> GuardedProgram<N> {
    p.with_unit_reduction_memo(Some(Default::default()))
        .unwrap()
}
#[test]
fn native_memo_terminal_promotion_preserves_ordered_guards_and_decode() {
    for cancel in [true, false] {
        let p = enabled(guarded_pending_diamond(cancel));
        for t in [[5], [2], [6]] {
            let got = p.reduce_memoized(t, Default::default()).unwrap();
            assert!(!got.cache_hit);
            assert_equal(&got.reduction, &p.reduce(t, Default::default()).unwrap());
        }
        let p = p
            .with_promoted_terminals_replayed([[1]], usize::MAX)
            .unwrap();
        let got = p.reduce_memoized([6], Default::default()).unwrap();
        assert!(got.cache_hit);
        assert_equal(&got.reduction, &p.reduce([6], Default::default()).unwrap());
        let c = CoefficientContext::try_new(["pending_n", "pending_x", "pending_y"]).unwrap();
        assert!(
            got.reduction
                .nonzero_conditions
                .contains(&c.parameter("pending_x").unwrap().numerator)
        );
        assert_eq!(
            got.reduction
                .nonzero_conditions
                .contains(&c.parameter("pending_y").unwrap().numerator),
            !cancel
        );
        let decoded = GuardedProgram::decode_generated(
            &p.encode_native(Default::default()).unwrap(),
            p.sources().clone(),
            Default::default(),
        )
        .unwrap();
        assert!(decoded.unit_reduction_memo_usage().unwrap().is_none());
        assert_equal(
            &got.reduction,
            &decoded.reduce([6], Default::default()).unwrap(),
        );
    }
}
#[test]
fn native_memo_budget_grid_matches_cold_reduction() {
    for applications in 0..=5 {
        for pending in 0..=4 {
            for entries in [0, 1, 4] {
                let p = capacity_fixture()
                    .with_unit_reduction_memo(Some(GuardedUnitMemoLimits {
                        max_entries: entries,
                        ..Default::default()
                    }))
                    .unwrap();
                let limits = GuardedReductionLimits {
                    max_rule_applications: applications,
                    max_pending_integrals: pending,
                };
                assert_equal(
                    &p.reduce_memoized([8], limits).unwrap().reduction,
                    &p.reduce([8], limits).unwrap(),
                );
                let p = p
                    .with_promoted_terminals_replayed([[1], [2], [3]], usize::MAX)
                    .unwrap();
                for _ in 0..2 {
                    assert_equal(
                        &p.reduce_memoized([8], limits).unwrap().reduction,
                        &p.reduce([8], limits).unwrap(),
                    );
                }
            }
        }
    }
}
#[test]
fn native_memo_work_limit_cannot_be_reclassified_as_complete() {
    let p = enabled(diamond_program(false));
    let limits = GuardedReductionLimits {
        max_rule_applications: 3,
        max_pending_integrals: 100,
    };
    let old = p.reduce_memoized([4], limits).unwrap();
    assert!(
        old.reduction
            .unresolved
            .iter()
            .any(|r| r.reason == GuardedApplicationFailure::WorkLimit)
    );
    assert_eq!(p.unit_reduction_memo_usage().unwrap().unwrap().entries, 0);
    let p = p
        .with_promoted_terminals_replayed([[1]], usize::MAX)
        .unwrap();
    let new = p.reduce_memoized([4], limits).unwrap();
    assert!(!new.cache_hit);
    assert!(new.reduction.unresolved.is_empty());
    assert_equal(&new.reduction, &p.reduce([4], limits).unwrap());
}
#[test]
fn native_memo_rebind_union_reset_and_removal_invalidate() {
    let p = enabled(diamond_program(false));
    p.reduce_memoized([4], Default::default()).unwrap();
    let fallback = GuardedProgram::new(p.sources().clone(), vec![], []).unwrap();
    let p = p.union_verified(fallback, [], usize::MAX).unwrap();
    assert!(p.unit_reduction_memo_usage().unwrap().is_none());
    assert!(
        !p.reduce_memoized([4], Default::default())
            .unwrap()
            .cache_hit
    );
    let p = enabled(p);
    p.reduce_memoized([4], Default::default()).unwrap();
    let p = p.with_terminals_verified([], usize::MAX).unwrap();
    assert!(p.unit_reduction_memo_usage().unwrap().is_none());
    assert!(
        enabled(diamond_program(false))
            .with_promoted_terminals_replayed([[3]], usize::MAX)
            .is_err()
    );
    assert!(
        enabled(diamond_program(false))
            .with_promoted_terminals_replayed([[-1]], usize::MAX)
            .is_err()
    );
    assert!(
        enabled(diamond_program(false))
            .with_promoted_terminals_replayed([[1]], 0)
            .is_err()
    );
    let p = diamond_program(false)
        .with_terminals_replayed([[1]], usize::MAX)
        .unwrap();
    assert!(
        enabled(p)
            .with_promoted_terminals_replayed([], usize::MAX)
            .is_err()
    );
}
#[test]
fn native_memo_quotas_limits_eviction_and_owned_bytes_are_explicit() {
    let caps = GuardedUnitMemoLimits {
        max_entries: 1,
        ..Default::default()
    };
    let p = capacity_fixture()
        .with_unit_reduction_memo(Some(caps))
        .unwrap();
    p.reduce_memoized([8], Default::default()).unwrap();
    p.reduce_memoized([7], Default::default()).unwrap();
    assert!(
        !p.reduce_memoized([8], Default::default())
            .unwrap()
            .cache_hit
    );
    let u = p.unit_reduction_memo_usage().unwrap().unwrap();
    assert!(
        u.entries <= 1
            && u.output_terms <= caps.max_output_terms
            && u.conditions <= caps.max_conditions
            && u.polynomial_terms <= caps.max_polynomial_terms
            && u.owned_bytes <= caps.max_owned_bytes
    );
    let limits = GuardedReductionLimits {
        max_pending_integrals: 99999,
        ..Default::default()
    };
    assert!(!p.reduce_memoized([8], limits).unwrap().cache_hit);
    for c in [
        GuardedUnitMemoLimits {
            max_output_terms: 0,
            ..Default::default()
        },
        GuardedUnitMemoLimits {
            max_conditions: 0,
            ..Default::default()
        },
        GuardedUnitMemoLimits {
            max_polynomial_terms: 0,
            ..Default::default()
        },
    ] {
        let p = diamond_program(false)
            .with_unit_reduction_memo(Some(c))
            .unwrap();
        for _ in 0..2 {
            assert!(
                !p.reduce_memoized([4], Default::default())
                    .unwrap()
                    .cache_hit
            );
        }
    }
    assert!(
        capacity_fixture()
            .with_unit_reduction_memo(Some(GuardedUnitMemoLimits {
                max_owned_bytes: 0,
                ..Default::default()
            }))
            .is_err()
    );
}
#[test]
fn native_memo_encoded_payload_cap_declines_without_changing_the_result() {
    let caps = GuardedUnitMemoLimits {
        max_entries: 1,
        ..Default::default()
    };
    let baseline = diamond_program(false);
    let expected = baseline.reduce([4], Default::default()).unwrap();
    let empty = baseline.with_unit_reduction_memo(Some(caps)).unwrap();
    let metadata_bytes = empty
        .unit_reduction_memo_usage()
        .unwrap()
        .unwrap()
        .owned_bytes;
    // The metadata fits but even the smallest encoded nonzero entry does not.
    let p = empty
        .with_unit_reduction_memo(Some(GuardedUnitMemoLimits {
            max_owned_bytes: metadata_bytes + 1,
            ..caps
        }))
        .unwrap();
    for _ in 0..2 {
        let got = p.reduce_memoized([4], Default::default()).unwrap();
        assert!(!got.cache_hit);
        assert_equal(&got.reduction, &expected);
        let usage = p.unit_reduction_memo_usage().unwrap().unwrap();
        assert_eq!(usage.entries, 0);
        assert_eq!(usage.owned_bytes, metadata_bytes);
    }
}
#[test]
fn native_memo_cannot_cross_source_context_or_survive_reset() {
    let p = enabled(diamond_program(false));
    p.reduce_memoized([4], Default::default()).unwrap();
    // Same names/roles/measure, different exact source row: union must reject.
    assert!(
        p.union_verified(diamond_program(true), [], usize::MAX)
            .is_err()
    );
    let p = enabled(diamond_program(false));
    p.reduce_memoized([4], Default::default()).unwrap();
    let p = p
        .with_unit_reduction_memo(Some(Default::default()))
        .unwrap();
    assert_eq!(p.unit_reduction_memo_usage().unwrap().unwrap().entries, 0);
    assert!(
        !p.reduce_memoized([4], Default::default())
            .unwrap()
            .cache_hit
    );
    let p = p.with_terminals_replayed([], usize::MAX).unwrap();
    assert!(p.unit_reduction_memo_usage().unwrap().is_none());
    let p = enabled(p);
    p.reduce_memoized([4], Default::default()).unwrap();
    let fallback = GuardedProgram::new(p.sources().clone(), vec![], []).unwrap();
    let p = p.union_replayed(fallback, [], usize::MAX).unwrap();
    assert!(p.unit_reduction_memo_usage().unwrap().is_none());
}
#[test]
fn native_memo_owned_decode_error_propagates() {
    let p = enabled(diamond_program(false));
    p.reduce_memoized([4], Default::default()).unwrap();
    {
        let mut memo = p.unit_memo.as_ref().unwrap().lock().unwrap();
        let entry = memo.slots.iter_mut().flatten().next().unwrap();
        let mut bytes = entry.atoms.to_vec();
        bytes.push(0);
        entry.atoms = bytes.into_boxed_slice();
    }
    // Only outer framing is corrupted. A hit must fail, not silently recompute.
    assert!(p.reduce_memoized([4], Default::default()).is_err());
    assert!(p.reduce([4], Default::default()).is_ok());
}
#[test]
fn native_memo_mixed_promoted_and_uncovered_outputs_match() {
    let p = enabled(capacity_fixture());
    for t in [[8], [7], [6]] {
        p.reduce_memoized(t, Default::default()).unwrap();
    }
    let p = p
        .with_promoted_terminals_replayed([[1], [2]], usize::MAX)
        .unwrap();
    for t in [[8], [7], [6]] {
        let r = p.reduce_memoized(t, Default::default()).unwrap();
        assert!(r.cache_hit);
        assert_equal(&r.reduction, &p.reduce(t, Default::default()).unwrap());
    }
    let p = enabled(diamond_program(true));
    let before = p.reduce_memoized([4], Default::default()).unwrap();
    assert!(before.reduction.terms.is_empty() && before.reduction.unresolved.is_empty());
    let p = p
        .with_promoted_terminals_replayed([[1]], usize::MAX)
        .unwrap();
    assert_equal(
        &before.reduction,
        &p.reduce_memoized([4], Default::default())
            .unwrap()
            .reduction,
    );
}
#[test]
fn native_memo_coupled_failure_and_corrupt_proof_cannot_be_promoted() {
    let c = CoefficientContext::try_new(["memo_guard_n", "memo_guard_m"]).unwrap();
    let n = c.parameter("memo_guard_n").unwrap();
    let m = c.parameter("memo_guard_m").unwrap();
    let domain = IndexDomain::new([
        IndexBounds::new(Some(1), None).unwrap(),
        IndexBounds::new(Some(1), None).unwrap(),
    ])
    .unwrap();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "memo-coupled-guard",
            [IndexRole::Occupation; 2],
            [0, 1],
            vec![GuardedSource::new(
                "coupled",
                vec![
                    Term {
                        integral: Integral::symbolic([0, 0]).unwrap(),
                        coefficient: (&n - &m).numerator,
                    },
                    Term {
                        integral: Integral::symbolic([-1, 0]).unwrap(),
                        coefficient: (-c.one()).numerator,
                    },
                ],
                domain.clone(),
            )],
        )
        .unwrap(),
    );
    let found = sources
        .solve_domains(
            vec![domain],
            SearchOptions {
                max_depth: Some(1),
                sample_seed: 0,
                ..Default::default()
            },
            1,
        )
        .unwrap();
    let p = enabled(GuardedProgram::new(sources, found.rules, []).unwrap());
    for _ in 0..2 {
        let r = p.reduce_memoized([2, 2], Default::default()).unwrap();
        assert!(!r.cache_hit);
        assert!(matches!(
            r.reduction.unresolved[0].reason,
            GuardedApplicationFailure::ConditionVanished { .. }
        ));
    }
    assert!(
        p.with_promoted_terminals_replayed([[2, 2]], usize::MAX)
            .is_err()
    );
    let mut p = enabled(capacity_fixture());
    p.reduce_memoized([8], Default::default()).unwrap();
    p.rules[0].candidate.rhs[0].coefficient = -p.rules[0].candidate.rhs[0].coefficient.clone();
    assert!(matches!(
        p.with_promoted_terminals_replayed([[1], [2], [3]], usize::MAX),
        Err(SolverError::Certification(_))
    ));
}
