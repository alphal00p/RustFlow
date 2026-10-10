use super::*;
fn batch_limits(n: usize) -> GuardedUnitBatchLimits {
    GuardedUnitBatchLimits {
        max_workers: n.min(3),
        max_in_flight: n,
    }
}
fn compare_items<const N: usize>(
    batch: &GuardedUnitBatch<N>,
    serial: &[Result<GuardedMemoizedReduction<N>, SolverError>],
) {
    assert_eq!(batch.items.len(), serial.len());
    for (got, expected) in batch.items.iter().zip(serial) {
        match (got, expected) {
            (Ok(g), Ok(e)) => {
                assert_eq!(g.cache_hit, e.cache_hit);
                assert_equal(&g.reduction, &e.reduction);
            }
            (Err(g), Err(e)) => assert_eq!(g.to_string(), e.to_string()),
            _ => panic!("ordered native result differs: {got:?} vs {expected:?}"),
        }
    }
    assert_eq!(
        batch.stats.returned_logical_rule_applications,
        batch.stats.returned_uncached_rule_applications
            + batch.stats.returned_memoized_rule_applications
    );
}
fn entry_order<const N: usize>(p: &GuardedProgram<N>) -> Vec<([i64; N], usize, usize)> {
    p.unit_memo
        .as_ref()
        .unwrap()
        .lock()
        .unwrap()
        .slots
        .iter()
        .flatten()
        .map(|e| {
            (
                e.target,
                e.limits.max_rule_applications,
                e.limits.max_pending_integrals,
            )
        })
        .collect()
}
#[test]
fn bounded_batch_preserves_ordered_conditions_and_canceled_unused_guards() {
    for cancel in [true, false] {
        let serial = enabled(guarded_pending_diamond(cancel));
        let parallel = enabled(guarded_pending_diamond(cancel));
        let labels = [[6], [5], [2], [6]];
        let limits = Default::default();
        let expected = labels
            .iter()
            .map(|x| serial.reduce_memoized(*x, limits))
            .collect::<Vec<_>>();
        let got = parallel
            .reduce_batch_memoized(&labels, limits, batch_limits(4))
            .unwrap();
        compare_items(&got, &expected);
        assert_eq!(entry_order(&serial), entry_order(&parallel));
        assert_same_memo(&serial, &parallel);
        assert_eq!(got.stats.unused_prefetch_known_rule_applications, 0);
        assert_eq!(
            got.stats.actual_known_rule_applications,
            got.stats.returned_uncached_rule_applications
        );
        assert!(got.stats.actual_application_count_complete);
    }
}
#[test]
fn bounded_batch_snapshot_hit_eviction_and_duplicate_miss_match_serial_lru() {
    let make = || {
        super::super::super::tests::sample("batch-lru")
            .with_unit_reduction_memo(Some(GuardedUnitMemoLimits {
                max_entries: 1,
                ..Default::default()
            }))
            .unwrap()
    };
    let serial = make();
    let parallel = make();
    let limits = Default::default();
    serial.reduce_memoized([3], limits).unwrap();
    parallel.reduce_memoized([3], limits).unwrap();
    let before = parallel.encode_native(Default::default()).unwrap();
    let labels = [[4], [3], [4], [4]];
    let expected = labels
        .iter()
        .map(|x| serial.reduce_memoized(*x, limits))
        .collect::<Vec<_>>();
    let got = parallel
        .reduce_batch_memoized(&labels, limits, batch_limits(4))
        .unwrap();
    compare_items(&got, &expected);
    assert_eq!(entry_order(&serial), entry_order(&parallel));
    assert_same_memo(&serial, &parallel);
    assert_eq!(got.stats.prefetched_native_calls, 1);
    assert_eq!(got.stats.fallback_native_calls, 2);
    assert_eq!(got.stats.cache_hits, 1);
    assert_eq!(before, parallel.encode_native(Default::default()).unwrap());
}
#[test]
fn bounded_batch_unretained_duplicates_execute_again_and_plain_calls_do_not_deduplicate() {
    for memo in [
        None,
        Some(GuardedUnitMemoLimits {
            max_entries: 0,
            ..Default::default()
        }),
    ] {
        let serial = diamond_program(false)
            .with_unit_reduction_memo(memo)
            .unwrap();
        let parallel = diamond_program(false)
            .with_unit_reduction_memo(memo)
            .unwrap();
        let labels = [[4]; 4];
        let limits = Default::default();
        let expected = labels
            .iter()
            .map(|x| serial.reduce_memoized(*x, limits))
            .collect::<Vec<_>>();
        let got = parallel
            .reduce_batch_memoized(&labels, limits, batch_limits(4))
            .unwrap();
        compare_items(&got, &expected);
        assert_eq!(got.stats.cache_hits, 0);
        assert_eq!(
            got.stats.prefetched_native_calls + got.stats.fallback_native_calls,
            4
        );
        assert_eq!(
            got.stats.actual_known_rule_applications,
            got.stats.returned_logical_rule_applications
        );
    }
}
#[test]
fn bounded_batch_keeps_work_limits_invalid_occupations_and_norule_distinct() {
    for applications in 0..=4 {
        for pending in 0..=3 {
            let serial = enabled(capacity_fixture());
            let parallel = enabled(capacity_fixture());
            let labels = [[8], [-1], [1], [8]];
            let limits = GuardedReductionLimits {
                max_rule_applications: applications,
                max_pending_integrals: pending,
            };
            let expected = labels
                .iter()
                .map(|x| serial.reduce_memoized(*x, limits))
                .collect::<Vec<_>>();
            let got = parallel
                .reduce_batch_memoized(&labels, limits, batch_limits(4))
                .unwrap();
            compare_items(&got, &expected);
            assert_eq!(entry_order(&serial), entry_order(&parallel));
            assert_same_memo(&serial, &parallel);
        }
    }
}
#[test]
fn bounded_batch_exact_limit_keys_and_memo_caps_survive_ordered_commit() {
    let make = || enabled(super::super::super::tests::sample("batch-limit-keys"));
    let serial = make();
    let parallel = make();
    let low = GuardedReductionLimits {
        max_rule_applications: 2,
        max_pending_integrals: 1,
    };
    for p in [&serial, &parallel] {
        let first = p.reduce_memoized([2], low).unwrap();
        assert!(first.reduction.unresolved.is_empty());
        assert!(p.reduce_memoized([2], low).unwrap().cache_hit);
    }
    let labels = [[2], [1], [2]];
    let limits = Default::default();
    let expected = labels
        .iter()
        .map(|x| serial.reduce_memoized(*x, limits))
        .collect::<Vec<_>>();
    let got = parallel
        .reduce_batch_memoized(&labels, limits, batch_limits(3))
        .unwrap();
    compare_items(&got, &expected);
    assert_eq!(entry_order(&serial), entry_order(&parallel));
    assert_same_memo(&serial, &parallel);
    let keys = entry_order(&parallel);
    assert!(keys.contains(&([2], low.max_rule_applications, low.max_pending_integrals)));
    assert!(keys.contains(&(
        [2],
        limits.max_rule_applications,
        limits.max_pending_integrals
    )));
    assert!(!same_limits(low, limits));
    let caps = GuardedUnitMemoLimits::default();
    let usage = parallel.unit_reduction_memo_usage().unwrap().unwrap();
    assert!(
        usage.entries <= caps.max_entries
            && usage.owned_bytes <= caps.max_owned_bytes
            && usage.output_terms <= caps.max_output_terms
            && usage.conditions <= caps.max_conditions
            && usage.polynomial_terms <= caps.max_polynomial_terms
    );
}
#[test]
fn bounded_batch_reports_owned_codec_error_per_item_and_later_executed_work() {
    let make = || {
        let p = enabled(super::super::super::tests::sample("batch-codec"));
        p.reduce_memoized([2], Default::default()).unwrap();
        p.unit_memo
            .as_ref()
            .unwrap()
            .lock()
            .unwrap()
            .slots
            .iter_mut()
            .flatten()
            .find(|e| e.target == [2])
            .unwrap()
            .state = vec![255].into_boxed_slice();
        p
    };
    let serial = make();
    let parallel = make();
    let limits = Default::default();
    let labels = [[2], [3]];
    let expected = labels
        .iter()
        .map(|x| serial.reduce_memoized(*x, limits))
        .collect::<Vec<_>>();
    let got = parallel
        .reduce_batch_memoized(&labels, limits, batch_limits(2))
        .unwrap();
    compare_items(&got, &expected);
    assert!(got.items[0].is_err());
    assert!(got.items[1].is_ok());
    assert_eq!(got.stats.prefetched_native_calls, 1);
    assert_eq!(got.stats.actual_known_rule_applications, 3);
    // This is an explicit non-fail-fast batch. A caller returning the first
    // logical error must still report the later three executed applications.
}
#[test]
fn bounded_batch_rejects_bounds_before_touching_cache_or_work() {
    let p = enabled(diamond_program(false));
    let before = p.unit_reduction_memo_usage().unwrap().unwrap().entries;
    for b in [
        GuardedUnitBatchLimits {
            max_workers: 0,
            max_in_flight: 2,
        },
        GuardedUnitBatchLimits {
            max_workers: 3,
            max_in_flight: 2,
        },
        GuardedUnitBatchLimits {
            max_workers: 1,
            max_in_flight: 1,
        },
    ] {
        assert!(
            p.reduce_batch_memoized(&[[4], [4]], Default::default(), b)
                .is_err()
        );
    }
    assert_eq!(
        p.unit_reduction_memo_usage().unwrap().unwrap().entries,
        before
    );
    let empty = p
        .reduce_batch_memoized(&[], Default::default(), Default::default())
        .unwrap();
    assert!(empty.items.is_empty());
    assert_eq!(empty.stats.actual_known_rule_applications, 0);
}

fn assert_same_memo<const N: usize>(a: &GuardedProgram<N>, b: &GuardedProgram<N>) {
    let a = a.unit_memo.as_ref().unwrap().lock().unwrap();
    let b = b.unit_memo.as_ref().unwrap().lock().unwrap();
    assert_eq!(a.slots.len(), b.slots.len());
    for (a, b) in a.slots.iter().zip(b.slots.iter()) {
        match (a, b) {
            (None, None) => {}
            (Some(a), Some(b)) => {
                assert_eq!(a.target, b.target);
                assert!(same_limits(a.limits, b.limits));
                assert_eq!(a.state, b.state);
                assert_eq!(a.atoms, b.atoms);
                assert_eq!(format!("{:?}", a.records), format!("{:?}", b.records));
                assert_eq!(format!("{:?}", a.conditions), format!("{:?}", b.conditions));
                assert_eq!(a.applications, b.applications);
                assert_eq!(a.polynomial_terms, b.polynomial_terms);
            }
            _ => panic!("ordered empty/occupied memo slot differs"),
        }
    }
}
#[test]
fn bounded_batch_duplicate_coupled_guard_failures_are_not_cached_or_canceled() {
    let make = || {
        let c = CoefficientContext::try_new(["memo_guard_n", "memo_guard_m"]).unwrap();
        let n = c.parameter("memo_guard_n").unwrap();
        let m = c.parameter("memo_guard_m").unwrap();
        let domain = IndexDomain::new([IndexBounds::new(Some(1), None).unwrap(); 2]).unwrap();
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
        enabled(GuardedProgram::new(sources, found.rules, []).unwrap())
    };
    let serial = make();
    let parallel = make();
    let labels = [[2, 2], [2, 2], [3, 2]];
    let limits = Default::default();
    let expected = labels
        .iter()
        .map(|x| serial.reduce_memoized(*x, limits))
        .collect::<Vec<_>>();
    let got = parallel
        .reduce_batch_memoized(&labels, limits, batch_limits(3))
        .unwrap();
    compare_items(&got, &expected);
    assert_same_memo(&serial, &parallel);
    for result in &got.items[..2] {
        let r = result.as_ref().unwrap();
        assert!(!r.cache_hit);
        assert!(matches!(
            r.reduction.unresolved[0].reason,
            GuardedApplicationFailure::ConditionVanished { .. }
        ));
    }
    assert!(got.stats.fallback_native_calls >= 1);
}

#[test]
fn bounded_batch_native_error_marks_actual_count_incomplete_and_continues() {
    // Deliberately corrupt a private verified field to exercise the error path.
    // No public constructor admits this invalid polynomial/index map.
    let make = || {
        let mut p = enabled(super::super::super::tests::sample("batch-native-error"));
        let wrong = CoefficientContext::try_new(Vec::<String>::new()).unwrap();
        p.rules[0].nonzero_conditions = vec![wrong.one().numerator];
        p
    };
    let serial = make();
    let parallel = make();
    let labels = [[2], [2], [0]];
    let limits = Default::default();
    let expected = labels
        .iter()
        .map(|x| serial.reduce_memoized(*x, limits))
        .collect::<Vec<_>>();
    let got = parallel
        .reduce_batch_memoized(&labels, limits, batch_limits(3))
        .unwrap();
    compare_items(&got, &expected);
    assert_same_memo(&serial, &parallel);
    for item in &got.items[..2] {
        assert!(
            item.as_ref()
                .unwrap_err()
                .to_string()
                .contains("scalar seed index map is invalid")
        );
    }
    assert!(
        got.items[2]
            .as_ref()
            .unwrap()
            .reduction
            .unresolved
            .is_empty()
    );
    assert_eq!(got.stats.native_calls_returning_error, 2);
    assert_eq!(got.stats.prefetched_native_calls, 2);
    assert_eq!(got.stats.fallback_native_calls, 1);
    assert!(!got.stats.actual_application_count_complete);
    assert_eq!(got.stats.actual_known_rule_applications, 0);
}
