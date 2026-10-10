use super::*;
use crate::algebra::CoefficientContext;
use crate::solver::guarded::{GuardedApplicationStatus, GuardedProgram, GuardedSource};
use crate::solver::{SearchOptions, Term};
use std::sync::Arc;

fn hidden_zero() -> Arc<GuardedSourceSystem<2>> {
    let c = CoefficientContext::try_new(["dz_a", "dz_b", "dz_x"]).unwrap();
    let a = c.parameter("dz_a").unwrap();
    let x = c.parameter("dz_x").unwrap();
    let source = |id, row| {
        GuardedSource::new(
            id,
            row,
            IndexDomain::new([
                IndexBounds::new(Some(1), None).unwrap(),
                IndexBounds::new(None, Some(0)).unwrap(),
            ])
            .unwrap(),
        )
    };
    Arc::new(
        GuardedSourceSystem::new(
            "direct-zero-first-recurrence",
            [IndexRole::Ordinary; 2],
            [0, 1],
            vec![
                source(
                    "earlier-recurrence",
                    vec![
                        Term {
                            integral: Integral::symbolic([0, 0]).unwrap(),
                            coefficient: c.one().numerator,
                        },
                        Term {
                            integral: Integral::symbolic([-1, 0]).unwrap(),
                            coefficient: (-x).numerator,
                        },
                    ],
                ),
                source(
                    "later-aligned-zero",
                    vec![Term {
                        integral: Integral::symbolic([1, -1]).unwrap(),
                        coefficient: (&c.integer(-2) * &a).numerator,
                    }],
                ),
            ],
        )
        .unwrap(),
    )
}

#[test]
fn direct_zero_point_replays_existing_row_hidden_by_first_recurrence() {
    let sources = hidden_zero();
    let target = [2, -1];
    let domain = IndexDomain::new(target.map(IndexBounds::fixed)).unwrap();
    let first = sources
        .solve_domains(
            vec![domain],
            SearchOptions {
                max_depth: Some(3),
                ..Default::default()
            },
            1,
        )
        .unwrap();
    let first = GuardedProgram::new(sources.clone(), first.rules, []).unwrap();
    assert!(!first.apply(&target).unwrap().terms.is_empty());

    let zero = sources.direct_zero_rules_at_points(&[target], 3).unwrap();
    assert_eq!(zero.attempted_rows, 3);
    assert_eq!(zero.completed_points, [target]);
    assert!(zero.solution.unresolved.is_empty());
    assert_eq!(zero.solution.rules.len(), 1);
    let rule = &zero.solution.rules[0];
    assert!(rule.candidate().rhs.is_empty());
    assert_eq!(rule.candidate().sources[0].basis_row, 1);
    assert_eq!(
        rule.candidate().sources[0].seed.integral,
        Integral::numeric([1, 0]).unwrap()
    );
    sources.replay_rule(rule).unwrap();
    let program = GuardedProgram::new(sources.clone(), zero.solution.rules, []).unwrap();
    let applied = program.apply(&target).unwrap();
    assert!(matches!(
        applied.status,
        GuardedApplicationStatus::Applied { .. }
    ));
    assert!(applied.terms.is_empty());
    let bytes = program.encode_native(Default::default()).unwrap();
    let replayed = GuardedProgram::decode_generated(&bytes, sources, Default::default()).unwrap();
    assert!(
        replayed
            .reduce(target, Default::default())
            .unwrap()
            .unresolved
            .is_empty()
    );
}

#[test]
fn direct_zero_budget_is_global_and_completed_points_are_distinct_from_gaps() {
    let sources = hidden_zero();
    let points = [[2, -1], [2, -1], [3, -1]];
    let result = sources.direct_zero_rules_at_points(&points, 3).unwrap();
    assert_eq!(result.attempted_rows, 3);
    assert_eq!(result.solution.rules.len(), 1);
    assert_eq!(result.completed_points, [[2, -1]]);
    assert_eq!(result.solution.unresolved.len(), 1);
    assert_eq!(
        result.solution.unresolved[0].reason,
        GuardedUnresolvedReason::DomainBudget
    );
    assert!(result.solution.unresolved[0].domain.contains(&[3, -1]));
    let zero = sources.direct_zero_rules_at_points(&points, 0).unwrap();
    assert_eq!(zero.attempted_rows, 0);
    assert!(zero.completed_points.is_empty());
    assert_eq!(zero.solution.unresolved.len(), 2);
}

fn occupation_source(minimum: i64, vanished_condition: bool) -> Arc<GuardedSourceSystem<2>> {
    let c = CoefficientContext::try_new(["dz_n", "dz_h", "dz_condition"]).unwrap();
    let condition = if vanished_condition {
        c.parameter("dz_n").unwrap() - c.one()
    } else {
        c.parameter("dz_condition").unwrap()
    };
    Arc::new(
        GuardedSourceSystem::new(
            "direct-zero-occupation",
            [IndexRole::Ordinary, IndexRole::Occupation],
            [0, 1],
            vec![
                GuardedSource::new(
                    "surface-shift",
                    vec![Term {
                        integral: Integral::symbolic([1, -1]).unwrap(),
                        coefficient: c.one().numerator,
                    }],
                    IndexDomain::new([
                        IndexBounds::new(Some(1), None).unwrap(),
                        IndexBounds::new(Some(minimum), None).unwrap(),
                    ])
                    .unwrap(),
                )
                .with_nonzero_conditions(vec![condition.numerator]),
            ],
        )
        .unwrap(),
    )
}

#[test]
fn direct_zero_occupation_guards_conditions_and_upfront_validation() {
    let sources = occupation_source(1, false);
    let target = [2, 0];
    let result = sources.direct_zero_rules_at_points(&[target], 1).unwrap();
    assert_eq!(result.solution.rules.len(), 1);
    assert_eq!(result.solution.rules[0].nonzero_conditions().len(), 1);
    let program = GuardedProgram::new(sources.clone(), result.solution.rules, []).unwrap();
    assert_eq!(program.apply(&target).unwrap().nonzero_conditions.len(), 1);
    // Invalid occupation must fail even if an earlier point exhausts all work.
    assert!(
        sources
            .direct_zero_rules_at_points(&[target, [2, -1]], 0)
            .is_err()
    );

    let guard = occupation_source(2, false)
        .direct_zero_rules_at_points(&[target], 1)
        .unwrap();
    assert!(guard.solution.rules.is_empty());
    assert_eq!(guard.completed_points, [target]);
    assert_eq!(
        guard.solution.unresolved[0].reason,
        GuardedUnresolvedReason::SearchExhausted
    );
    assert_eq!(
        guard.skipped_seeds[0].reason,
        GuardedDirectZeroSkipReason::SourceGuard
    );
    let vanished = occupation_source(1, true)
        .direct_zero_rules_at_points(&[target], 1)
        .unwrap();
    assert!(vanished.solution.rules.is_empty());
    assert_eq!(
        vanished.skipped_seeds[0].reason,
        GuardedDirectZeroSkipReason::RejectedProof
    );
}

#[test]
fn direct_zero_never_promotes_unsupported_seeds_or_existing_zero_evidence() {
    let c = CoefficientContext::try_new(["dz_boundary"]).unwrap();
    let sources = GuardedSourceSystem::new(
        "direct-zero-overflow",
        [IndexRole::Ordinary],
        [0],
        vec![GuardedSource::new(
            "negative-offset",
            vec![Term {
                integral: Integral::symbolic([-1]).unwrap(),
                coefficient: c.one().numerator,
            }],
            IndexDomain::unrestricted(),
        )],
    )
    .unwrap();
    let result = sources
        .direct_zero_rules_at_points(&[[i64::from(Power::MAX)]], 1)
        .unwrap();
    assert!(result.solution.rules.is_empty());
    assert_eq!(
        result.skipped_seeds[0].reason,
        GuardedDirectZeroSkipReason::UnsupportedPower
    );
    let result = sources
        .direct_zero_rules_at_points(&[[i64::MAX]], 1)
        .unwrap();
    assert_eq!(result.attempted_rows, 0);
    assert!(result.completed_points.is_empty());
    assert_eq!(
        result.solution.unresolved[0].reason,
        GuardedUnresolvedReason::UnsupportedPower
    );

    let cut = GuardedSourceSystem::new(
        "direct-zero-required-cut",
        [IndexRole::RequiredCut],
        [0],
        vec![GuardedSource::new(
            "cut-source",
            vec![Term {
                integral: Integral::symbolic([0]).unwrap(),
                coefficient: c.one().numerator,
            }],
            IndexDomain::unrestricted(),
        )],
    )
    .unwrap();
    let result = cut.direct_zero_rules_at_points(&[[0], [1]], 0).unwrap();
    assert_eq!(result.already_zero_points, [[0]]);
    assert_eq!(result.completed_points, [[0]]);
    assert!(result.solution.rules.is_empty());
    assert_eq!(result.solution.unresolved.len(), 1);
}

#[test]
fn direct_zero_original_source_replay_rejects_corrupted_seed() {
    let sources = hidden_zero();
    let mut result = sources.direct_zero_rules_at_points(&[[2, -1]], 3).unwrap();
    let mut rule = result.solution.rules.pop().unwrap();
    rule.candidate.sources[0].seed.integral = Integral::numeric([2, 0]).unwrap();
    assert!(sources.replay_rule(&rule).is_err());
}

#[test]
fn direct_zero_rejects_invalid_occupation_images_and_sector_crossings() {
    let c = CoefficientContext::try_new(["dz_image_n", "dz_image_h"]).unwrap();
    let sources = GuardedSourceSystem::new(
        "direct-zero-invalid-image",
        [IndexRole::Ordinary, IndexRole::Occupation],
        [0, 1],
        vec![GuardedSource::new(
            "two-images",
            vec![
                Term {
                    integral: Integral::symbolic([1, -1]).unwrap(),
                    coefficient: c.one().numerator,
                },
                Term {
                    integral: Integral::symbolic([1, -2]).unwrap(),
                    coefficient: c.one().numerator,
                },
            ],
            IndexDomain::new([
                IndexBounds::new(Some(1), None).unwrap(),
                IndexBounds::new(Some(0), None).unwrap(),
            ])
            .unwrap(),
        )],
    )
    .unwrap();
    let result = sources.direct_zero_rules_at_points(&[[2, 0]], 2).unwrap();
    assert!(result.solution.rules.is_empty());
    assert!(
        result
            .skipped_seeds
            .iter()
            .any(|s| s.reason == GuardedDirectZeroSkipReason::InvalidOccupationImage)
    );
    assert_eq!(result.completed_points, [[2, 0]]);
    let result = sources.direct_zero_rules_at_points(&[[1, 0]], 2).unwrap();
    assert!(result.solution.rules.is_empty());
    assert!(
        result
            .skipped_seeds
            .iter()
            .all(|s| s.reason == GuardedDirectZeroSkipReason::SectorCrossing)
    );
}
