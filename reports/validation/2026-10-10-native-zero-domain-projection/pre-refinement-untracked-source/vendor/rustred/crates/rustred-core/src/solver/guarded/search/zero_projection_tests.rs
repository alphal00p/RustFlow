use super::*;
use crate::algebra::CoefficientContext;
use crate::solver::guarded::{GuardedApplicationStatus, GuardedProgram, GuardedSource};
use crate::solver::{Integral, SeedSource};
use std::sync::Arc;

fn point<const N: usize>(p: [i64; N]) -> IndexDomain<N> {
    IndexDomain::new(p.map(IndexBounds::fixed)).unwrap()
}

fn conditional_source(with_zero: bool) -> (CoefficientContext, Arc<GuardedSourceSystem<1>>) {
    let c = CoefficientContext::new(["zp_n", "zp_x", "zp_y"]);
    let source = GuardedSource::new(
        "conditional-recurrence",
        vec![
            Term {
                integral: Integral::symbolic([0]).unwrap(),
                coefficient: c.parameter("zp_x").unwrap().numerator,
            },
            Term {
                integral: Integral::symbolic([-1]).unwrap(),
                coefficient: -c.one().numerator,
            },
        ],
        IndexDomain::new([IndexBounds::new(Some(2), None).unwrap()]).unwrap(),
    )
    .with_nonzero_conditions(vec![c.parameter("zp_y").unwrap().numerator]);
    let system = GuardedSourceSystem::new(
        "uniform-measure-zero-projection",
        [IndexRole::Ordinary],
        [0],
        vec![source],
    )
    .unwrap()
    .with_zero_domains(if with_zero { vec![point([1])] } else { vec![] })
    .unwrap();
    (c, Arc::new(system))
}

#[test]
fn finite_point_zero_replays_conditions_and_binds_zero_context() {
    let (c, system) = conditional_source(true);
    let found = system.direct_zero_rules_at_points(&[[2]], 2).unwrap();
    assert_eq!(found.solution.rules.len(), 1);
    assert_eq!(found.completed_points, [[2]]);
    let rule = &found.solution.rules[0];
    assert!(rule.candidate().rhs.is_empty());
    assert!(
        rule.nonzero_conditions()
            .contains(&c.parameter("zp_x").unwrap().numerator)
    );
    assert!(
        rule.nonzero_conditions()
            .contains(&c.parameter("zp_y").unwrap().numerator)
    );
    system.replay_rule(rule).unwrap();
    let (_, missing_zero) = conditional_source(false);
    assert!(missing_zero.replay_rule(rule).is_err());
    let program = GuardedProgram::new(system.clone(), found.solution.rules, []).unwrap();
    let bytes = program.encode_native(Default::default()).unwrap();
    let decoded = GuardedProgram::decode_generated(&bytes, system, Default::default()).unwrap();
    let applied = decoded.apply(&[2]).unwrap();
    assert!(matches!(
        applied.status,
        GuardedApplicationStatus::Applied { .. }
    ));
    assert!(applied.terms.is_empty());
    assert_eq!(applied.nonzero_conditions.len(), 2);
    assert!(GuardedProgram::decode_generated(&bytes, missing_zero, Default::default()).is_err());
}

#[test]
fn partial_zero_box_does_not_erase_a_symbolic_source_column() {
    let (_, system) = conditional_source(true);
    let order = IntegralOrder::new([true], [false])
        .with_roles([IndexRole::Ordinary])
        .unwrap();
    let seed = Seed {
        integral: Integral::symbolic([0]).unwrap(),
        shifts: [0],
    };
    let partial = GuardedSearchScope::new(
        &system,
        IndexDomain::new([IndexBounds::new(Some(2), Some(4)).unwrap()]).unwrap(),
    );
    let row = partial
        .instantiate(0, &system.system.rows()[0], &seed, &order)
        .unwrap()
        .unwrap();
    // I(n-1) is zero only at n=2; the image [1,3] is not contained in {1}.
    assert_eq!(row.len(), 2);
    let complete = GuardedSearchScope::new(&system, point([2]));
    assert_eq!(
        complete
            .instantiate(0, &system.system.rows()[0], &seed, &order)
            .unwrap()
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn occupation_bulk_survives_and_invalid_images_are_checked_before_projection() {
    let c = CoefficientContext::new(["zp_a", "zp_h"]);
    let roles = [IndexRole::Ordinary, IndexRole::Occupation];
    let rows = vec![
        GuardedSource::new(
            "valid-bulk-and-surface",
            vec![Term {
                integral: Integral::symbolic([0, 0]).unwrap(),
                coefficient: c.one().numerator,
            }],
            IndexDomain::for_roles(&roles),
        ),
        GuardedSource::new(
            "invalid-image-must-not-be-hidden",
            vec![
                Term {
                    integral: Integral::symbolic([0, 1]).unwrap(),
                    coefficient: c.one().numerator,
                },
                Term {
                    integral: Integral::symbolic([0, -1]).unwrap(),
                    coefficient: c.one().numerator,
                },
            ],
            IndexDomain::for_roles(&roles),
        ),
    ];
    let system = GuardedSourceSystem::new("occupation-projection", roles, [0, 1], rows)
        .unwrap()
        .with_zero_domains(vec![
            IndexDomain::new([
                IndexBounds::unbounded(),
                IndexBounds::new(Some(1), None).unwrap(),
            ])
            .unwrap(),
        ])
        .unwrap();
    let order = IntegralOrder::new([true, false], [false; 2])
        .with_roles(roles)
        .unwrap();
    let scope = GuardedSearchScope::new(&system, point([2, 0]));
    let seed = Seed {
        integral: Integral::numeric([2, 0]).unwrap(),
        shifts: [0; 2],
    };
    assert_eq!(
        scope
            .instantiate(0, &system.system.rows()[0], &seed, &order)
            .unwrap()
            .unwrap()
            .len(),
        1
    );
    assert!(
        scope
            .instantiate(1, &system.system.rows()[1], &seed, &order)
            .unwrap()
            .is_none()
    );
}

#[test]
fn recentering_keeps_the_projected_zero_discovery_domain_guard() {
    let c = CoefficientContext::new(["zp_shift_n", "zp_shift_p"]);
    let roles = [IndexRole::Ordinary; 2];
    let source = GuardedSource::new(
        "shifted-measure-zero",
        vec![
            Term {
                integral: Integral::symbolic([0, 0]).unwrap(),
                coefficient: c.one().numerator,
            },
            Term {
                integral: Integral::symbolic([-1, -1]).unwrap(),
                coefficient: -c.one().numerator,
            },
        ],
        IndexDomain::unrestricted(),
    );
    let discovery = IndexDomain::new([
        IndexBounds::new(Some(2), Some(4)).unwrap(),
        IndexBounds::fixed(0),
    ])
    .unwrap();
    let system = GuardedSourceSystem::new("recentered-zero", roles, [0, 1], vec![source])
        .unwrap()
        .with_zero_domains(vec![
            IndexDomain::new([
                IndexBounds::new(Some(2), Some(4)).unwrap(),
                IndexBounds::fixed(-1),
            ])
            .unwrap(),
        ])
        .unwrap();
    let case = CoordinateCase::new([None, Some(0)]).unwrap();
    let candidate = RuleCandidate {
        target: case.integral(),
        case: case.clone().into(),
        rhs: vec![],
        sources: vec![SeedSource {
            basis_row: 0,
            seed: Seed {
                integral: case.integral().shifted([1, 0]).unwrap(),
                shifts: [1, 0],
            },
        }],
        stats: Default::default(),
    };
    let order = IntegralOrder::new([true, false], [false; 2])
        .with_roles(roles)
        .unwrap();
    let mut rule = system
        .seal_candidate(candidate, order, discovery.clone())
        .unwrap();
    assert_eq!(
        rule.domain.bounds()[0],
        IndexBounds::new(Some(3), Some(4)).unwrap()
    );
    system.replay_rule(&rule).unwrap();
    rule.domain = discovery;
    assert!(system.replay_rule(&rule).is_err());
}

#[test]
fn required_cut_zeros_and_explicit_measure_zeros_coexist() {
    let c = CoefficientContext::new(["zp_plain", "zp_cut", "zp_bulk"]);
    let roles = [
        IndexRole::Ordinary,
        IndexRole::RequiredCut,
        IndexRole::Occupation,
    ];
    let source = GuardedSource::new(
        "both-kinds-of-zero",
        [[0, 0, 0], [-1, 0, 0], [0, -1, 0]]
            .into_iter()
            .map(|shift| Term {
                integral: Integral::symbolic(shift).unwrap(),
                coefficient: c.one().numerator,
            })
            .collect(),
        IndexDomain::for_roles(&roles),
    );
    let system = Arc::new(
        GuardedSourceSystem::new("both-zero-evidence", roles, [0, 1, 2], vec![source])
            .unwrap()
            .with_zero_domains(vec![point([1, 1, 0])])
            .unwrap(),
    );
    let found = system.direct_zero_rules_at_points(&[[2, 1, 0]], 3).unwrap();
    assert_eq!(found.solution.rules.len(), 1);
    let program = GuardedProgram::new(system, found.solution.rules, []).unwrap();
    assert!(program.apply(&[2, 1, 0]).unwrap().terms.is_empty());
}
