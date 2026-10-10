//! Isolated cut-then-occupation-degree experimental-order invariants; no production source or proof edits.
use rustred::algebra::CoefficientContext;
use rustred::solver::guarded::{
    GuardedApplicationFailure, GuardedApplicationStatus, GuardedProgram, GuardedSource,
    GuardedSourceSystem, IndexDomain, IndexRole,
};
use rustred::solver::{Integral, IntegralOrder, Power, Term};
use std::cmp::Ordering;
use std::sync::Arc;

fn numeric(p: [i16; 4]) -> Integral<4> {
    Integral::numeric(p).unwrap()
}

fn main() {
    use IndexRole::{Occupation, Ordinary, RequiredCut};
    let roles = [RequiredCut, Ordinary, Occupation, Occupation];
    let order = IntegralOrder::new([true, true, false, false], [false; 4])
        .with_roles(roles)
        .unwrap();
    assert_eq!(order.deltas(), &[true, false, false, false]);
    // Physical sector and required-cut degree precede occupation degree;
    // occupation degree precedes the remaining physical complexity.
    assert_eq!(order.compare(&numeric([1, 1, 0, 0]), &numeric([1, -1, 0, 2])), Ordering::Less);
    assert_eq!(order.compare(&numeric([2, 1, 0, 0]), &numeric([1, 1, 0, 1])), Ordering::Less);
    assert_eq!(order.compare(&numeric([1, 10, 0, 0]), &numeric([1, 1, 0, 1])), Ordering::Greater);
    assert_eq!(order.compare(&numeric([1, 2, 0, 2]), &numeric([1, 1, 1, 1])), Ordering::Less);
    assert_eq!(order.compare(&numeric([1, 1, 2, 0]), &numeric([1, 1, 1, 1])), Ordering::Less);

    let mut points = Vec::new();
    for cut in 1..=2 {
        for ordinary in -2..=2 {
            for upper in 0..=2 {
                for lower in 0..=2 {
                    points.push(numeric([cut, ordinary, upper, lower]));
                }
            }
        }
    }
    let mut pairs = 0usize;
    let mut triples = 0usize;
    for a in &points {
        for b in &points {
            let ab = order.compare(a, b);
            assert_eq!(ab, order.compare(b, a).reverse());
            assert_eq!(ab == Ordering::Equal, a == b);
            pairs += 1;
            for c in &points {
                let bc = order.compare(b, c);
                if ab != Ordering::Greater && bc != Ordering::Greater {
                    assert_ne!(order.compare(a, c), Ordering::Greater);
                }
                triples += 1;
            }
        }
    }

    let mut affine_checks = 0usize;
    for positive in [false, true] {
        for occupation_sector in [false, true] {
            let mixed_order = IntegralOrder::new(
                [true, positive, occupation_sector, occupation_sector],
                [false; 4],
            )
            .with_roles(roles)
            .unwrap();
            for pattern in 0u8..16 {
                let mut keys = Vec::new();
                for a in -1i16..=1 {
                    for b in -1i16..=1 {
                        for c in -1i16..=1 {
                            for d in -1i16..=1 {
                                let values = [a, b, c, d];
                                keys.push(Integral::new(std::array::from_fn(|axis| {
                                    let symbolic = pattern & (1 << axis) != 0;
                                    let value = if symbolic {
                                        values[axis]
                                    } else {
                                        values[axis] + match axis {
                                            0 => 2,
                                            1 if positive => 2,
                                            1 => -2,
                                            _ => 1,
                                        }
                                    };
                                    Power::new(symbolic, value).unwrap()
                                })));
                            }
                        }
                    }
                }
                for base in [[3, if positive { 3 } else { -3 }, 2, 4], [7, if positive { 7 } else { -7 }, 4, 2]] {
                    let concrete = |key: &Integral<4>| numeric(std::array::from_fn(|axis| {
                        key[axis].value() + if key[axis].is_symbolic() { base[axis] } else { 0 }
                    }));
                    for left in &keys {
                        for right in &keys {
                            assert_eq!(
                                mixed_order.compare(left, right),
                                mixed_order.compare(&concrete(left), &concrete(right)),
                                "symbolic/numeric mismatch: pattern={pattern}, base={base:?}, left={left:?}, right={right:?}"
                            );
                            affine_checks += 1;
                        }
                    }
                }
            }
        }
    }

    let ordinary = IntegralOrder::new([true, false, true, false], [false; 4]);
    let typed = ordinary.clone().with_roles([Ordinary; 4]).unwrap();
    let mut no_occupation_checks = 0;
    for left in &points {
        for right in &points {
            assert_eq!(ordinary.compare(left, right), typed.compare(left, right));
            no_occupation_checks += 1;
        }
    }

    // No rules are installed. These results test only native role semantics,
    // so a theta H0 must remain an uncovered bulk integral, never a cut zero.
    let coefficients = CoefficientContext::try_new(["inv_a", "inv_b", "inv_h", "inv_j"]).unwrap();
    let source = GuardedSource::new(
        "formal-invariant-fixture",
        vec![Term { integral: Integral::symbolic([0; 4]).unwrap(), coefficient: coefficients.one().numerator }],
        IndexDomain::for_roles(&roles),
    );
    let context = Arc::new(GuardedSourceSystem::new(
        "isolated occupation-order role invariant", roles, [0, 1, 2, 3], vec![source],
    ).unwrap());
    let program = GuardedProgram::new(context, vec![], []).unwrap();
    let mut zero_checks = 0;
    for cut in [-1, 0] {
        for upper in 0..=2 {
            for lower in 0..=2 {
                assert!(matches!(program.apply(&[cut, 1, upper, lower]).unwrap().status, GuardedApplicationStatus::Zero));
                zero_checks += 1;
            }
        }
    }
    assert!(matches!(program.apply(&[1, 1, 0, 0]).unwrap().status,
        GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)));
    assert!(matches!(program.apply(&[1, 1, 1, 0]).unwrap().status,
        GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)));
    for cut in [0, 1] {
        assert!(matches!(program.apply(&[cut, 1, -1, 0]).unwrap().status,
            GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::InvalidOccupation { .. })));
    }
    let result = serde_json::json!({
        "status":"pass",
        "scope":"isolated native cut-then-occupation-degree comparator and role invariants; no closure or amplitude claim",
        "explicit_priority_witnesses":5,
        "numeric_points":points.len(),
        "antisymmetry_equality_pairs":pairs,
        "transitivity_triples":triples,
        "symbolic_common_base_comparisons":affine_checks,
        "no_occupation_invariance_pairs":no_occupation_checks,
        "required_cut_zero_checks":zero_checks,
        "bulk_surface_not_zero_checks":2,
        "invalid_occupation_checks":2,
    });
    std::fs::write(std::env::args().nth(1).unwrap(), serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    println!("{result}");
}
