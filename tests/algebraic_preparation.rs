use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicSystem, RootSeed, SquareRoot};
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem};
use symbolica_amflow::{Error, FlowOptions, Precision, RunContext};

fn hierarchy() -> (AlgebraicSystem, Symbol) {
    let x = symbol!("prepared_algebraic_tests::x");
    let r = symbol!("prepared_algebraic_tests::r");
    (
        AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![
                    vec![vec![Atom::new()]],
                    vec![vec![Atom::num(1) / Atom::var(r)]],
                    vec![vec![Atom::new()]],
                    vec![vec![Atom::new()]],
                ],
            },
            roots: vec![SquareRoot {
                symbol: r,
                radicand: Atom::num(1) + Atom::var(x),
            }],
            nonzero_conditions: Vec::new(),
        },
        r,
    )
}

#[test]
fn immutable_preparation_preserves_both_sheets_and_fresh_epsilon_profiles() {
    let (mut system, root) = hierarchy();
    let prepared = system.prepare().unwrap();
    // Mutating the caller-owned object cannot change an existing preparation.
    system.system.matrices[1][0][0] = Atom::num(7);
    system.roots[0].radicand = Atom::num(9);
    system
        .nonzero_conditions
        .push(Atom::var(system.system.variable));
    assert_eq!(
        prepared.source().system.matrices[1][0][0],
        Atom::num(1) / Atom::var(root)
    );
    assert!(prepared.source().nonzero_conditions.is_empty());
    for (working, digits, order) in [(60, 35, 64), (100, 75, 112)] {
        let p = Precision::decimal(working).unwrap();
        let compiled = prepared.compile(p).unwrap();
        let point = p.rational(&Rational::from((1, 16)));
        let positive = symbolica_amflow::ComplexFloat::new(
            p.rational(&Rational::from((17, 16))).re.sqrt(),
            p.real(0),
        );
        for sign in [1, -1] {
            let seeds = BTreeMap::from([(
                root,
                if sign == 1 {
                    RootSeed::Principal
                } else {
                    RootSeed::Opposite
                },
            )]);
            let result = compiled
                .transport(
                    &EpsilonBoundary {
                        point: p.zero(),
                        leading: -2,
                        coefficients: vec![
                            vec![p.i(1)],
                            vec![p.zero()],
                            vec![p.zero()],
                            vec![p.zero()],
                        ],
                    },
                    std::slice::from_ref(&point),
                    &seeds,
                    &FlowOptions {
                        digits,
                        guard_digits: working - digits,
                        series_order: order,
                        ..Default::default()
                    },
                    &RunContext::default(),
                    true,
                )
                .unwrap();
            assert_eq!(result.solution.leading, -2);
            let primitive = p.scale(&p.sub(&positive, &p.i(1)), 2 * sign, 1);
            let mut expected = p.i(1);
            for n in 0..4 {
                if n > 0 {
                    expected = p.scale(&p.mul(&expected, &primitive), 1, n as i64);
                }
                assert!(
                    p.norm(&p.sub(&result.solution.coefficients[n][0], &expected))
                        < p.tolerance(digits)
                );
            }
            assert!(
                p.norm(&p.sub(&result.branches.roots[0].1, &p.scale(&positive, sign, 1)))
                    < p.tolerance(digits)
            );
            assert_eq!(result.solution.diagnostics.working_bits, p.bits);
            assert_eq!(result.branches.roots[0].1.re.as_raw().prec(), p.bits);
        }
    }
}

#[test]
fn high_precision_root_values_are_recomputed_from_exact_source() {
    let (mut system, root) = hierarchy();
    system.roots[0].radicand = Atom::num((1, 3)) + Atom::var(system.system.variable);
    let prepared = system.prepare().unwrap();
    let low = Precision { bits: 80 };
    let high = Precision { bits: 400 };
    let seeds = BTreeMap::from([(root, RootSeed::Principal)]);
    let low_value = prepared
        .compile(low)
        .unwrap()
        .branch_state_at(&low.zero(), &seeds)
        .unwrap()
        .roots[0]
        .1
        .clone();
    let high_value = prepared
        .compile(high)
        .unwrap()
        .branch_state_at(&high.zero(), &seeds)
        .unwrap()
        .roots[0]
        .1
        .clone();
    let expected = symbolica_amflow::ComplexFloat::new(
        high.rational(&Rational::from((1, 3))).re.sqrt(),
        high.real(0),
    );
    assert!(high.norm(&high.sub(&high_value, &expected)) < high.tolerance(110));
    assert!(high.norm(&high.sub(&high.round(&low_value), &expected)) > high.tolerance(30));
    assert_eq!(high_value.re.as_raw().prec(), high.bits);
}

#[test]
fn original_domain_holes_survive_repeated_preparation_compilation() {
    let (mut system, root) = hierarchy();
    let x = system.system.variable;
    system.nonzero_conditions = vec![Atom::var(x) - Atom::num((1, 32))];
    let prepared = system.prepare().unwrap();
    system.nonzero_conditions.clear();
    for working in [40, 90] {
        let p = Precision::decimal(working).unwrap();
        let compiled = prepared.compile(p).unwrap();
        let hole = p.rational(&Rational::from((1, 32)));
        assert!(
            compiled
                .singularities()
                .iter()
                .any(|point| p.close(point, &hole, working - 5))
        );
        assert!(
            compiled
                .branch_state_at(&hole, &BTreeMap::from([(root, RootSeed::Principal)]))
                .is_err()
        );
    }
}

#[test]
fn preparation_and_fresh_compilation_preserve_typed_validation_and_cancellation() {
    let (system, _) = hierarchy();
    let prepared = system.prepare().unwrap();
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        system.prepare_with_context(&context),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        prepared.compile_with_context(Precision { bits: 200 }, &context),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        prepared.compile(Precision { bits: 1 }),
        Err(Error::InvalidInput(_))
    ));
    let mut invalid = system;
    invalid.roots[0].radicand = Atom::new();
    assert!(matches!(invalid.prepare(), Err(Error::InvalidInput(_))));
}
