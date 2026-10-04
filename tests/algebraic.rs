use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicSystem, RootSeed, SquareRoot};
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem};
use symbolica_amflow::{
    BoundaryData, DifferentialSystem, Error, FlowOptions, Precision, Prescription, RunContext,
};

fn definitions() -> (Symbol, Symbol, Vec<SquareRoot>) {
    let x = symbol!("algebraic_tests::x");
    let r = symbol!("algebraic_tests::r");
    (
        x,
        r,
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(x),
        }],
    )
}
fn options(working: u32, order: usize) -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: working - 20,
        series_order: order,
        ..Default::default()
    }
}
fn winding(p: Precision) -> Vec<symbolica_amflow::ComplexFloat> {
    vec![
        p.complex(1, 1),
        p.complex(-1, 1),
        p.complex(-1, -1),
        p.complex(1, -1),
        p.i(1),
    ]
}
#[test]
fn ordinary_winding_tracks_sheets_and_refines() {
    let (x, r, roots) = definitions();
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::num(1) / Atom::var(r)]],
        },
        roots,
    );
    let mut previous = None;
    for (working, order) in [(50, 64), (70, 88)] {
        let p = Precision::decimal(working).unwrap();
        let compiled = system.compile(p).unwrap();
        let boundary = BoundaryData {
            point: p.i(1),
            values: vec![p.i(1)],
        };
        let seeds = BTreeMap::from([(r, RootSeed::Principal)]);
        let one = compiled
            .transport_ordinary(
                &boundary,
                &winding(p),
                &seeds,
                &options(working, order),
                &RunContext::default(),
                true,
            )
            .unwrap();
        let expected = p.exp(&p.i(-4));
        assert!(p.norm(&p.sub(&one.solution.coefficients[0][0], &expected)) < p.tolerance(35));
        assert_eq!(one.branches.point, p.i(1));
        assert_eq!(one.branches.roots[0].0, r);
        assert!(p.norm(&p.add(&one.branches.roots[0].1, &p.i(1))) < p.tolerance(35));
        assert!(one.solution.verified_digits.is_none());
        assert!(one.solution.diagnostics.rejected_steps > 0);
        assert_eq!(one.solution.segments.len(), one.solution.diagnostics.steps);
        if let Some(previous) = previous {
            assert!(p.norm(&p.sub(&one.solution.coefficients[0][0], &previous)) < p.tolerance(35));
        }
        previous = Some(one.solution.coefficients[0][0].clone());
        let two = winding(p).into_iter().chain(winding(p)).collect::<Vec<_>>();
        let result = compiled
            .transport_ordinary(
                &boundary,
                &two,
                &seeds,
                &options(working, order),
                &RunContext::default(),
                false,
            )
            .unwrap();
        assert!(p.norm(&p.sub(&result.solution.coefficients[0][0], &p.i(1))) < p.tolerance(35));
        assert!(p.norm(&p.sub(&result.branches.roots[0].1, &p.i(1))) < p.tolerance(35));
    }
}
#[test]
fn direct_epsilon_hierarchy_and_concurrent_sheets() {
    let (x, r, roots) = definitions();
    let p = Precision::decimal(50).unwrap();
    let system = AlgebraicSystem {
        system: EpsilonSystem {
            variable: x,
            matrices: vec![
                vec![vec![Atom::new()]],
                vec![vec![Atom::num(1) / Atom::var(r)]],
                vec![vec![Atom::new()]],
                vec![vec![Atom::new()]],
            ],
        },
        roots,
    };
    let compiled = system.compile(p).unwrap();
    let boundary = EpsilonBoundary {
        point: p.i(1),
        leading: -2,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
    };
    std::thread::scope(|scope| {
        let mut workers = Vec::new();
        for (seed, exponent) in [(RootSeed::Principal, -4), (RootSeed::Opposite, 4)] {
            let compiled = &compiled;
            let boundary = &boundary;
            workers.push(scope.spawn(move || {
                let result = compiled
                    .transport(
                        boundary,
                        &winding(p),
                        &BTreeMap::from([(r, seed)]),
                        &options(50, 64),
                        &RunContext::default(),
                        false,
                    )
                    .unwrap();
                let mut expected = p.i(1);
                for k in 0..4 {
                    assert!(
                        p.norm(&p.sub(&result.solution.coefficients[k][0], &expected))
                            < p.tolerance(32)
                    );
                    expected = p.scale(&expected, exponent, k as i64 + 1);
                }
                assert_eq!(result.solution.leading, -2);
            }));
        }
        for worker in workers {
            worker.join().unwrap();
        }
    });
}
#[test]
fn root_powers_normalize_exactly_and_regular_rational_systems_remain_valid() {
    let (x, r, roots) = definitions();
    let p = Precision::decimal(45).unwrap();
    let a = Atom::var(x);
    let b = Atom::var(r);
    // r^3/x - r + r^-1 = r^-1 exactly on either declared sheet.
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![b.clone().pow(3) / a - &b + b.clone().pow(-1)]],
        },
        roots,
    );
    let compiled = system.compile(p).unwrap();
    let result = compiled
        .transport_ordinary(
            &BoundaryData {
                point: p.i(1),
                values: vec![p.i(1)],
            },
            &[p.i(4)],
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &options(45, 64),
            &RunContext::default(),
            false,
        )
        .unwrap();
    assert!(p.norm(&p.sub(&result.solution.coefficients[0][0], &p.exp(&p.i(2)))) < p.tolerance(28));
    let reduced_denominator = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::num(1) / (Atom::var(r).pow(2) + 1)]],
        },
        definitions().2,
    )
    .compile(p)
    .unwrap();
    let result = reduced_denominator
        .transport_ordinary(
            &BoundaryData {
                point: p.i(1),
                values: vec![p.i(1)],
            },
            &[p.i(3)],
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &options(45, 64),
            &RunContext::default(),
            false,
        )
        .unwrap();
    assert!(p.norm(&p.sub(&result.solution.coefficients[0][0], &p.i(2))) < p.tolerance(28));
    let zero_denominator = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::num(1) / (Atom::var(r).pow(2) - Atom::var(x))]],
        },
        definitions().2,
    );
    assert!(matches!(
        zero_denominator.compile(p),
        Err(Error::InvalidInput(_))
    ));
    let rational = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::num(1)]],
        },
        vec![],
    )
    .compile(p)
    .unwrap();
    let result = rational
        .transport_ordinary(
            &BoundaryData {
                point: p.i(0),
                values: vec![p.i(1)],
            },
            &[p.i(1)],
            &BTreeMap::new(),
            &options(45, 64),
            &RunContext::default(),
            false,
        )
        .unwrap();
    assert!(p.norm(&p.sub(&result.solution.coefficients[0][0], &p.exp(&p.i(1)))) < p.tolerance(28));
}
#[test]
fn explicit_radicand_i0_and_invalid_registry_seeds() {
    let (x, r, roots) = definitions();
    let p = Precision::decimal(40).unwrap();
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::new()]],
        },
        roots.clone(),
    );
    let compiled = system.compile(p).unwrap();
    let boundary = BoundaryData {
        point: p.i(-1),
        values: vec![p.i(1)],
    };
    for (side, sign) in [(Prescription::PlusI0, 1), (Prescription::MinusI0, -1)] {
        let result = compiled
            .transport_ordinary(
                &boundary,
                &[],
                &BTreeMap::from([(r, RootSeed::I0(side))]),
                &options(40, 32),
                &RunContext::default(),
                false,
            )
            .unwrap();
        assert_eq!(result.branches.roots[0].1, p.complex(0, sign));
    }
    for value in [p.zero(), p.i(1), p.parse("inf", "0").unwrap()] {
        assert!(
            compiled
                .transport_ordinary(
                    &boundary,
                    &[],
                    &BTreeMap::from([(r, RootSeed::Value(value))]),
                    &options(40, 32),
                    &RunContext::default(),
                    false
                )
                .is_err()
        );
    }
    let hint = compiled
        .transport_ordinary(
            &boundary,
            &[],
            &BTreeMap::from([(r, RootSeed::Value(p.parse("0", "1.001").unwrap()))]),
            &options(40, 32),
            &RunContext::default(),
            false,
        )
        .unwrap();
    assert_eq!(hint.branches.roots[0].1, p.complex(0, 1));
    assert!(matches!(
        system.compile(Precision { bits: 1 }),
        Err(Error::InvalidInput(_))
    ));
    let mut invalid = system.clone();
    invalid.roots[0].symbol = x;
    assert!(invalid.validate().is_err());
    invalid = system.clone();
    invalid.roots[0].radicand = Atom::var(r);
    assert!(invalid.validate().is_err());
    invalid = system.clone();
    invalid.system.matrices[0][0][0] = Atom::num(1) / (Atom::var(r) + 1);
    assert!(matches!(invalid.compile(p), Err(Error::Unsupported(_))));
    invalid = system.clone();
    invalid.system.matrices[0][0][0] = parse!("sin(algebraic_tests::x)");
    assert!(invalid.compile(p).is_err());
    assert!(
        compiled
            .transport_ordinary(
                &boundary,
                &[],
                &BTreeMap::new(),
                &options(40, 32),
                &RunContext::default(),
                false
            )
            .is_err()
    );
}

#[test]
fn two_root_product_uses_independent_sheets_and_rational_radicands() {
    let x = symbol!("algebraic_product::x");
    let r = symbol!("algebraic_product::r");
    let s = symbol!("algebraic_product::s");
    let p = Precision::decimal(45).unwrap();
    let first = Atom::var(r);
    let second = Atom::var(s);
    let compiled = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![&first * &second + first.pow(2) * second.pow(2) - 1]],
        },
        vec![
            SquareRoot {
                symbol: r,
                radicand: Atom::var(x),
            },
            SquareRoot {
                symbol: s,
                radicand: Atom::num(1) / Atom::var(x),
            },
        ],
    )
    .compile(p)
    .unwrap();
    let boundary = BoundaryData {
        point: p.i(1),
        values: vec![p.i(1)],
    };
    for (seed, sign) in [(RootSeed::Principal, 1), (RootSeed::Opposite, -1)] {
        let result = compiled
            .transport_ordinary(
                &boundary,
                &[p.i(2)],
                &BTreeMap::from([(r, RootSeed::Principal), (s, seed)]),
                &options(45, 64),
                &RunContext::default(),
                false,
            )
            .unwrap();
        assert!(
            p.norm(&p.sub(&result.solution.coefficients[0][0], &p.exp(&p.i(sign))))
                < p.tolerance(28)
        );
        assert!(
            p.norm(&p.sub(
                &p.mul(&result.branches.roots[0].1, &result.branches.roots[1].1),
                &p.i(sign)
            )) < p.tolerance(35)
        );
    }
}
