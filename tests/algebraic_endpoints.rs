use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicSystem, RootSeed, SquareRoot};
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem};
use symbolica_amflow::*;

fn atom(s: &str) -> Atom {
    Atom::parse(s, "algebraic_endpoint_tests", Default::default()).unwrap()
}
fn ordinary(matrix: &[&[&str]], radicands: &[(&str, &str)]) -> AlgebraicSystem {
    AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: symbol!("algebraic_endpoint_tests::x"),
            matrix: matrix
                .iter()
                .map(|row| row.iter().map(|s| atom(s)).collect())
                .collect(),
        },
        radicands
            .iter()
            .map(|(name, expression)| SquareRoot {
                symbol: atom(name).as_view().as_var_view().unwrap().get_symbol(),
                radicand: atom(expression),
            })
            .collect(),
    )
}
fn seeds() -> BTreeMap<Symbol, RootSeed> {
    BTreeMap::from([(symbol!("algebraic_endpoint_tests::r"), RootSeed::Principal)])
}
fn boundary(p: Precision, values: Vec<Vec<ComplexFloat>>) -> EpsilonBoundary {
    EpsilonBoundary {
        point: p.rational(&Rational::from((1, 16))),
        leading: 0,
        coefficients: values,
    }
}

#[test]
fn puiseux_endpoint_matches_both_sheets_and_refines() {
    let system = ordinary(&[&["1/r"]], &[("r", "x")]);
    let prepared = system
        .prepare_frobenius(16, &RunContext::default())
        .unwrap();
    assert_eq!(
        prepared.lifted_system().system().matrix,
        vec![
            vec![atom("0"), atom("1/x")],
            vec![atom("1"), atom("1/(2*x)")]
        ]
    );
    let mut previous = None;
    for (digits, order) in [(55, 40), (80, 64)] {
        let p = Precision::decimal(digits).unwrap();
        for sign in [1, -1] {
            let choices = BTreeMap::from([(
                symbol!("algebraic_endpoint_tests::r"),
                if sign == 1 {
                    RootSeed::Principal
                } else {
                    RootSeed::Opposite
                },
            )]);
            let start = boundary(p, vec![vec![p.exp(&p.scale(&p.i(sign), 1, 2))]]);
            let expansion = prepared
                .match_boundary(&start, &choices, p, order, 0, &RunContext::default())
                .unwrap();
            let limit = expansion.coefficient_limits().unwrap();
            assert!(p.close(&limit[0][0], &p.i(1), 45));
            let point = p.rational(&Rational::from((1, 9)));
            for winding in [-1, 0, 1, 2] {
                let expected_sign = if winding % 2 == 0 { sign } else { -sign };
                let actual = expansion.evaluate(&point, winding).unwrap();
                assert!(p.close(
                    &actual.coefficients[0][0],
                    &p.exp(&p.scale(&p.i(expected_sign), 2, 3)),
                    45
                ));
            }
            if sign == 1 {
                if let Some(previous) = &previous {
                    assert!(p.close(&limit[0][0], previous, 45));
                }
                previous = Some(limit[0][0].clone());
            }
        }
    }
}

#[test]
fn resonant_coupled_blocks_produce_half_power_logarithms() {
    // y = [sqrt(x) log(x), 1, 2 sqrt(x) (log(x)-2)]. The rational
    // restriction has resonant sectors, despite the half-integral physical
    // powers. All physical components have a finite limit.
    let system = ordinary(
        &[
            &["1/(2*x)", "1/r", "0"],
            &["0", "0", "0"],
            &["1/x", "0", "0"],
        ],
        &[("r", "x")],
    );
    let prepared = system
        .prepare_frobenius(16, &RunContext::default())
        .unwrap();
    let p = Precision::decimal(70).unwrap();
    let x = p.rational(&Rational::from((1, 16)));
    let log = p.log(&x);
    let start = boundary(
        p,
        vec![vec![
            p.scale(&log, 1, 4),
            p.i(1),
            p.scale(&p.sub(&log, &p.i(2)), 1, 2),
        ]],
    );
    let expansion = prepared
        .match_boundary(&start, &seeds(), p, 12, 0, &RunContext::default())
        .unwrap();
    assert!(
        expansion
            .basis()
            .columns
            .iter()
            .any(|column| column.coefficients.iter().any(|logs| logs.len() > 1))
    );
    for winding in [0, 1] {
        let point = p.rational(&Rational::from((1, 9)));
        let log = p.add(&p.log(&point), &p.scale(&p.log(&p.i(-1)), 2 * winding, 1));
        let sign = if winding == 0 { 1 } else { -1 };
        let expected = [
            p.scale(&log, sign, 3),
            p.i(1),
            p.scale(&p.sub(&log, &p.i(2)), 2 * sign, 3),
        ];
        let actual = expansion.evaluate(&point, winding as i32).unwrap();
        for (a, b) in actual.coefficients[0].iter().zip(&expected) {
            assert!(p.close(a, b, 55));
        }
    }
    let limit = expansion.coefficient_limits().unwrap();
    for (a, b) in limit[0].iter().zip([p.zero(), p.i(1), p.zero()]) {
        assert!(p.close(a, &b, 55));
    }
}

#[test]
fn epsilon_hierarchy_keeps_laurent_leading_power_and_endpoint() {
    let mut system = ordinary(&[&["0"]], &[("r", "x")]);
    system.system = EpsilonSystem {
        variable: system.system.variable,
        matrices: vec![
            vec![vec![atom("0")]],
            vec![vec![atom("1/r")]],
            vec![vec![atom("0")]],
        ],
    };
    let prepared = system
        .prepare_frobenius(16, &RunContext::default())
        .unwrap();
    let p = Precision::decimal(60).unwrap();
    let mut start = boundary(p, vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]]);
    start.leading = -2;
    let expansion = prepared
        .match_boundary(&start, &seeds(), p, 8, 0, &RunContext::default())
        .unwrap();
    let at = expansion
        .evaluate(&p.rational(&Rational::from((1, 9))), 0)
        .unwrap();
    assert_eq!(at.leading, -2);
    for (row, expected) in
        at.coefficients
            .iter()
            .zip([p.i(1), p.scale(&p.i(1), 1, 6), p.scale(&p.i(1), 1, 72)])
    {
        assert!(p.close(&row[0], &expected, 45));
    }
    let limit = expansion.coefficient_limits().unwrap();
    for (row, expected) in
        limit
            .iter()
            .zip([p.i(1), p.scale(&p.i(-1), 1, 2), p.scale(&p.i(1), 1, 8)])
    {
        assert!(p.close(&row[0], &expected, 45));
    }
}

#[test]
fn dependent_roots_retain_independent_sheet_constraints() {
    let system = ordinary(&[&["1/r+1/s"]], &[("r", "x"), ("s", "4*x")]);
    let prepared = system.prepare_frobenius(8, &RunContext::default()).unwrap();
    let p = Precision::decimal(60).unwrap();
    for (seed, exponent) in [(RootSeed::Principal, 3), (RootSeed::Opposite, 1)] {
        let mut choices = seeds();
        choices.insert(symbol!("algebraic_endpoint_tests::s"), seed);
        let start = boundary(p, vec![vec![p.exp(&p.scale(&p.i(exponent), 1, 4))]]);
        let expansion = prepared
            .match_boundary(&start, &choices, p, 44, 0, &RunContext::default())
            .unwrap();
        let actual = expansion
            .evaluate(&p.rational(&Rational::from((1, 9))), 0)
            .unwrap();
        assert!(p.close(
            &actual.coefficients[0][0],
            &p.exp(&p.scale(&p.i(exponent), 1, 3)),
            45
        ));
        assert!(p.close(&expansion.coefficient_limits().unwrap()[0][0], &p.i(1), 45));
    }
    let nonunit = ordinary(&[&["1/(r-s)"]], &[("r", "x"), ("s", "x")]);
    assert!(matches!(
        nonunit.rational_lift(8, &RunContext::default()),
        Err(Error::Unsupported(_) | Error::InvalidInput(_))
    ));
}

#[test]
fn root_poles_and_infinity_use_the_same_projection() {
    // y(x)=1+1/sqrt(x), normalized with exact boundary. At infinity the
    // auxiliary r*y diverges, but the physical y has limit one.
    let system = ordinary(&[&["-1/(2*x*(r+1))"]], &[("r", "x")]);
    let transformed = system
        .pullback_coordinate(system.system.variable, &atom("1/x"))
        .unwrap();
    let prepared = transformed
        .prepare_frobenius(8, &RunContext::default())
        .unwrap();
    let p = Precision::decimal(60).unwrap();
    let start = boundary(p, vec![vec![p.rational(&Rational::from((5, 4)))]]);
    let expansion = prepared
        .match_boundary(&start, &seeds(), p, 36, 0, &RunContext::default())
        .unwrap();
    assert!(p.close(&expansion.coefficient_limits().unwrap()[0][0], &p.i(1), 40));
    let actual = expansion
        .evaluate(&p.rational(&Rational::from((1, 9))), 0)
        .unwrap();
    assert!(p.close(
        &actual.coefficients[0][0],
        &p.rational(&Rational::from((4, 3))),
        40
    ));
}

#[test]
fn domain_cancellation_limits_and_cancellation_are_explicit() {
    let system = ordinary(&[&["(x^2-1)/(r*(x-1))"]], &[("r", "x")]);
    assert!(matches!(
        system.rational_lift(1, &RunContext::default()),
        Err(Error::Limit(_))
    ));
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        system.prepare_frobenius(8, &context),
        Err(Error::Cancelled)
    ));
    let p = Precision::decimal(40).unwrap();
    let shifted = system
        .pullback_coordinate(system.system.variable, &atom("x+2"))
        .unwrap();
    let lifted = shifted.rational_lift(8, &RunContext::default()).unwrap();
    let at_hole = EpsilonBoundary {
        point: p.i(-1),
        leading: 0,
        coefficients: vec![vec![p.i(1)]],
    };
    assert!(lifted.lift_boundary(&at_hole, &seeds(), p).is_err());
    let prepared_hole = system.prepare_frobenius(8, &RunContext::default()).unwrap();
    let expansion_hole = prepared_hole
        .match_boundary(
            &boundary(p, vec![vec![p.i(1)]]),
            &seeds(),
            p,
            12,
            0,
            &RunContext::default(),
        )
        .unwrap();
    assert!(matches!(
        expansion_hole.evaluate(&p.i(1), 0),
        Err(Error::InvalidInput(_))
    ));
    let divergent = ordinary(&[&["-1/(2*x)"]], &[("r", "x")]);
    let prepared = divergent
        .prepare_frobenius(8, &RunContext::default())
        .unwrap();
    let expansion = prepared
        .match_boundary(
            &boundary(p, vec![vec![p.i(4)]]),
            &seeds(),
            p,
            8,
            0,
            &RunContext::default(),
        )
        .unwrap();
    assert!(matches!(
        expansion.coefficient_limits(),
        Err(Error::Numerical(_))
    ));
}

#[test]
fn root_product_subgroup_avoids_unused_radicals_and_full_extension_growth() {
    let system = ordinary(
        &[&["1/(r*s)"]],
        &[("r", "x"), ("s", "1+x"), ("unused", "2+x")],
    );
    let prepared = system.prepare_frobenius(2, &RunContext::default()).unwrap();
    assert_eq!(prepared.lifted_system().system().matrix.len(), 2);
    let p = Precision::decimal(70).unwrap();
    let mut choices = seeds();
    choices.insert(symbol!("algebraic_endpoint_tests::s"), RootSeed::Principal);
    choices.insert(
        symbol!("algebraic_endpoint_tests::unused"),
        RootSeed::Principal,
    );
    let half = p.rational(&Rational::from((1, 2)));
    let start_value = p.powi(
        &p.add(
            &p.rational(&Rational::from((1, 4))),
            &p.pow(&p.rational(&Rational::from((17, 16))), &half),
        ),
        2,
    );
    let expansion = prepared
        .match_boundary(
            &boundary(p, vec![vec![start_value]]),
            &choices,
            p,
            64,
            0,
            &RunContext::default(),
        )
        .unwrap();
    assert!(p.close(&expansion.coefficient_limits().unwrap()[0][0], &p.i(1), 55));
    let point = p.rational(&Rational::from((1, 9)));
    let expected = p.powi(
        &p.add(
            &p.rational(&Rational::from((1, 3))),
            &p.pow(&p.rational(&Rational::from((10, 9))), &half),
        ),
        2,
    );
    assert!(p.close(
        &expansion.evaluate(&point, 0).unwrap().coefficients[0][0],
        &expected,
        55
    ));
}

#[test]
fn lower_lip_matching_agrees_with_regular_algebraic_continuation() {
    let system = ordinary(&[&["1/r"]], &[("r", "x")]);
    let p = Precision::decimal(65).unwrap();
    let choices = BTreeMap::from([(
        symbol!("algebraic_endpoint_tests::r"),
        RootSeed::I0(Prescription::MinusI0),
    )]);
    let start = EpsilonBoundary {
        point: p.rational(&Rational::from((-1, 16))),
        leading: 0,
        coefficients: vec![vec![p.exp(&p.scale(&p.complex(0, -1), 1, 2))]],
    };
    let prepared = system.prepare_frobenius(8, &RunContext::default()).unwrap();
    let expansion = prepared
        .match_boundary(&start, &choices, p, 52, -1, &RunContext::default())
        .unwrap();
    let point = p.rational(&Rational::from((1, 16)));
    let regular = system
        .compile(p)
        .unwrap()
        .transport(
            &start,
            &[p.scale(&p.complex(0, -1), 1, 16), point.clone()],
            &choices,
            &FlowOptions {
                digits: 35,
                series_order: 80,
                ..Default::default()
            },
            &RunContext::default(),
            false,
        )
        .unwrap();
    let actual = expansion.evaluate(&point, 0).unwrap();
    assert!(p.close(
        &actual.coefficients[0][0],
        &p.exp(&p.scale(&p.i(1), 1, 2)),
        50
    ));
    assert!(p.close(
        &actual.coefficients[0][0],
        &regular.solution.coefficients[0][0],
        33
    ));
}

#[test]
fn exact_complex_radicands_keep_their_branch_and_precision() {
    let mut system = ordinary(&[&["1/r"]], &[("r", "x")]);
    // The reserved exact imaginary coefficient is also accepted by the public
    // algebraic transport input. The rational owner must receive native i.
    system.roots[0].radicand = Atom::var(symbol!("symbolica_amflow::imaginary_unit")) * atom("2*x");
    let prepared = system.prepare_frobenius(8, &RunContext::default()).unwrap();
    let p = Precision::decimal(70).unwrap();
    // sqrt(2 i x)=(1+i)sqrt(x), integral dx/sqrt(2 i x)=(1-i)sqrt(x).
    let start = boundary(p, vec![vec![p.exp(&p.scale(&p.complex(1, -1), 1, 4))]]);
    let expansion = prepared
        .match_boundary(&start, &seeds(), p, 48, 0, &RunContext::default())
        .unwrap();
    let at = expansion
        .evaluate(&p.rational(&Rational::from((1, 9))), 0)
        .unwrap();
    assert!(p.close(
        &at.coefficients[0][0],
        &p.exp(&p.scale(&p.complex(1, -1), 1, 3)),
        55
    ));
    assert!(p.close(&expansion.coefficient_limits().unwrap()[0][0], &p.i(1), 55));
}

#[test]
fn lifted_boundaries_roundtrip_through_the_existing_rational_solver() {
    let system = ordinary(&[&["1/r"]], &[("r", "x")]);
    let lifted = system.rational_lift(2, &RunContext::default()).unwrap();
    assert_eq!(lifted.root_monomials(), vec![atom("1"), atom("r")]);
    let p = Precision::decimal(60).unwrap();
    let mut start = boundary(p, vec![vec![p.i(1)]]);
    start.leading = -3;
    let (start_lift, _) = lifted.lift_boundary(&start, &seeds(), p).unwrap();
    let target = p.rational(&Rational::from((1, 9)));
    let result = lifted
        .system()
        .compile(p, &Default::default())
        .unwrap()
        .transport(
            &start_lift,
            std::slice::from_ref(&target),
            &FlowOptions {
                digits: 30,
                series_order: 72,
                ..Default::default()
            },
            &RunContext::default(),
        )
        .unwrap();
    let projected = lifted
        .project_boundary(
            &BoundaryData {
                point: result.point,
                values: result.values,
            },
            start.leading,
        )
        .unwrap();
    assert_eq!(projected.leading, -3);
    assert_eq!(projected.point, target);
    assert!(p.close(
        &projected.coefficients[0][0],
        &p.exp(&p.rational(&Rational::from((1, 6)))),
        28
    ));
}

#[test]
fn expanded_epsilon_sectors_are_not_misidentified_as_dimensional_limits() {
    let mut system = ordinary(&[&["0"]], &[]);
    system.system.matrices.push(vec![vec![atom("1/x")]]);
    let prepared = system.prepare_frobenius(2, &RunContext::default()).unwrap();
    let p = Precision::decimal(60).unwrap();
    let start = boundary(p, vec![vec![p.i(1)], vec![p.zero()]]);
    let expansion = prepared
        .match_boundary(&start, &BTreeMap::new(), p, 8, 0, &RunContext::default())
        .unwrap();
    let at = expansion
        .evaluate(&p.rational(&Rational::from((1, 9))), 0)
        .unwrap();
    assert!(p.close(
        &at.coefficients[1][0],
        &p.log(&p.rational(&Rational::from((16, 9)))),
        45
    ));
    assert!(
        matches!(expansion.coefficient_limits(), Err(Error::Numerical(message)) if message.contains("divergence"))
    );
}
