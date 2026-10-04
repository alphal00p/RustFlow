use symbolica::prelude::*;
use symbolica_amflow::frobenius::FrobeniusBasis;
use symbolica_amflow::*;

fn logarithmic_basis() -> FrobeniusBasis {
    DifferentialSystem {
        variable: symbol!("partial_x"),
        matrix: vec![
            vec![Atom::new(), parse!("1/partial_x")],
            vec![Atom::new(), Atom::new()],
        ],
    }
    .frobenius(Precision::decimal(60).unwrap(), &Default::default(), 4)
    .unwrap()
}

fn constraint(component: usize, power: Atom, log_power: usize, value: i64) -> AsymptoticConstraint {
    AsymptoticConstraint {
        component,
        power,
        log_power,
        coefficient: Precision::decimal(60).unwrap().i(value),
    }
}

fn evaluate(
    basis: &FrobeniusBasis,
    constants: &[ComplexFloat],
    point: &ComplexFloat,
) -> Vec<ComplexFloat> {
    let p = basis.precision;
    basis
        .evaluate(point, &Default::default())
        .unwrap()
        .iter()
        .map(|row| {
            row.iter()
                .zip(constants)
                .fold(p.zero(), |sum, (a, c)| p.add(&sum, &p.mul(a, c)))
        })
        .collect()
}

#[test]
fn partial_log_coefficients_determine_unlisted_components() {
    let basis = logarithmic_basis();
    let p = basis.precision;
    let constants = basis
        .match_constraints(&[
            constraint(0, Atom::new(), 0, 2),
            constraint(0, Atom::new(), 1, 3),
        ])
        .unwrap();
    let point = p.rational(&Rational::from((1, 4)));
    let actual = evaluate(&basis, &constants, &point);
    assert!(p.close(
        &actual[0],
        &p.add(&p.i(2), &p.scale(&p.log(&point), 3, 1)),
        45
    ));
    assert!(p.close(&actual[1], &p.i(3), 45));
}

#[test]
fn partial_constraints_report_missing_constants_and_inconsistency() {
    let basis = logarithmic_basis();
    assert!(matches!(
        basis.match_constraints(&[constraint(0, Atom::new(), 0, 2)]),
        Err(Error::IncompleteReduction(_))
    ));
    assert!(matches!(
        basis.match_constraints(&[
            constraint(0, Atom::new(), 0, 2),
            constraint(0, Atom::new(), 1, 3),
            constraint(0, Atom::new(), 0, 4),
        ]),
        Err(Error::Numerical(_))
    ));
}

#[test]
fn fractional_powers_remain_distinct_exact_sectors() {
    let p = Precision::decimal(60).unwrap();
    let basis = DifferentialSystem {
        variable: symbol!("fractional_x"),
        matrix: vec![
            vec![parse!("1/(2*fractional_x)"), Atom::new()],
            vec![Atom::new(), parse!("1/(3*fractional_x)")],
        ],
    }
    .frobenius(p, &Default::default(), 4)
    .unwrap();
    let constants = basis
        .match_constraints(&[
            constraint(0, Atom::num(Rational::from((1, 2))), 0, 2),
            constraint(1, Atom::num(Rational::from((1, 3))), 0, 3),
        ])
        .unwrap();
    let actual = evaluate(&basis, &constants, &p.rational(&Rational::from((1, 64))));
    assert!(p.close(&actual[0], &p.rational(&Rational::from((1, 4))), 45));
    assert!(p.close(&actual[1], &p.rational(&Rational::from((3, 4))), 45));
}

#[test]
fn invalid_constraints_and_uncomputed_series_are_rejected() {
    let basis = logarithmic_basis();
    let mut nonfinite = constraint(0, Atom::new(), 0, 0);
    nonfinite.coefficient = basis.precision.parse("NaN", "0").unwrap();
    assert!(matches!(
        basis.match_constraints(&[nonfinite]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        basis.match_constraints(&[constraint(2, Atom::new(), 0, 1)]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        basis.match_constraints(&[]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        basis.match_constraints(&[constraint(0, Atom::num(50), 0, 1)]),
        Err(Error::Accuracy(_))
    ));
    assert!(matches!(
        basis.match_constraints(&[constraint(
            0,
            parse!("1000000000000000000000000000000000000"),
            0,
            0
        )]),
        Err(Error::Accuracy(_))
    ));
}

#[test]
fn integer_resonance_can_be_fixed_by_one_components_logarithm() {
    let p = Precision::decimal(60).unwrap();
    let basis = DifferentialSystem {
        variable: symbol!("resonant_partial_x"),
        matrix: vec![
            vec![Atom::new(), Atom::new()],
            vec![Atom::num(1), parse!("1/resonant_partial_x")],
        ],
    }
    .frobenius(p, &Default::default(), 4)
    .unwrap();
    let constants = basis
        .match_constraints(&[
            constraint(1, Atom::num(1), 0, 3),
            constraint(1, Atom::num(1), 1, 2),
        ])
        .unwrap();
    let point = p.rational(&Rational::from((1, 4)));
    let actual = evaluate(&basis, &constants, &point);
    assert!(p.close(&actual[0], &p.i(2), 45));
    assert!(p.close(
        &actual[1],
        &p.mul(&point, &p.add(&p.i(3), &p.scale(&p.log(&point), 2, 1))),
        45
    ));
}

#[test]
fn symbolic_regulator_powers_are_matched_before_specialization() {
    use symbolica_amflow::frobenius::FrobeniusColumn;
    let p = Precision::decimal(60).unwrap();
    let epsilon = parse!("partial_epsilon");
    let basis = FrobeniusBasis {
        precision: p,
        columns: vec![
            FrobeniusColumn {
                exponent: epsilon.clone(),
                coefficients: vec![vec![vec![p.i(1), p.zero()]]],
            },
            FrobeniusColumn {
                exponent: Atom::num(2) * &epsilon,
                coefficients: vec![vec![vec![p.zero(), p.i(1)]]],
            },
        ],
    };
    let constants = basis
        .match_constraints(&[
            constraint(0, epsilon.clone(), 0, 2),
            constraint(1, Atom::num(2) * &epsilon, 0, 3),
        ])
        .unwrap();
    assert!(p.close(&constants[0], &p.i(2), 45));
    assert!(p.close(&constants[1], &p.i(3), 45));
    // Even when both powers evaluate to zero at epsilon=0, the formal sectors
    // are not merged while extracting the coefficient constraints.
    assert!(matches!(
        basis.match_constraints(&[
            constraint(0, Atom::new(), 0, 2),
            constraint(1, Atom::new(), 0, 3),
        ]),
        Err(Error::IncompleteReduction(_))
    ));
}
