use symbolica::prelude::*;
use symbolica_amflow::{
    AsymptoticConstraint, BoundaryData, ComplexFloat, DifferentialSystem, Error, FlowOptions,
    Precision, Result, RunContext,
};
fn evaluate(kind: &str, digits: u32, order: usize) -> Result<Vec<ComplexFloat>> {
    let requested = if kind == "weight20" {
        if digits >= 140 { 70 } else { 50 }
    } else {
        30
    };
    let p = Precision::decimal(digits)?;
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let parameters = Default::default();
    let letters: Vec<Atom> = match kind {
        "trailing" => vec![Atom::num(1), Atom::num(-10), Atom::new()],
        "complex" => vec![
            Atom::num(10),
            Atom::num(-10) + imaginary,
            Atom::num(Rational::from((-1, 2))),
            Atom::num(-50),
        ],
        "weight20" => (1..=20).map(Atom::num).collect(),
        _ => return Err(Error::InvalidInput("unknown profile".into())),
    };
    let dimension = letters.len() + 1;
    let parameter = symbol!("mpl_additional::t");
    let mut matrix = vec![vec![Atom::new(); dimension]; dimension];
    for (i, a) in letters.iter().enumerate() {
        matrix[i][i + 1] = Atom::one() / (Atom::var(parameter) - a);
    }
    let system = DifferentialSystem {
        variable: parameter,
        matrix,
    };
    let boundary = if kind == "trailing" {
        let boundary_order = digits as usize;
        let basis = system.frobenius(p, &parameters, boundary_order)?;
        let constants = basis.match_constraints(
            &(0..dimension)
                .map(|component| AsymptoticConstraint {
                    component,
                    power: Atom::new(),
                    log_power: 0,
                    coefficient: if component + 1 == dimension {
                        p.i(1)
                    } else {
                        p.zero()
                    },
                })
                .collect::<Vec<_>>(),
        )?;
        let point = p.rational(&Rational::from((1, 8)));
        let values = basis
            .evaluate(&point, &Default::default())?
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&constants)
                    .fold(p.zero(), |sum, (a, c)| p.add(&sum, &p.mul(a, c)))
            })
            .collect();
        BoundaryData { point, values }
    } else {
        let mut values = vec![p.zero(); dimension];
        values[dimension - 1] = p.i(1);
        BoundaryData {
            point: p.zero(),
            values,
        }
    };
    let compiled = system.compile(p, &parameters)?;
    let endpoint = if kind == "complex" {
        1
    } else if kind == "trailing" {
        4
    } else {
        21
    };
    let waypoints = if kind == "complex" {
        vec![p.i(endpoint)]
    } else {
        vec![
            p.parse("0", "-0.5")?,
            p.parse(&endpoint.to_string(), "-0.5")?,
            p.i(endpoint),
        ]
    };
    let options = FlowOptions {
        digits: requested,
        guard_digits: digits - requested,
        series_order: order,
        max_steps: 20000,
        ..Default::default()
    };
    let result = compiled.transport(&boundary, &waypoints, &options, &RunContext::default())?;
    Ok(result.values)
}

fn check_profile(
    case: &str,
    low_digits: u32,
    low_order: usize,
    high_digits: u32,
    high_order: usize,
) {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../fixtures/diffexp/mpl-additional.json")).unwrap();
    let low = evaluate(case, low_digits, low_order).unwrap();
    let high = evaluate(case, high_digits, high_order).unwrap();
    assert_eq!(low.len(), high.len());
    let p = Precision::decimal(high_digits).unwrap();
    for (a, b) in low.iter().zip(&high) {
        assert!(p.norm(&p.sub(a, b)) < p.tolerance(20));
    }
    // Check genuine relative accuracy of the small, nonzero requested MPL too.
    assert!(p.norm(&p.sub(&low[0], &high[0])) < p.tolerance(20) * p.norm(&high[0]));
    let mut compared = 0;
    for profile in fixture["profiles"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|v| v["case"] == case)
    {
        let reference = p
            .parse(
                profile["value"]["real"].as_str().unwrap(),
                profile["value"]["imaginary"].as_str().unwrap(),
            )
            .unwrap();
        let difference = p.sub(&high[0], &reference);
        assert!(p.norm(&difference) < p.tolerance(20));
        assert!(p.norm(&difference) < p.tolerance(20) * p.norm(&high[0]));
        // Source Accuracy, including finite-accuracy zero, limits comparison.
        let real_bound = p
            .parse(
                profile["real_absolute_uncertainty_estimate"]
                    .as_str()
                    .unwrap(),
                "0",
            )
            .unwrap();
        let imaginary_bound = p
            .parse(
                profile["imaginary_absolute_uncertainty_estimate"]
                    .as_str()
                    .unwrap(),
                "0",
            )
            .unwrap();
        assert!(
            p.norm(&ComplexFloat::new(difference.re.clone(), p.zero().im)) <= p.norm(&real_bound)
        );
        assert!(
            p.norm(&ComplexFloat::new(difference.im.clone(), p.zero().im))
                <= p.norm(&imaginary_bound)
        );
        compared += 1;
    }
    assert_eq!(compared, if case == "weight20" { 1 } else { 2 });
}

#[test]
fn trailing_zero_polylogarithm_uses_independent_logarithmic_boundary() {
    check_profile("trailing", 80, 50, 100, 75);
}

#[test]
fn complex_letter_polylogarithm_matches_both_original_orders() {
    check_profile("complex", 80, 50, 100, 75);
}

#[test]
fn weight_twenty_polylogarithm_matches_original_finite_accuracy_zero() {
    check_profile("weight20", 100, 100, 140, 140);
}
