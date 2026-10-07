use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::*;

fn tadpole(mass_squared: i64) -> IntegralFamily {
    IntegralFamily {
        name: format!("combination_tadpole_{mass_squared}"),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: Atom::num(-mass_squared),
            scalar_products: vec![Atom::num(1)],
        }],
        physical_propagators: 1,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}

#[test]
fn epsilon_pole_weights_and_cross_family_cancellation_precede_fitting() {
    // [I(m²=1)-I(m²=4)/4]/eps = Gamma(eps)/(eps*(1-eps))*(1-4^-eps).
    // Higher poles cancel between distinct families; weights must be evaluated
    // before fitting, with their extra epsilon pole included in the fit range.
    let result = solve_integral_combinations(
        &[
            (
                tadpole(1),
                BTreeMap::from([(Integral(vec![1]), parse!("1/eps"))]),
            ),
            (
                tadpole(4),
                BTreeMap::from([(Integral(vec![1]), parse!("-1/(4*eps)"))]),
            ),
        ],
        &KinematicPoint::default(),
        0,
        &FlowOptions {
            workers: 2,
            ..Default::default()
        },
        &RustRedBackend::default(),
        &RunContext::default(),
    )
    .unwrap();
    let p = Precision::decimal(80).unwrap();
    let gamma = ComplexFloat::new(p.real(0).euler(), p.real(0));
    let logarithm = p.log(&p.i(4));
    let finite = p.sub(
        &p.mul(&logarithm, &p.sub(&p.i(1), &gamma)),
        &p.scale(&p.mul(&logarithm, &logarithm), 1, 2),
    );
    assert_eq!(result.verified_digits, Some(20));
    for pole in [-3, -2] {
        assert!(p.close(&result.coefficients[&pole], &p.zero(), 20));
    }
    assert!(p.close(&result.coefficients[&-1], &logarithm, 20));
    assert!(p.close(&result.coefficients[&0], &finite, 20));
}

#[test]
fn zero_combinations_keep_the_familys_laurent_range_and_inexact_weights_are_rejected() {
    let evaluate = |weight| {
        solve_integral_combinations(
            &[(tadpole(1), BTreeMap::from([(Integral(vec![1]), weight)]))],
            &KinematicPoint::default(),
            -1,
            &FlowOptions::default(),
            &RustRedBackend::default(),
            &RunContext::default(),
        )
    };
    let zero = evaluate(Atom::new()).unwrap();
    assert_eq!(
        zero.coefficients.keys().copied().collect::<Vec<_>>(),
        [-2, -1]
    );
    assert!(zero.coefficients.values().all(|value| value.is_zero()));
    assert!(matches!(
        evaluate(Atom::num(Float::with_val(80, 1))),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn shared_integral_projections_preserve_cancellations_and_epsilon_prefactors() {
    let rows = vec![
        BTreeMap::from([
            (Integral(vec![1]), parse!("1-eps")),
            (Integral(vec![2]), Atom::num(-1)),
        ]),
        BTreeMap::from([(Integral(vec![2]), parse!("eps"))]),
        BTreeMap::from([(Integral(vec![1]), parse!("eps"))]),
    ];
    let values = solve_integral_projections(
        &[(tadpole(1), rows)],
        &KinematicPoint::default(),
        1,
        &FlowOptions::default(),
        &RustRedBackend::default(),
        &RunContext::default(),
    )
    .unwrap();
    let p = Precision::decimal(80).unwrap();
    let gamma = ComplexFloat::new(p.real(0).euler(), p.real(0));
    assert_eq!(values.len(), 3);
    for result in &values {
        assert_eq!(result.verified_digits, Some(20));
    }
    // The first row is an exact IBP cancellation; the other two have finite
    // terms supplied by their epsilon prefactors multiplying divergent masters.
    assert!(
        values[0]
            .coefficients
            .values()
            .all(|c| p.close(c, &p.zero(), 20))
    );
    for result in &values[1..] {
        assert!(p.close(&result.coefficients[&0], &p.i(1), 20));
        assert!(p.close(&result.coefficients[&-1], &p.zero(), 20));
    }
    assert!(p.close(&values[1].coefficients[&1], &p.neg(&gamma), 20));
    assert!(p.close(&values[2].coefficients[&1], &p.sub(&p.i(1), &gamma), 20));
}

#[test]
fn projection_rows_sum_across_families_and_keep_zero_outputs() {
    let groups = [
        (
            tadpole(1),
            vec![
                BTreeMap::from([(Integral(vec![1]), parse!("eps"))]),
                BTreeMap::new(),
            ],
        ),
        (
            tadpole(4),
            vec![
                BTreeMap::from([(Integral(vec![1]), parse!("-eps/4"))]),
                BTreeMap::new(),
            ],
        ),
    ];
    let values = solve_integral_projections(
        &groups,
        &KinematicPoint::default(),
        1,
        &FlowOptions::default(),
        &RustRedBackend::default(),
        &RunContext::default(),
    )
    .unwrap();
    let p = Precision::decimal(80).unwrap();
    assert!(p.close(&values[0].coefficients[&0], &p.zero(), 20));
    assert!(p.close(&values[0].coefficients[&1], &p.log(&p.i(4)), 20));
    assert!(values[1].coefficients.values().all(|v| v.is_zero()));
    let bad = [(tadpole(1), vec![BTreeMap::new()]), (tadpole(4), vec![])];
    assert!(matches!(
        solve_integral_projections(
            &bad,
            &KinematicPoint::default(),
            0,
            &FlowOptions::default(),
            &RustRedBackend::default(),
            &RunContext::default(),
        ),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn algebraic_projection_and_gamma_normalization_precede_fitting() {
    use symbolica_amflow::algebraic::SquareRoot;
    use symbolica_amflow::transport_cache::{RootGerm, RootSheet};
    let root = symbol!("projection_r");
    let mut factors = ProjectionFactors::with_roots(
        vec![SquareRoot {
            symbol: root,
            radicand: Atom::num(-2),
        }],
        RootGerm {
            sheets: BTreeMap::from([(root, RootSheet::Opposite)]),
        },
    );
    factors.normalization = vec![
        SampleNormalization::Gamma {
            offset: Rational::from(1),
            slope: Rational::from(1),
            power: -1,
        },
        SampleNormalization::EulerGammaExponential(Rational::from(2)),
        SampleNormalization::Power {
            base: Atom::num(-3),
            slope: Rational::from(1),
            prescription: Prescription::MinusI0,
        },
    ];
    // eps*I_2 = Gamma(1+eps). After cancelling that gamma, the result is
    // -i sqrt(2) exp(eps*(2 EulerGamma + log(3) - i pi)).
    let result = solve_integral_projections_normalized(
        &[(
            tadpole(1),
            vec![BTreeMap::from([(
                Integral(vec![2]),
                Atom::var(root) * parse!("eps"),
            )])],
        )],
        &KinematicPoint::default(),
        1,
        &FlowOptions::default(),
        &RustRedBackend::default(),
        &RunContext::default(),
        &factors,
    )
    .unwrap()
    .remove(0);
    let p = Precision::decimal(100).unwrap();
    let r = ComplexFloat::new(p.real(0), -p.real(2).sqrt());
    let gamma = ComplexFloat::new(p.real(0).euler(), p.real(0));
    let log = ComplexFloat::new(p.log(&p.i(3)).re, -p.real(0).pi());
    assert_eq!(result.verified_digits, Some(20));
    assert!(p.close(&result.coefficients[&0], &r, 20));
    assert!(p.close(
        &result.coefficients[&1],
        &p.mul(&r, &p.add(&p.scale(&gamma, 2, 1), &log)),
        20
    ));
    factors.sheets.clear();
    assert!(matches!(
        solve_integral_projections_normalized(
            &[(tadpole(1), vec![BTreeMap::new()])],
            &KinematicPoint::default(),
            0,
            &FlowOptions::default(),
            &RustRedBackend::default(),
            &RunContext::default(),
            &factors,
        ),
        Err(Error::InvalidInput(_))
    ));
}
