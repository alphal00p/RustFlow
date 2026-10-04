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
    let gamma = ComplexFloat::new(
        Float::from_raw(rug::Float::with_val(p.bits, rug::float::Constant::Euler)),
        p.real(0),
    );
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
