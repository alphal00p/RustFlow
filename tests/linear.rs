use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::linear::{LinearDeformation, PreparedLinearFlow};
use symbolica_amflow::*;

fn equal(a: &Atom, b: &Atom) {
    assert!((a - b).together().cancel().is_zero(), "{a} != {b}");
}

fn hqet() -> IntegralFamily {
    IntegralFamily {
        name: "one_loop_hqet".into(),
        loops: vec!["l".into()],
        external: vec!["v".into()],
        external_gram: vec![vec![Atom::num(1)]],
        propagators: vec![
            Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            Propagator {
                constant: Atom::num(-2),
                scalar_products: vec![Atom::new(), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("linear_test::eps"),
        dimension: 4,
    }
}

// Complete a chosen loop-linear denominator by unchanged coordinate ISPs.
fn linear_basis(loops: usize, external: usize, matrix: Vec<Atom>) -> IntegralFamily {
    let offset = loops * (loops + 1) / 2;
    let count = offset + loops * external;
    assert_eq!(matrix.len(), loops * external);
    let omitted = offset + matrix.iter().position(|a| !a.is_zero()).unwrap();
    let mut products = vec![Atom::new(); offset];
    products.extend(matrix);
    let mut propagators = vec![Propagator {
        constant: Atom::num(3),
        scalar_products: products,
    }];
    for coordinate in 0..count {
        if coordinate == omitted {
            continue;
        }
        let mut scalar_products = vec![Atom::new(); count];
        scalar_products[coordinate] = Atom::num(1);
        propagators.push(Propagator {
            constant: Atom::new(),
            scalar_products,
        });
    }
    IntegralFamily {
        name: "linear_branch_basis".into(),
        loops: (0..loops).map(|i| format!("l{i}")).collect(),
        external: (0..external).map(|i| format!("p{i}")).collect(),
        external_gram: (0..external)
            .map(|i| {
                (0..external)
                    .map(|j| Atom::num(i64::from(i == j)))
                    .collect()
            })
            .collect(),
        propagators,
        physical_propagators: 1,
        epsilon: symbol!("linear_test::eps"),
        dimension: 4,
    }
}

#[test]
fn rank_one_branch_deformation_preserves_slots_numerators_and_exact_endpoint() {
    let family = linear_basis(
        2,
        2,
        vec![Atom::num(2), Atom::num(6), Atom::num(-4), Atom::num(-12)],
    );
    let x = symbol!("linear_test::x");
    let deformation = LinearDeformation::new(&family, x).unwrap().unwrap();
    assert_eq!(deformation.lines.len(), 1);
    assert_eq!(deformation.lines[0].propagator, 0);
    assert_eq!(
        deformation.lines[0].loop_coefficients,
        vec![Atom::num(1), Atom::num(-2)]
    );
    let quadratic = &deformation.family.propagators[0].scalar_products;
    equal(&quadratic[0], &Atom::var(x));
    equal(&quadratic[1], &(Atom::num(-4) * Atom::var(x)));
    equal(&quadratic[2], &(Atom::num(4) * Atom::var(x)));
    assert_eq!(&quadratic[3..], &family.propagators[0].scalar_products[3..]);
    let original = deformation.family.at(&KinematicPoint(BTreeMap::from([(
        Atom::var(x),
        Atom::new(),
    )])));
    for (i, (a, b)) in original
        .propagators
        .iter()
        .zip(&family.propagators)
        .enumerate()
    {
        equal(&a.constant, &b.constant);
        for (actual, expected) in a.scalar_products.iter().zip(&b.scalar_products) {
            equal(actual, expected);
        }
        if i > 0 {
            assert_eq!(
                deformation.family.propagators[i].scalar_products,
                b.scalar_products
            );
        }
    }
    // Rescaling a linear denominator changes its coefficients, while the
    // normalized auxiliary branch and its x=0 denominator remain exact.
    let mut rescaled = family.clone();
    for coefficient in &mut rescaled.propagators[0].scalar_products {
        *coefficient = &*coefficient * Atom::num(7);
    }
    rescaled.propagators[0].constant *= Atom::num(7);
    let rescaled = LinearDeformation::new(&rescaled, x).unwrap().unwrap();
    assert_eq!(
        rescaled.lines[0].loop_coefficients,
        deformation.lines[0].loop_coefficients
    );
}

#[test]
fn unsupported_branch_geometry_and_deformation_symbol_collisions_are_rejected() {
    let x = symbol!("linear_test::x");
    let rank_two = linear_basis(
        2,
        2,
        vec![Atom::num(1), Atom::new(), Atom::new(), Atom::num(1)],
    );
    assert!(matches!(
        LinearDeformation::new(&rank_two, x),
        Err(Error::Unsupported(_))
    ));
    let imaginary = Atom::num(symbolica::domains::float::Complex::new(
        Rational::from(0),
        Rational::from(1),
    ));
    let complex_branch = linear_basis(2, 1, vec![Atom::num(1), imaginary]);
    assert!(matches!(
        LinearDeformation::new(&complex_branch, x),
        Err(Error::Unsupported(_))
    ));
    let mut family = hqet();
    assert!(matches!(
        LinearDeformation::new(&family, family.epsilon),
        Err(Error::InvalidInput(_))
    ));
    family.propagators[1].constant = Atom::var(x);
    assert!(matches!(
        LinearDeformation::new(&family, x),
        Err(Error::InvalidInput(_))
    ));
    family.propagators[1].scalar_products[0] = Atom::num(1);
    family.propagators[1].constant = Atom::num(-2);
    assert!(LinearDeformation::new(&family, x).unwrap().is_none());
}

#[test]
fn original_denominator_domain_is_retained_after_branch_normalization() {
    let u_symbol = symbol!("linear_test::u");
    let u = Atom::var(u_symbol);
    let family = linear_basis(2, 1, vec![Atom::num(1) / &u, Atom::num(2) / &u]);
    let deformation = LinearDeformation::new(&family, symbol!("linear_test::x"))
        .unwrap()
        .unwrap();
    assert_eq!(
        deformation.lines[0].loop_coefficients,
        vec![Atom::num(1), Atom::num(2)]
    );
    let zero = KinematicPoint(BTreeMap::from([(u, Atom::new())]));
    assert!(deformation.nonzero_conditions.iter().any(|a| {
        !a.derivative(u_symbol).is_zero() && zero.apply(a).together().cancel().is_zero()
    }));
}

#[test]
fn deformation_derivative_is_the_exact_eikonal_propagator_identity() {
    let x = symbol!("linear_test::x");
    let deformation = LinearDeformation::new(&hqet(), x).unwrap().unwrap();
    // d/dx [D1^-a (D2+x D1)^-b] = -b D1^-(a-1) (D2+x D1)^-(b+1).
    for (a, b) in [(1, 1), (3, 2), (-2, 3)] {
        let derivative =
            reduction::parameter_derivative(&deformation.family, &Integral(vec![a, b]), x).unwrap();
        assert_eq!(derivative.len(), 1);
        assert_eq!(derivative[0].0, Integral(vec![a - 1, b + 1]));
        equal(&derivative[0].1, &Atom::num(-i64::from(b)));
    }
}

#[test]
fn automatic_hqet_limit_matches_gamma_formula_and_precision_refinement() {
    let family = hqet();
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let flow = PreparedLinearFlow::new(
        &family,
        &[Integral(vec![1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let refined = FlowOptions {
        guard_digits: 60,
        series_order: 112,
        ..options.clone()
    };
    let p = Precision::decimal(80).unwrap();
    // For D=4-2eps and v^2=1, positive Feynman parameter t gives
    // (l^2+t(2v.l-2)) = (l+tv)^2-t(t+2).
    // Its tadpole integral and the beta integral yield
    // 2^(1-2eps) Gamma(1-eps) Gamma(2eps-1), in d^Dl/(i pi^(D/2)).
    // eps=2/3 lies in the convergence strip; eps=1/10 checks continuation.
    for epsilon in [Rational::from((2, 3)), Rational::from((1, 10))] {
        let first = flow
            .evaluate(&epsilon, &options, &backend, &context)
            .unwrap();
        let second = flow
            .evaluate(&epsilon, &refined, &backend, &context)
            .unwrap();
        assert!(first.reference_point > 0);
        assert_eq!(first.reference_point, second.reference_point);
        assert!(first.reference_basis_size > 0);
        let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re).unwrap();
        let expected = p.mul(
            &p.pow(
                &p.i(2),
                &p.rational(&(Rational::from(1) - &epsilon * &Rational::from(2))),
            ),
            &p.mul(
                &gamma(Rational::from(1) - &epsilon),
                &gamma(&epsilon * &Rational::from(2) - Rational::from(1)),
            ),
        );
        assert!(
            p.close(&first.values[0], &expected, 20),
            "epsilon={epsilon}: {} != {expected}",
            first.values[0]
        );
        assert!(
            p.close(&second.values[0], &expected, 20),
            "epsilon={epsilon}: {} != {expected}",
            second.values[0]
        );
        assert!(p.close(&first.values[0], &second.values[0], 20));
    }
}
