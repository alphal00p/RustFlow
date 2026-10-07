use symbolica::prelude::*;
use symbolica_amflow::linear::PreparedLinearFlow;
use symbolica_amflow::*;

// q², 2q.v-2, 2q.w-2, with v²=w²=0 and v.w=1.
fn family() -> IntegralFamily {
    IntegralFamily {
        name: "independent_eikonal_lines".into(),
        loops: vec!["k".into()],
        external: vec!["v".into(), "w".into()],
        external_gram: vec![
            vec![Atom::num(0), Atom::num(1)],
            vec![Atom::num(1), Atom::num(0)],
        ],
        propagators: vec![
            Propagator {
                constant: Atom::num(0),
                scalar_products: vec![Atom::num(1), Atom::num(0), Atom::num(0)],
            },
            Propagator {
                constant: Atom::num(-2),
                scalar_products: vec![Atom::num(0), Atom::num(2), Atom::num(0)],
            },
            Propagator {
                constant: Atom::num(-2),
                scalar_products: vec![Atom::num(0), Atom::num(0), Atom::num(2)],
            },
        ],
        physical_propagators: 3,
        epsilon: symbol!("linear_two_directions::eps"),
        dimension: 4,
    }
}

#[test]
fn two_independent_eikonal_directions() {
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let prepared = PreparedLinearFlow::new(
        &family(),
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let epsilon = Rational::from((1, 3));
    let actual = prepared
        .evaluate(&epsilon, &options, &backend, &context)
        .unwrap();
    let p = Precision::decimal(80).unwrap();
    let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re).unwrap();
    let expected = p.neg(&p.mul(
        &p.pow(&p.i(2), &p.rational(&(-Rational::from(1) - &epsilon))),
        &p.mul(
            &p.mul(&gamma(epsilon.clone()), &gamma(epsilon.clone())),
            &gamma(Rational::from(1) - &epsilon),
        ),
    ));
    assert!(
        p.close(&actual.values[0], &expected, 20),
        "{} != {expected}",
        actual.values[0]
    );
}

#[test]
fn ft_two_independent_eikonal_directions() {
    let options = FlowOptions {
        recursion: RecursionMode::Ft,
        ..FlowOptions::default()
    };
    let refined = FlowOptions {
        guard_digits: 60,
        series_order: 112,
        ..options.clone()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let p = Precision::decimal(80).unwrap();
    let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re).unwrap();
    let oracle: serde_json::Value = serde_json::from_str(include_str!(
        "../reports/validation/2026-10-06-linear-ft/oracle-metadata.json"
    ))
    .unwrap();
    assert_eq!(oracle["requested_digits"], 24);
    let stored = oracle["coefficients"][0].as_str().unwrap();
    let reference = p.parse(stored.split('`').next().unwrap(), "0").unwrap();
    for epsilon in [Rational::from((1, 3)), Rational::from((2, 5))] {
        let expected = p.neg(&p.mul(
            &p.pow(&p.i(2), &p.rational(&(-Rational::from(1) - &epsilon))),
            &p.mul(
                &p.powi(&gamma(epsilon.clone()), 2),
                &gamma(Rational::from(1) - &epsilon),
            ),
        ));
        if epsilon == (1, 3) {
            assert!(p.close(&expected, &reference, 20));
        }
        let mut original = None;
        for settings in [&options, &refined] {
            // Fresh evaluator: refinement cannot reuse a numerical memo entry.
            let evaluator = ft::FtEvaluator::new(&backend, &context);
            let actual = evaluator
                .evaluate(&family(), &Integral(vec![1, 1, 1]), &epsilon, settings)
                .unwrap();
            assert!(p.close(&actual, &expected, 20), "{actual} != {expected}");
            if let Some(first) = &original {
                assert!(p.close(&actual, first, 20));
            }
            original = Some(actual);
            // Raising the first linear power differentiates its residual energy.
            let raised = evaluator
                .evaluate(&family(), &Integral(vec![1, 2, 1]), &epsilon, settings)
                .unwrap();
            let raised_expected = p.neg(&p.scale(&p.mul(&p.rational(&epsilon), &expected), 1, 2));
            assert!(
                p.close(&raised, &raised_expected, 20),
                "{raised} != {raised_expected}"
            );
        }
    }
}

#[test]
fn ft_linear_recursion_preserves_the_euclidean_domain_guard() {
    let mut family = family();
    family.external_gram[0][1] = Atom::num(-1);
    family.external_gram[1][0] = Atom::num(-1);
    let options = FlowOptions {
        recursion: RecursionMode::Ft,
        ..Default::default()
    };
    // Here M²=2t+2u-2tu has a positive-parameter zero. The high-level
    // dispatcher must retain FT's contour guard, rather than silently using
    // the AMF algorithm selected for other linear evaluations.
    assert!(
        matches!(solve_integrals(&family, &[Integral(vec![1,1,1])], &KinematicPoint::default(),0,&options,&RustRedBackend::default(),&RunContext::default()),Err(Error::Unsupported(message)) if message.contains("Euclidean F"))
    );
}

#[test]
fn ft_mixed_loop_linear_denominator() {
    // l1², l1.l2, l2², l1.v, l1.w, l2.v, l2.w.
    let mut propagators = vec![
        Propagator {
            constant: Atom::num(0),
            scalar_products: [1, 0, 0, 0, 0, 0, 0].map(Atom::num).to_vec(),
        },
        Propagator {
            constant: Atom::num(0),
            scalar_products: [0, 0, 1, 0, 0, 0, 0].map(Atom::num).to_vec(),
        },
        Propagator {
            constant: Atom::num(-2),
            scalar_products: [0, 0, 0, 2, 0, 0, 2].map(Atom::num).to_vec(),
        },
    ];
    for coordinate in [1, 4, 5, 6] {
        let mut scalar_products = vec![Atom::num(0); 7];
        scalar_products[coordinate] = Atom::num(1);
        propagators.push(Propagator {
            constant: Atom::num(0),
            scalar_products,
        });
    }
    let family = IntegralFamily {
        name: "mixed_loop_linear".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec!["v".into(), "w".into()],
        external_gram: vec![
            vec![Atom::num(1), Atom::num(0)],
            vec![Atom::num(0), Atom::num(1)],
        ],
        propagators,
        physical_propagators: 3,
        epsilon: symbol!("linear_two_directions::eps"),
        dimension: 4,
    };
    // This linear coefficient matrix has rank two. AMF's squared-branch
    // deformation must remain unavailable, while FT does not need it.
    assert!(matches!(
        linear::LinearDeformation::new(&family, symbol!("linear_two_directions::x")),
        Err(Error::Unsupported(_))
    ));
    let options = FlowOptions {
        recursion: RecursionMode::Ft,
        ..FlowOptions::default()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let epsilon = Rational::from((7, 8));
    let evaluator = ft::FtEvaluator::new(&backend, &context);
    let actual = evaluator
        .evaluate(
            &family,
            &Integral(vec![1, 1, 1, 0, 0, 0, 0]),
            &epsilon,
            &options,
        )
        .unwrap();
    let p = Precision::decimal(80).unwrap();
    let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re).unwrap();
    let expected = p.neg(&p.mul(
        &p.pow(
            &p.i(2),
            &p.rational(&(Rational::from(3) - &epsilon * &Rational::from(4))),
        ),
        &p.mul(
            &p.powi(&gamma(Rational::from(1) - &epsilon), 2),
            &gamma(&epsilon * &Rational::from(4) - Rational::from(3)),
        ),
    ));
    assert!(p.close(&actual, &expected, 20), "{actual} != {expected}");
    let refined = FlowOptions {
        guard_digits: 60,
        series_order: 112,
        ..options.clone()
    };
    let repeated = ft::FtEvaluator::new(&backend, &context)
        .evaluate(
            &family,
            &Integral(vec![1, 1, 1, 0, 0, 0, 0]),
            &epsilon,
            &refined,
        )
        .unwrap();
    assert!(p.close(&repeated, &expected, 20));
    assert!(p.close(&repeated, &actual, 20));
    // A mixed routing with |det T|=2 tests the D-dimensional integration
    // measure, not a family-name or denominator-pattern identity.
    let routed = family
        .transform_loops(&[
            vec![Atom::num(1), Atom::num(1)],
            vec![Atom::num(1), Atom::num(-1)],
        ])
        .unwrap();
    let routed_value = ft::FtEvaluator::new(&backend, &context)
        .evaluate(
            &routed,
            &Integral(vec![1, 1, 1, 0, 0, 0, 0]),
            &epsilon,
            &refined,
        )
        .unwrap();
    let dimension = Rational::from(4) - &epsilon * &Rational::from(2);
    let jacobian = p.pow(&p.i(2), &p.rational(&dimension));
    assert!(p.close(&p.mul(&routed_value, &jacobian), &expected, 20));
}

#[test]
fn public_ft_linear_samples_are_fitted_after_exact_projection() {
    let family = family();
    let options = FlowOptions {
        recursion: RecursionMode::Ft,
        workers: 2,
        ..FlowOptions::default()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let integral = Integral(vec![1, 1, 1]);
    let result = solve_integrals(
        &family,
        std::slice::from_ref(&integral),
        &KinematicPoint::default(),
        0,
        &options,
        &backend,
        &context,
    )
    .unwrap();
    let p = Precision::decimal(80).unwrap();
    let gamma = ComplexFloat::new(p.real(0).euler(), p.real(0));
    let pi = ComplexFloat::new(p.real(0).pi(), p.real(0));
    let logarithm = p.add(&gamma, &p.log(&p.i(2)));
    let finite = p.neg(&p.add(
        &p.scale(&p.powi(&logarithm, 2), 1, 4),
        &p.scale(&p.powi(&pi, 2), 1, 8),
    ));
    assert_eq!(result[0].verified_digits, Some(20));
    for (power, expected) in [
        (-2, p.scale(&p.i(-1), 1, 2)),
        (-1, p.scale(&logarithm, 1, 2)),
        (0, finite),
    ] {
        assert!(p.close(&result[0].coefficients[&power], &expected, 20));
    }
    let projected = solve_integral_projections(
        &[(
            family.clone(),
            vec![std::collections::BTreeMap::from([(
                integral,
                Atom::var(family.epsilon).pow(2),
            )])],
        )],
        &KinematicPoint::default(),
        0,
        &options,
        &backend,
        &context,
    )
    .unwrap();
    assert_eq!(projected[0].verified_digits, Some(20));
    assert!(p.close(&projected[0].coefficients[&0], &p.scale(&p.i(-1), 1, 2), 20));
}
