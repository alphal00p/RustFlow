use symbolica::prelude::*;
use symbolica_amflow::*;

fn bubble() -> IntegralFamily {
    IntegralFamily {
        name: "massless_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(1)]],
        propagators: vec![
            Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            Propagator {
                constant: Atom::num(1),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}

#[test]
fn bubble_automatic_flow_matches_gamma_formula() {
    let options = FlowOptions::default();
    let context = RunContext::default();
    let prepared = PreparedFlow::new(
        &bubble(),
        &[Integral(vec![1, 1])],
        &KinematicPoint::default(),
        &RustRedBackend::default(),
        &options,
        &context,
    )
    .unwrap();
    let eps = Rational::from((1, 10));
    let answer = prepared
        .evaluate(&eps, &options, &boundary::OneLoopBoundary, &context)
        .unwrap();
    let p = Precision::decimal(60).unwrap();
    let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re).unwrap();
    let g = p.div(
        &p.mul(
            &gamma(eps.clone()),
            &p.powi(&gamma(Rational::from(1) - &eps), 2),
        ),
        &gamma(Rational::from(2) - &eps * &Rational::from(2)),
    );
    let phase = p.exp(&p.mul(
        &p.complex(0, 1),
        &p.mul(
            &p.rational(&eps),
            &ComplexFloat::new(p.real(1).pi(), p.real(0)),
        ),
    ));
    let expected = p.mul(&g, &phase);
    assert!(
        p.close(&answer[0], &expected, 20),
        "{} != {expected}",
        answer[0]
    );
}

#[test]
fn laurent_fit_exact_polynomial() {
    let p = Precision::decimal(70).unwrap();
    let samples = epsilon::epsilon_samples(8, 100).unwrap();
    let values = samples
        .iter()
        .map(|eps| {
            let e = p.rational(eps);
            p.add(
                &p.scale(&p.powi(&e, -2), 7, 3),
                &p.add(
                    &p.scale(&p.powi(&e, -1), -5, 1),
                    &p.add(&p.i(11), &p.scale(&e, 13, 2)),
                ),
            )
        })
        .collect::<Vec<_>>();
    let fit = fit_epsilon(&samples, &values, -2, 1, p).unwrap();
    assert!(p.close(&fit.coefficients[&-2], &p.scale(&p.i(1), 7, 3), 50));
    assert!(p.close(&fit.coefficients[&-1], &p.i(-5), 50));
    assert!(p.close(&fit.coefficients[&0], &p.i(11), 50));
    assert!(p.close(&fit.coefficients[&1], &p.scale(&p.i(1), 13, 2), 50));
    assert_eq!(fit.verified_digits, None);
}

#[test]
fn ft_and_retained_targets_agree_with_amf_for_euclidean_massive_bubbles() {
    let mut family = bubble();
    family.external_gram[0][0] = Atom::num(-2);
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-5); // (l+p)^2 - 3
    let target = Integral(vec![1, 1]);
    let eps = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let context = RunContext::default();
    let backend = RustRedBackend::default();
    let options = FlowOptions::default();
    let ft = ft::evaluate_bubble(&family, &target, &eps, &options, &context).unwrap();
    for skip_reduction in [false, true] {
        let options = FlowOptions {
            skip_reduction,
            ..options.clone()
        };
        let prepared = PreparedFlow::new(
            &family,
            std::slice::from_ref(&target),
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )
        .unwrap();
        let amf = prepared
            .evaluate(&eps, &options, &boundary::OneLoopBoundary, &context)
            .unwrap();
        assert!(p.close(&ft, &amf[0], 20), "FT {ft}, AMF {}", amf[0]);
    }
}

#[test]
fn massless_bubble_laurent_twenty_digits() {
    let options = FlowOptions::default();
    let values = solve_integrals(
        &bubble(),
        &[Integral(vec![1, 1])],
        &KinematicPoint::default(),
        0,
        &options,
        &RustRedBackend::default(),
        &RunContext::default(),
    )
    .unwrap();
    let p = Precision::decimal(60).unwrap();
    let result = &values[0];
    assert_eq!(result.verified_digits, Some(20));
    assert!(p.close(&result.coefficients[&-2], &p.zero(), 20));
    assert!(p.close(&result.coefficients[&-1], &p.i(1), 20));
    let expected = p
        .parse(
            "1.4227843350984671393934879099175975689578406640600764",
            "3.1415926535897932384626433832795028841971693993751058",
        )
        .unwrap();
    assert!(p.close(&result.coefficients[&0], &expected, 20));
}

#[test]
fn exact_complex_tadpole_and_opposite_prescription() {
    let p = Precision::decimal(60).unwrap();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let family = IntegralFamily {
        name: "complex_tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: Atom::num(symbolica::domains::float::Complex::new(
                Rational::from(-2),
                Rational::from(1),
            )),
            scalar_products: vec![Atom::num(1)],
        }],
        physical_propagators: 1,
        epsilon: symbol!("complex_eps"),
        dimension: 4,
    };
    let options = FlowOptions::default();
    let eps = Rational::from((1, 10));
    let prepared = PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let value = prepared
        .evaluate(&eps, &options, &boundary::OneLoopBoundary, &context)
        .unwrap();
    let expected = vacuum::tadpole(1, &p.complex(2, -1), &Rational::from((19, 5)), p).unwrap();
    assert!(
        p.close(&value[0], &expected, 20),
        "{} != {expected}",
        value[0]
    );
    let family = bubble();
    let prepared = PreparedFlow::new(
        &family,
        &[Integral(vec![1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let plus = prepared
        .evaluate(&eps, &options, &boundary::OneLoopBoundary, &context)
        .unwrap();
    let minus = prepared
        .evaluate(
            &eps,
            &FlowOptions {
                prescription: Prescription::MinusI0,
                ..options
            },
            &boundary::OneLoopBoundary,
            &context,
        )
        .unwrap();
    assert!(p.close(
        &plus[0],
        &ComplexFloat::new(minus[0].re.clone(), -minus[0].im.clone()),
        20
    ));
}

#[test]
fn physical_massless_triangle_matches_independent_gamma_formula() {
    let gram = vec![
        vec![Atom::new(), Atom::num(15)],
        vec![Atom::num(15), Atom::new()],
    ];
    let family = IntegralFamily {
        name: "physical_triangle".into(),
        loops: vec!["l".into()],
        external: vec!["p1".into(), "p2".into()],
        external_gram: gram.clone(),
        propagators: [[0, 0], [1, 0], [1, 1]]
            .iter()
            .map(|external| Propagator::quadratic(&[1], external, Atom::new(), &gram).unwrap())
            .collect(),
        physical_propagators: 3,
        epsilon: symbol!("triangle_eps"),
        dimension: 4,
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = FlowOptions::default();
    let flow = PreparedFlow::new(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let p = Precision::decimal(60).unwrap();
    let eps = Rational::from((1, 10));
    let value = flow.evaluate(&eps, &options, &boundary, &context).unwrap();
    let gamma = |n, d| {
        p.gamma_real(&p.rational(&Rational::from((n, d))).re)
            .unwrap()
    };
    let prefactor = p.div(
        &p.mul(&gamma(11, 10), &p.powi(&gamma(-1, 10), 2)),
        &gamma(4, 5),
    );
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let expected = p.mul(
        &p.mul(
            &prefactor,
            &p.pow(&p.i(30), &p.rational(&Rational::from((-11, 10)))),
        ),
        &p.exp(&p.mul(&p.complex(0, 1), &p.scale(&pi, 1, 10))),
    );
    assert!(
        p.close(&value[0], &expected, 20),
        "{} != {expected}",
        value[0]
    );
}

#[test]
fn alternative_mass_placements_agree_for_euclidean_bubble() {
    let mut family = bubble();
    family.external_gram[0][0] = Atom::num(-2);
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-5);
    let target = Integral(vec![1, 1]);
    let epsilon = Rational::from((1, 10));
    let precision = Precision::decimal(60).unwrap();
    let context = RunContext::default();
    let backend = RustRedBackend::default();
    let expected = ft::evaluate_bubble(
        &family,
        &target,
        &epsilon,
        &FlowOptions::default(),
        &context,
    )
    .unwrap();
    for mass_mode in [
        MassMode::Mass,
        MassMode::Propagator,
        MassMode::Branch,
        MassMode::Loop,
    ] {
        let options = FlowOptions {
            mass_mode,
            ..Default::default()
        };
        let prepared = PreparedFlow::new(
            &family,
            std::slice::from_ref(&target),
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )
        .unwrap();
        let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
        let actual = prepared
            .evaluate(&epsilon, &options, &boundary, &context)
            .unwrap();
        assert!(
            precision.close(&actual[0], &expected, 20),
            "{:?}: {} != {expected}",
            options.mass_mode,
            actual[0]
        );
    }
}

#[test]
fn incomparable_top_sectors_keep_their_independent_normalizations() {
    let mut family = bubble();
    family.external_gram[0][0] = Atom::num(-2);
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-5);
    let targets = [Integral(vec![1, 0]), Integral(vec![0, 1])];
    let options = FlowOptions::default();
    let context = RunContext::default();
    let backend = RustRedBackend::default();
    let epsilon = Rational::from((1, 10));
    let precision = Precision::decimal(60).unwrap();
    let prepared = PreparedFlow::new(
        &family,
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    assert_eq!(prepared.reduced.basis.len(), 2);
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let actual = prepared
        .evaluate(&epsilon, &options, &boundary, &context)
        .unwrap();
    for (value, mass_squared) in actual.iter().zip([1, 3]) {
        let expected = vacuum::tadpole(
            1,
            &precision.i(mass_squared),
            &Rational::from((19, 5)),
            precision,
        )
        .unwrap();
        assert!(precision.close(value, &expected, 20));
    }
}
