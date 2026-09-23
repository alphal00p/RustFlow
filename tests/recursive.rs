use symbolica::prelude::*;
use symbolica_amflow::*;
fn sunset() -> IntegralFamily {
    IntegralFamily {
        name: "sunset".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::num(3), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::new(), &[]).unwrap(),
            Propagator::quadratic(&[1, 1], &[], Atom::new(), &[]).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}
#[test]
fn genuine_two_loop_flow_generates_vacuum_boundary() {
    let family = sunset();
    let eps = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let options = FlowOptions {
        mass_mode: MassMode::Propagators(vec![0]),
        ..Default::default()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let prepared = PreparedFlow::new(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let value = prepared
        .evaluate(&eps, &options, &provider, &context)
        .unwrap();
    let expected =
        vacuum::single_mass_sunset([1, 1, 1], &p.i(3), &Rational::from((19, 5)), p).unwrap();
    assert!(
        p.close(&value[0], &expected, 20),
        "{} != {}",
        value[0],
        expected
    );
}

#[test]
fn product_of_tadpoles_uses_recursive_soft_and_hard_families() {
    let mut family = sunset();
    family.propagators[1].constant = Atom::num(-2);
    family.propagators[2] = Propagator {
        constant: Atom::new(),
        scalar_products: vec![Atom::new(), Atom::num(1), Atom::new()],
    };
    family.physical_propagators = 2;
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let prepared = PreparedFlow::new(
        &family,
        &[Integral(vec![1, 1, 0])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let eps = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let value = prepared
        .evaluate(&eps, &options, &provider, &context)
        .unwrap();
    let d = Rational::from((19, 5));
    let expected = p.mul(
        &vacuum::tadpole(1, &p.i(3), &d, p).unwrap(),
        &vacuum::tadpole(1, &p.i(2), &d, p).unwrap(),
    );
    assert!(
        p.close(&value[0], &expected, 20),
        "{} != {expected}",
        value[0]
    );
}

#[test]
fn two_mass_sunset_recurses_and_matches_independent_beta_integral() {
    let mut family = sunset();
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-1);
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let target = Integral(vec![1, 1, 1]);
    let eps = Rational::from((3, 4));
    let p = Precision::decimal(60).unwrap();
    let prepared = PreparedFlow::new(
        &family,
        std::slice::from_ref(&target),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let value = prepared
        .evaluate(&eps, &options, &provider, &context)
        .unwrap();
    // Schwinger parameters: -Gamma(3-D)*Gamma(2-D/2)^2 /
    // [(D/2-1)*Gamma(4-D)], for two unit masses and one zero mass.
    let gamma = |n, d| {
        p.gamma_real(&p.rational(&Rational::from((n, d))).re)
            .unwrap()
    };
    let expected = p.neg(&p.div(
        &p.mul(&gamma(1, 2), &p.powi(&gamma(3, 4), 2)),
        &p.scale(&gamma(3, 2), 1, 4),
    ));
    assert!(
        p.close(&value[0], &expected, 20),
        "{} != {expected}",
        value[0]
    );
}

#[test]
fn exceptional_dimension_is_rejected_instead_of_misclassified() {
    let mut family = sunset();
    for d in &mut family.propagators {
        d.constant = Atom::num(-1);
    }
    let options = FlowOptions {
        dimension: 2,
        ..Default::default()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let prepared = PreparedFlow::new(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    assert!(
        matches!(prepared.evaluate(&Rational::from((1,2)),&options,&provider,&context),Err(Error::Unsupported(message)) if message.contains("indicial resonance"))
    );
}

#[test]
fn exact_sampled_ibp_preserves_dimensional_sectors() {
    let family = sunset();
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let epsilon = Rational::from((1, 100));
    let prepared = PreparedFlow::new_at_epsilon(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
        &epsilon,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let result = prepared
        .evaluate(&epsilon, &options, &provider, &context)
        .unwrap();
    let p = Precision::decimal(60).unwrap();
    let expected =
        vacuum::single_mass_sunset([1, 1, 1], &p.i(3), &Rational::from((199, 50)), p).unwrap();
    assert!(p.close(&result[0], &expected, 20));
    assert!(matches!(
        prepared.evaluate(&Rational::from((1, 101)), &options, &provider, &context),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn two_mass_sunset_laurent_coefficients_twenty_digits() {
    let mut family = sunset();
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-1);
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let expansion = solve_integrals(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        0,
        &options,
        &backend,
        &context,
    )
    .unwrap()
    .remove(0);
    let p = Precision::decimal(60).unwrap();
    // Independent beta-integral identity: Gamma(eps)^2/[(1-eps)(1-2eps)].
    let gamma = p
        .parse(
            "0.577215664901532860606512090082402431042159335939923598805767",
            "0",
        )
        .unwrap();
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let finite = p.add(
        &p.sub(&p.i(7), &p.scale(&gamma, 6, 1)),
        &p.add(
            &p.scale(&p.powi(&gamma, 2), 2, 1),
            &p.scale(&p.powi(&pi, 2), 1, 6),
        ),
    );
    for power in [-4, -3] {
        assert!(p.close(&expansion.coefficients[&power], &p.zero(), 20));
    }
    assert!(p.close(&expansion.coefficients[&-2], &p.i(1), 20));
    assert!(p.close(
        &expansion.coefficients[&-1],
        &p.sub(&p.i(3), &p.scale(&gamma, 2, 1)),
        20
    ));
    assert!(p.close(&expansion.coefficients[&0], &finite, 20));
    assert_eq!(expansion.verified_digits, Some(20));
}

#[test]
fn linear_numerator_boundary_keeps_the_odd_half_power_series() {
    let family = IntegralFamily {
        name: "shifted_tadpole_numerator".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-2)]],
        propagators: vec![
            Propagator {
                constant: Atom::num(-5),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
            Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::new(), Atom::num(1)],
            },
        ],
        physical_propagators: 1,
        epsilon: symbol!("numerator_eps"),
        dimension: 4,
    };
    let options = FlowOptions {
        skip_reduction: true,
        ..Default::default()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let eps = Rational::from((1, 10));
    let flow = PreparedFlow::new(
        &family,
        &[Integral(vec![1, -1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let result = flow.evaluate(&eps, &options, &provider, &context).unwrap();
    let p = Precision::decimal(60).unwrap();
    let expected = p.scale(
        &vacuum::tadpole(1, &p.i(3), &Rational::from((19, 5)), p).unwrap(),
        2,
        1,
    );
    assert!(
        p.close(&result[0], &expected, 20),
        "{} != {expected}",
        result[0]
    );
}

#[test]
fn unimodular_routing_and_denominator_permutation_preserve_flow() {
    let family = sunset();
    let mut routed = family
        .transform_loops(&[
            vec![Atom::num(1), Atom::num(1)],
            vec![Atom::new(), Atom::num(1)],
        ])
        .unwrap();
    routed.propagators = vec![
        routed.propagators[1].clone(),
        routed.propagators[2].clone(),
        routed.propagators[0].clone(),
    ];
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = FlowOptions::default();
    let eps = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let prepared = PreparedFlow::new(
        &routed,
        &[Integral(vec![1, 1, 2])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let value = prepared
        .evaluate(&eps, &options, &boundary, &context)
        .unwrap();
    let expected =
        vacuum::single_mass_sunset([2, 1, 1], &p.i(3), &Rational::from((19, 5)), p).unwrap();
    assert!(p.close(&value[0], &expected, 20));
}

#[test]
fn single_mass_three_loop_boundary_terminates_through_ft() {
    let mut propagators = [
        ([1, 0, 0], 1),
        ([0, 1, 0], 0),
        ([0, 0, 1], 0),
        ([1, 1, 1], 0),
    ]
    .into_iter()
    .map(|(loops, m)| Propagator::quadratic(&loops, &[], Atom::num(m), &[]))
    .collect::<Result<Vec<_>>>()
    .unwrap();
    for coordinate in [1, 2] {
        let mut scalar_products = vec![Atom::new(); 6];
        scalar_products[coordinate] = Atom::num(1);
        propagators.push(Propagator {
            constant: Atom::new(),
            scalar_products,
        });
    }
    let family = IntegralFamily {
        name: "single_mass_banana".into(),
        loops: vec!["l1".into(), "l2".into(), "l3".into()],
        external: vec![],
        external_gram: vec![],
        propagators,
        physical_propagators: 4,
        epsilon: symbol!("banana_eps"),
        dimension: 4,
    };
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let p = Precision::decimal(60).unwrap();
    let value = provider
        .evaluate(
            &family,
            &Integral(vec![1, 1, 1, 1, 0, 0]),
            &Rational::from((3, 4)),
            p,
        )
        .unwrap();
    let gamma = |n, d| {
        p.gamma_real(&p.rational(&Rational::from((n, d))).re)
            .unwrap()
    };
    let expected = p.div(&p.mul(&gamma(1, 2), &p.powi(&gamma(1, 4), 4)), &gamma(5, 4));
    assert!(p.close(&value, &expected, 20), "{value} != {expected}");
}

#[test]
fn equivalent_denominator_orders_share_terminal_memo() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct Terminal(AtomicUsize);
    impl recursive::TerminalProvider for Terminal {
        fn evaluate(
            &self,
            family: &IntegralFamily,
            target: &Integral,
            eps: &Rational,
            p: Precision,
        ) -> Result<Option<ComplexFloat>> {
            self.0.fetch_add(1, Ordering::Relaxed);
            let index = family
                .propagators
                .iter()
                .position(|d| !d.constant.is_zero())
                .unwrap();
            let mut powers = vec![target.0[index] as u32];
            powers.extend(
                target
                    .0
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &n)| (i != index).then_some(n as u32)),
            );
            Ok(Some(vacuum::single_mass_sunset(
                powers.try_into().unwrap(),
                &p.i(3),
                &(Rational::from(4) - eps * &Rational::from(2)),
                p,
            )?))
        }
    }
    let family = sunset();
    let mut permuted = family.clone();
    permuted.name = "renamed_family".into();
    permuted.propagators.swap(0, 1);
    let terminal = Terminal(AtomicUsize::new(0));
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = FlowOptions::default();
    let provider =
        recursive::RecursiveBoundary::new(&backend, &options, &context).with_terminal(&terminal);
    let p = Precision::decimal(60).unwrap();
    let eps = Rational::from((1, 10));
    let first = provider
        .evaluate(&family, &Integral(vec![2, 1, 1]), &eps, p)
        .unwrap();
    let second = provider
        .evaluate(&permuted, &Integral(vec![1, 2, 1]), &eps, p)
        .unwrap();
    assert_eq!(first, second);
    assert_eq!(terminal.0.load(Ordering::Relaxed), 1);
    provider
        .evaluate(
            &permuted,
            &Integral(vec![1, 2, 1]),
            &Rational::from((1, 11)),
            p,
        )
        .unwrap();
    assert_eq!(terminal.0.load(Ordering::Relaxed), 2);
}

#[test]
fn native_factorized_replay_agrees_with_sparse_reduction_and_flow() {
    let mut family = sunset();
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-1);
    let context = RunContext::default();
    let options = FlowOptions::default();
    let sparse = RustRedBackend {
        factorized: false,
        ..Default::default()
    };
    let factorized = RustRedBackend {
        factorized: true,
        max_exact_frontier: 0,
        ..Default::default()
    };
    let (deformed, _) = family
        .deform(symbol!("factorized_eta"), &MassMode::Propagators(vec![0]))
        .unwrap();
    let targets = [
        Integral(vec![2, 1, 1]),
        Integral(vec![1, 2, 1]),
        Integral(vec![1, 1, 2]),
    ];
    let a = sparse.reduce(&deformed, &targets, &context).unwrap();
    let b = factorized.reduce(&deformed, &targets, &context).unwrap();
    for target in &targets {
        let a = a.expand(target).unwrap();
        let b = b.expand(target).unwrap();
        assert_eq!(a.keys().collect::<Vec<_>>(), b.keys().collect::<Vec<_>>());
        for (i, c) in &a {
            assert!((c - &b[i]).together().cancel().is_zero());
        }
    }
    let flow = PreparedFlow::new(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        &factorized,
        &options,
        &context,
    )
    .unwrap();
    let provider = recursive::RecursiveBoundary::new(&factorized, &options, &context);
    let eps = Rational::from((3, 4));
    let value = flow.evaluate(&eps, &options, &provider, &context).unwrap();
    let p = Precision::decimal(60).unwrap();
    let gamma = |n, d| {
        p.gamma_real(&p.rational(&Rational::from((n, d))).re)
            .unwrap()
    };
    let expected = p.neg(&p.div(
        &p.mul(&gamma(1, 2), &p.powi(&gamma(3, 4), 2)),
        &p.scale(&gamma(3, 2), 1, 4),
    ));
    assert!(
        p.close(&value[0], &expected, 20),
        "{} != {expected}",
        value[0]
    );
}

#[test]
fn nonvacuum_sunrise_matches_pinned_upstream_and_precision_refinement() {
    // AMFlow 2.0 examples/differential_equation_solver: s=1/2,
    // m^2=1, epsilon=10^-4. The original input and its precision annotations
    // are preserved in fixtures/amflow-2.0/sunrise_{input,blade_sol1}.wl.
    let gram = vec![vec![Atom::num((1, 2))]];
    let propagators = [
        ([1, 0], [0], 1),
        ([0, 1], [0], 0),
        ([1, 1], [1], 0),
        ([1, 0], [-1], 0),
        ([0, 1], [-1], 0),
    ]
    .iter()
    .map(|(l, e, m)| Propagator::quadratic(l, e, Atom::num(*m), &gram).unwrap())
    .collect();
    let family = IntegralFamily {
        name: "nonvacuum_sunrise".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec!["p".into()],
        external_gram: gram,
        propagators,
        physical_propagators: 3,
        epsilon: symbol!("sunrise_eps"),
        dimension: 4,
    };
    let targets = [Integral(vec![1, 1, 1, 0, 0]), Integral(vec![2, 1, 1, 0, 0])];
    let epsilon = Rational::from((1, 10000));
    let context = RunContext::default();
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let prepared = PreparedFlow::new_at_epsilon(
        &family,
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
        &epsilon,
    )
    .unwrap();
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let first = prepared
        .evaluate(&epsilon, &options, &boundary, &context)
        .unwrap();
    let refined = FlowOptions {
        guard_digits: 60,
        series_order: 112,
        ..options.clone()
    };
    let second = prepared
        .evaluate(&epsilon, &refined, &boundary, &context)
        .unwrap();
    let p = Precision::decimal(80).unwrap();
    let references = [
        "5.0007982346841555790321435863400882030491840852200593889e7",
        "4.9999230841935509579718210700268741734828138784040806673e7",
    ];
    for ((a, b), reference) in first.iter().zip(second).zip(references) {
        assert!(p.close(a, &b, 20), "precision convergence: {a} != {b}");
        assert!(
            p.close(&b, &p.parse(reference, "0").unwrap(), 20),
            "upstream comparison: {b}"
        );
        assert!(p.norm(&ComplexFloat::new(p.real(0), b.im.clone())) < p.tolerance(20));
    }
}
