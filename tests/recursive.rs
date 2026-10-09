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
        max_sector_batch: 1,
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

#[test]
fn physical_paper_subsector_closes_with_deeper_ibps_and_stable_phase() {
    // Depth two can leave a redundant direction with a nonphysical indicial
    // root in this sector. A deeper search removes it before region matching.
    let (family, _) = benchmarks::paper_two_loop().unwrap();
    let target = Integral(vec![1, 0, 1, 0, 1, 0, 1, 0, 0]);
    let backend = RustRedBackend {
        max_depth: 3,
        max_targets: 32768,
        max_sector_batch: 32,
        ..Default::default()
    };
    let epsilon = Rational::from((1, 2700));
    let options = FlowOptions::default();
    let context = RunContext::default();
    let prepared = PreparedFlow::new_at_epsilon(
        &family,
        &[target],
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
    assert!(p.close(&first[0], &second[0], 20));
    assert!(second[0].im > p.real(0));
    // Regression value from the Rust higher-precision calculation, not an
    // upstream reference and never used to supply a boundary condition.
    let expected = p
        .parse(
            "3641004.8899918837093934440970600173922862936154440",
            "8471.9396295136564478613884972476665495673226018302",
        )
        .unwrap();
    assert!(p.close(&second[0], &expected, 20));
}

#[test]
fn recursive_sampled_flows_keep_epsilon_in_memo_keys_and_handle_exceptional_samples() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct Recording {
        sampled: AtomicUsize,
        symbolic: AtomicUsize,
    }
    impl ReductionBackend for Recording {
        fn identity(&self) -> String {
            "recursive-sampling-test".into()
        }
        fn reduce(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            if family.loops.len() > 1 {
                self.symbolic.fetch_add(1, Ordering::Relaxed);
            }
            RustRedBackend::default().reduce(family, targets, context)
        }
        fn reduce_at_epsilon(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            epsilon: &Rational,
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            if family.loops.len() > 1 {
                self.sampled.fetch_add(1, Ordering::Relaxed);
            }
            RustRedBackend::default().reduce_at_epsilon(family, targets, epsilon, context)
        }
    }
    let backend = Recording {
        sampled: AtomicUsize::new(0),
        symbolic: AtomicUsize::new(0),
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let evaluator = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let mut family = sunset();
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-1);
    let target = Integral(vec![1, 1, 1]);
    let p = Precision::decimal(60).unwrap();
    let check = |epsilon: Rational| {
        let value = evaluator.evaluate(&family, &target, &epsilon, p).unwrap();
        let e = p.rational(&epsilon);
        let gamma = p.gamma_real(&e.re).unwrap();
        let expected = p.div(
            &p.powi(&gamma, 2),
            &p.mul(&p.sub(&p.i(1), &e), &p.sub(&p.i(1), &p.scale(&e, 2, 1))),
        );
        assert!(p.close(&value, &expected, 20));
    };
    check(Rational::from((1, 100)));
    let first_calls = backend.sampled.load(Ordering::Relaxed);
    assert!(first_calls > 0);
    check(Rational::from((1, 101)));
    assert!(backend.sampled.load(Ordering::Relaxed) > first_calls);
    let second_calls = backend.sampled.load(Ordering::Relaxed);
    check(Rational::from((1, 100)));
    assert_eq!(backend.sampled.load(Ordering::Relaxed), second_calls);
    assert_eq!(backend.symbolic.load(Ordering::Relaxed), 0);
    check(Rational::from((3, 4)));
    assert!(backend.symbolic.load(Ordering::Relaxed) > 0);
    assert_eq!(backend.sampled.load(Ordering::Relaxed), second_calls);
}

#[test]
fn boundary_batches_share_preparations_and_preserve_mass_choices() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    struct CountedBackend {
        backend: RustRedBackend,
        calls: AtomicUsize,
    }
    impl reduction::ReductionBackend for CountedBackend {
        fn identity(&self) -> String {
            self.backend.identity()
        }
        fn reduce(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            self.backend.reduce(family, targets, context)
        }
        fn reduce_at_epsilon(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            epsilon: &Rational,
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            self.backend
                .reduce_at_epsilon(family, targets, epsilon, context)
        }
    }
    let mut family = sunset();
    family.propagators[1].constant = Atom::num(-2);
    family.propagators[2].constant = Atom::num(-1);
    // Three distinct products of massive tadpoles. The second target must
    // shift line 1; the other two shift line 0 and can share one preparation.
    let targets = vec![
        Integral(vec![1, 1, 0]),
        Integral(vec![0, 1, 1]),
        Integral(vec![2, 0, 1]),
        Integral(vec![1, 1, 0]),
    ];
    let epsilon = Rational::from((1, 100));
    let options = FlowOptions::default();
    let backend = CountedBackend {
        backend: RustRedBackend::default(),
        calls: AtomicUsize::new(0),
    };
    let prepared = Arc::new(AtomicUsize::new(0));
    let counter = prepared.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if matches!(event, Progress::Prepared { .. }) {
                counter.fetch_add(1, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let p = Precision::decimal(60).unwrap();
    let batched = provider
        .evaluate_many(&family, &targets, &epsilon, p)
        .unwrap();
    assert_eq!(prepared.load(Ordering::Relaxed), 2);
    let batch_calls = backend.calls.load(Ordering::Relaxed);
    assert_eq!(batched[0], batched[3]);
    let scalar_backend = CountedBackend {
        backend: RustRedBackend::default(),
        calls: AtomicUsize::new(0),
    };
    let scalar_context = RunContext::default();
    let scalar = recursive::RecursiveBoundary::new(&scalar_backend, &options, &scalar_context);
    for (target, value) in targets.iter().zip(&batched) {
        let separate = scalar.evaluate(&family, target, &epsilon, p).unwrap();
        assert!(p.close(value, &separate, 30));
    }
    assert!(batch_calls < scalar_backend.calls.load(Ordering::Relaxed));
    let dimension = Rational::from(4) - &epsilon * &Rational::from(2);
    for (value, (a, ma, b, mb)) in
        batched
            .iter()
            .zip([(1, 3, 1, 2), (1, 2, 1, 1), (2, 3, 1, 1), (1, 3, 1, 2)])
    {
        let expected = p.mul(
            &vacuum::tadpole(a, &p.i(ma), &dimension, p).unwrap(),
            &vacuum::tadpole(b, &p.i(mb), &dimension, p).unwrap(),
        );
        assert!(p.close(value, &expected, 30));
    }
    let high = Precision::decimal(80).unwrap();
    let refined = provider
        .evaluate_many(&family, &targets, &epsilon, high)
        .unwrap();
    for (a, b) in batched.iter().zip(&refined) {
        assert!(high.close(a, b, 30));
    }
    assert_eq!(prepared.load(Ordering::Relaxed), 2);
    assert_eq!(backend.calls.load(Ordering::Relaxed), batch_calls);
    // Canonical scalar memo entries remain reusable after a batch, including
    // equivalent denominator permutations with a differently ordered family.
    let mut permuted = family.clone();
    permuted.propagators.swap(0, 1);
    let value = provider
        .evaluate(&permuted, &Integral(vec![0, 2, 1]), &epsilon, p)
        .unwrap();
    assert_eq!(value, batched[2]);
    assert_eq!(backend.calls.load(Ordering::Relaxed), batch_calls);
}

#[test]
fn connected_boundary_batch_preserves_recursive_siblings_and_refinement() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    let mut family = sunset();
    family.propagators[0].constant = Atom::num(-1);
    family.propagators[1].constant = Atom::num(-1);
    let targets = [
        Integral(vec![1, 1, 1]),
        Integral(vec![2, 1, 1]),
        Integral(vec![1, 2, 1]),
    ];
    let epsilon = Rational::from((1, 100));
    let backend = RustRedBackend::default();
    let options = FlowOptions::default();
    let preparations = Arc::new(AtomicUsize::new(0));
    let counter = preparations.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if matches!(event, Progress::Prepared { .. }) {
                counter.fetch_add(1, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let p = Precision::decimal(60).unwrap();
    let values = provider
        .evaluate_many(&family, &targets, &epsilon, p)
        .unwrap();
    assert_eq!(preparations.load(Ordering::Relaxed), 1);
    let scalar_context = RunContext::default();
    let scalar = recursive::RecursiveBoundary::new(&backend, &options, &scalar_context);
    for (target, value) in targets.iter().zip(&values) {
        let separate = scalar.evaluate(&family, target, &epsilon, p).unwrap();
        assert!(p.close(value, &separate, 20));
    }
    // An independent Schwinger-parameter beta integral gives I111. The
    // simultaneous mass derivative and equal-mass symmetry then give
    // I211 = I121 = (D-3)/2 * I111, at unit mass squared.
    let expected = |precision: Precision| {
        let e = precision.rational(&epsilon);
        let gamma = precision.gamma_real(&e.re).unwrap();
        let base = precision.div(
            &precision.powi(&gamma, 2),
            &precision.mul(
                &precision.sub(&precision.i(1), &e),
                &precision.sub(&precision.i(1), &precision.scale(&e, 2, 1)),
            ),
        );
        let raised = precision.mul(
            &base,
            &precision.sub(&precision.rational(&Rational::from((1, 2))), &e),
        );
        [base, raised.clone(), raised]
    };
    for (value, exact) in values.iter().zip(expected(p)) {
        assert!(p.close(value, &exact, 20));
    }
    // A fresh provider prevents reuse of the original numeric or symbolic
    // memo, while precision and expansion order both increase.
    let refined_options = FlowOptions {
        digits: 30,
        guard_digits: 50,
        series_order: 112,
        ..options
    };
    let refined_context = RunContext::default();
    let refined_provider =
        recursive::RecursiveBoundary::new(&backend, &refined_options, &refined_context);
    let high = Precision::decimal(80).unwrap();
    let refined = refined_provider
        .evaluate_many(&family, &targets, &epsilon, high)
        .unwrap();
    for ((value, refined), exact) in values.iter().zip(&refined).zip(expected(high)) {
        assert!(high.close(value, refined, 20));
        assert!(high.close(refined, &exact, 30));
    }
}

#[test]
fn boundary_batches_keep_terminal_providers_scaleless_values_and_duplicates() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct Terminal(AtomicUsize);
    impl recursive::TerminalProvider for Terminal {
        fn evaluate(
            &self,
            _: &IntegralFamily,
            target: &Integral,
            _: &Rational,
            p: Precision,
        ) -> Result<Option<ComplexFloat>> {
            self.0.fetch_add(1, Ordering::Relaxed);
            Ok((target.0 == [2]).then(|| p.i(42)))
        }
    }
    let mut family = sunset();
    family.loops.truncate(1);
    family.propagators.truncate(1);
    family.propagators[0].scalar_products = vec![Atom::num(1)];
    family.physical_propagators = 1;
    let terminal = Terminal(AtomicUsize::new(0));
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let provider =
        recursive::RecursiveBoundary::new(&backend, &options, &context).with_terminal(&terminal);
    let p = Precision::decimal(60).unwrap();
    let epsilon = Rational::from((1, 100));
    let targets = [
        Integral(vec![1]),
        Integral(vec![2]),
        Integral(vec![0]),
        Integral(vec![2]),
    ];
    let values = provider
        .evaluate_many(&family, &targets, &epsilon, p)
        .unwrap();
    assert_eq!(terminal.0.load(Ordering::Relaxed), 3);
    assert_eq!(values[1], p.i(42));
    assert_eq!(values[1], values[3]);
    assert_eq!(values[2], p.zero());
    let expected = vacuum::tadpole(1, &p.i(3), &Rational::from((199, 50)), p).unwrap();
    assert!(p.close(&values[0], &expected, 30));
    assert!(
        provider
            .evaluate_many(&family, &[], &epsilon, p)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        provider
            .evaluate_many(&family, &targets, &epsilon, p)
            .unwrap(),
        values
    );
    assert_eq!(terminal.0.load(Ordering::Relaxed), 3);
}

#[test]
fn tadpole_only_policy_recurses_when_custom_terminal_declines_sunset() {
    use recursive::RecursiveTerminalPolicy;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    struct Declines(AtomicUsize);
    impl recursive::TerminalProvider for Declines {
        fn evaluate(
            &self,
            _: &IntegralFamily,
            _: &Integral,
            _: &Rational,
            _: Precision,
        ) -> Result<Option<ComplexFloat>> {
            self.0.fetch_add(1, Ordering::Relaxed);
            Ok(None)
        }
    }
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let preparations = Arc::new(AtomicUsize::new(0));
    let observed = preparations.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if matches!(event, Progress::Prepared { .. }) {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let options = FlowOptions::default();
    let terminal = Declines(AtomicUsize::new(0));
    let family = sunset();
    let target = Integral(vec![1, 1, 1]);
    let epsilon = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let provider =
        recursive::RecursiveBoundary::new(&backend, &options, &context).with_terminal(&terminal);
    let analytic = provider.evaluate(&family, &target, &epsilon, p).unwrap();
    assert_eq!(preparations.load(Ordering::Relaxed), 0);
    assert_eq!(terminal.0.load(Ordering::Relaxed), 1);

    // Reusing a provider after changing its policy must not reuse the sunset
    // terminal's memoized value. The loop count decreases under radial peeling;
    // the remaining one-loop problem runs through ordinary AMF and tadpole seeds.
    let provider = provider.with_terminal_policy(RecursiveTerminalPolicy::TadpolesOnly);
    let recursive = provider.evaluate(&family, &target, &epsilon, p).unwrap();
    assert!(preparations.load(Ordering::Relaxed) > 0);
    assert!(terminal.0.load(Ordering::Relaxed) > 1);
    assert!(
        p.close(&recursive, &analytic, 20),
        "{recursive} != {analytic}"
    );

    let mut tadpole = family;
    tadpole.loops.truncate(1);
    tadpole.propagators.truncate(1);
    tadpole.propagators[0].scalar_products = vec![Atom::num(1)];
    tadpole.physical_propagators = 1;
    let before = preparations.load(Ordering::Relaxed);
    let seed = provider
        .evaluate(&tadpole, &Integral(vec![2]), &epsilon, p)
        .unwrap();
    let expected = vacuum::tadpole(2, &p.i(3), &Rational::from((19, 5)), p).unwrap();
    assert!(p.close(&seed, &expected, 40));
    assert_eq!(preparations.load(Ordering::Relaxed), before);
}

#[test]
fn replacing_a_terminal_provider_does_not_reuse_old_values() {
    struct Constant(i64);
    impl recursive::TerminalProvider for Constant {
        fn evaluate(
            &self,
            _: &IntegralFamily,
            _: &Integral,
            _: &Rational,
            p: Precision,
        ) -> Result<Option<ComplexFloat>> {
            Ok(Some(p.i(self.0)))
        }
    }
    let backend = RustRedBackend::default();
    let options = FlowOptions::default();
    let context = RunContext::default();
    let first = Constant(2);
    let second = Constant(3);
    let p = Precision::decimal(60).unwrap();
    let family = sunset();
    let target = Integral(vec![1, 1, 1]);
    let epsilon = Rational::from((1, 10));
    let provider =
        recursive::RecursiveBoundary::new(&backend, &options, &context).with_terminal(&first);
    assert_eq!(
        provider.evaluate(&family, &target, &epsilon, p).unwrap(),
        p.i(2)
    );
    let provider = provider.with_terminal(&second);
    assert_eq!(
        provider.evaluate(&family, &target, &epsilon, p).unwrap(),
        p.i(3)
    );
}
