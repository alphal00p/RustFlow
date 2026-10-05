//! Single-mass loop peeling and public automatic evaluation regressions.
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::prelude::*;
use symbolica_amflow::*;

fn six_line() -> IntegralFamily {
    IntegralFamily {
        name: "single_mass_connected_six_line".into(),
        loops: vec!["q".into(), "l".into(), "k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: [
            ([1, 0, 0], 3),
            ([0, 1, 0], 0),
            ([-1, 1, 0], 0),
            ([0, 0, 1], 0),
            ([-1, 0, 1], 0),
            ([0, 1, -1], 0),
        ]
        .into_iter()
        .map(|(v, m)| Propagator::quadratic(&v, &[], Atom::num(m), &[]).unwrap())
        .collect(),
        physical_propagators: 6,
        epsilon: symbol!("single_mass_recursive_test::epsilon"),
        dimension: 4,
    }
}
struct CountLoops {
    inner: RustRedBackend,
    calls: AtomicUsize,
    max_loops: AtomicUsize,
    peeled_children: AtomicUsize,
}
impl CountLoops {
    fn new() -> Self {
        Self {
            inner: RustRedBackend {
                max_depth: 3,
                max_targets: 65536,
                max_sector_batch: 32,
                max_backward_frontier: 1024,
                parametric_rules: true,
                symmetry_rules: true,
                ..Default::default()
            },
            calls: AtomicUsize::new(0),
            max_loops: AtomicUsize::new(0),
            peeled_children: AtomicUsize::new(0),
        }
    }
    fn record(&self, family: &IntegralFamily) {
        if family.name.ends_with("_single_mass_child") {
            assert_eq!(family.loops.len(), 2);
            self.peeled_children.fetch_add(1, Ordering::Relaxed);
        }
        self.calls.fetch_add(1, Ordering::Relaxed);
        self.max_loops
            .fetch_max(family.loops.len(), Ordering::Relaxed);
    }
}
impl ReductionBackend for CountLoops {
    fn identity(&self) -> String {
        format!("single-mass-loop-audit:{}", self.inner.identity())
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<reduction::Reduction> {
        self.record(family);
        self.inner.reduce(family, targets, context)
    }
    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<reduction::Reduction> {
        self.record(family);
        self.inner
            .reduce_at_epsilon(family, targets, epsilon, context)
    }
}

#[test]
fn single_mass_six_line_reuses_lower_loop_child_and_matches_parent_ft_with_refinement() {
    let family = six_line();
    let targets = [Integral(vec![1; 6]), Integral(vec![2, 1, 1, 1, 1, 1])];
    for epsilon in [Rational::from((1, 10)), Rational::from((1, 11))] {
        let mut previous: Option<Vec<ComplexFloat>> = None;
        for (digits, order) in [(60, 80), (80, 112)] {
            let p = Precision::decimal(digits).unwrap();
            let options = FlowOptions {
                guard_digits: digits - 20,
                series_order: order,
                max_steps: 1000,
                workers: 1,
                ..Default::default()
            };
            let context = RunContext::default();
            let backend = CountLoops::new();
            let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
            let scalar = provider
                .evaluate(&family, &targets[0], &epsilon, p)
                .unwrap();
            let calls = backend.calls.load(Ordering::Relaxed);
            assert!(calls > 0);
            assert_eq!(
                backend.max_loops.load(Ordering::Relaxed),
                2,
                "the new rule must remove one loop before any backend call"
            );
            let raised = provider
                .evaluate(&family, &targets[1], &epsilon, p)
                .unwrap();
            assert_eq!(
                backend.calls.load(Ordering::Relaxed),
                calls,
                "massive-power variants share the same child value"
            );
            assert!(p.close(
                &raised,
                &p.mul(&p.rational(&(-epsilon.clone())), &scalar),
                20
            ));
            let baseline_backend = CountLoops::new();
            let ft = ft::FtEvaluator::new(&baseline_backend, &context);
            let values = vec![scalar, raised];
            for (target, value) in targets.iter().zip(&values) {
                let independent = ft.evaluate(&family, target, &epsilon, &options).unwrap();
                assert!(p.close(value, &independent, 20), "{value} != {independent}");
            }
            if let Some(old) = previous {
                for (low, high) in old.iter().zip(&values) {
                    assert!(p.close(low, high, 20));
                }
            }
            previous = Some(values);
        }
    }
}

#[test]
fn complex_mass_boundaries_reuse_the_real_child_and_respect_cut_sides() {
    use symbolica::domains::float::Complex;
    let original = six_line();
    let epsilon = Rational::from((1, 10));
    let targets = [Integral(vec![1; 6]), Integral(vec![2, 1, 1, 1, 1, 1])];
    let tiny = Rational::one() / Rational::from(10).pow(30);
    let masses = [
        Complex::new(Rational::from(3), Rational::from((4, 3))),
        Complex::new(Rational::from(3), Rational::from((-4, 3))),
        Complex::new(Rational::from(-3), tiny.clone()),
        Complex::new(Rational::from(-3), -tiny),
    ];
    let mut previous: Option<Vec<Vec<ComplexFloat>>> = None;
    for (working, order) in [(60, 80), (80, 112)] {
        let p = Precision::decimal(working).unwrap();
        let options = FlowOptions {
            guard_digits: working - 20,
            series_order: order,
            ..Default::default()
        };
        let backend = CountLoops::new();
        let context = RunContext::default();
        let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
        let real = provider
            .evaluate(&original, &targets[0], &epsilon, p)
            .unwrap();
        let independent = ft::FtEvaluator::new(&CountLoops::new(), &context)
            .evaluate(&original, &targets[0], &epsilon, &options)
            .unwrap();
        assert!(p.close(&real, &independent, 20));
        let calls = backend.calls.load(Ordering::Relaxed);
        assert!(calls > 0);
        let mut results = Vec::new();
        for mass in &masses {
            let mut family = original.clone();
            family.propagators[0].constant = -Atom::num(mass.clone());
            let values = provider
                .evaluate_many(&family, &targets, &epsilon, p)
                .unwrap();
            // Homogeneity and the mass derivative are separate exact relations
            // on the analytically continued parent, including both cut sides.
            let numeric_mass = p
                .eval(&Atom::num(mass.clone()), &Default::default())
                .unwrap();
            let degree = p.rational(&(-Rational::from(3) * &epsilon));
            let expected = p.mul(
                &independent,
                &p.pow(&p.div(&numeric_mass, &p.i(3)), &degree),
            );
            assert!(p.close(&values[0], &expected, 20));
            assert!(p.close(
                &values[1],
                &p.mul(&p.div(&degree, &numeric_mass), &values[0]),
                20
            ));
            assert_eq!(
                backend.calls.load(Ordering::Relaxed),
                calls,
                "changing the mass or its power reuses the same massless child"
            );
            results.push(values);
        }
        assert_eq!(backend.max_loops.load(Ordering::Relaxed), 2);
        for pair in results.chunks_exact(2) {
            for (upper, lower) in pair[0].iter().zip(&pair[1]) {
                let conjugate = ComplexFloat::new(upper.re.clone(), -upper.im.clone());
                assert!(p.close(&conjugate, lower, 20));
            }
        }
        assert!(
            p.norm(&p.sub(&results[2][0], &results[3][0])) > p.real(1),
            "opposite sides of the mass cut must retain their discontinuity"
        );
        if let Some(old) = previous {
            for (low, high) in old.iter().flatten().zip(results.iter().flatten()) {
                assert!(p.close(low, high, 20));
            }
        }
        previous = Some(results);
    }
}

#[test]
fn public_complex_mass_amf_matches_homogeneity_and_precision_refinement() {
    use symbolica::domains::float::Complex;
    let original = six_line();
    let epsilon = Rational::from((1, 10));
    let targets = [Integral(vec![1; 6]), Integral(vec![2, 1, 1, 1, 1, 1])];
    let mut previous: Option<Vec<Vec<ComplexFloat>>> = None;
    for (working, order) in [(60, 80), (80, 112)] {
        let p = Precision::decimal(working).unwrap();
        let options = FlowOptions {
            guard_digits: working - 20,
            series_order: order,
            ..Default::default()
        };
        let context = RunContext::default();
        let independent = ft::FtEvaluator::new(&CountLoops::new(), &context)
            .evaluate(&original, &targets[0], &epsilon, &options)
            .unwrap();
        let mut results = Vec::new();
        for imaginary_part in [-4, 4] {
            let mass = Atom::num(Complex::new(
                Rational::from(3),
                Rational::from((imaginary_part, 3)),
            ));
            let mut family = original.clone();
            family.propagators[0].constant = -&mass;
            let backend = CountLoops::new();
            let prepared = PreparedFlow::new(
                &family,
                &targets,
                &KinematicPoint::default(),
                &backend,
                &options,
                &context,
            )
            .unwrap();
            let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
            let values = evaluate_samples(
                &prepared,
                std::slice::from_ref(&epsilon),
                &options,
                &provider,
                &context,
            )
            .unwrap()
            .remove(0);
            let mass = p.eval(&mass, &Default::default()).unwrap();
            let degree = p.rational(&(-Rational::from(3) * &epsilon));
            let expected = p.mul(&independent, &p.pow(&p.div(&mass, &p.i(3)), &degree));
            assert!(
                p.close(&values[0], &expected, 20),
                "{values:?} != {expected}"
            );
            assert!(p.close(&values[1], &p.mul(&p.div(&degree, &mass), &expected), 20));
            assert_eq!(backend.max_loops.load(Ordering::Relaxed), 3);
            assert!(backend.peeled_children.load(Ordering::Relaxed) > 0);
            results.push(values);
        }
        if let Some(old) = previous {
            for (low, high) in old.iter().flatten().zip(results.iter().flatten()) {
                assert!(p.close(low, high, 20));
            }
        }
        previous = Some(results);
    }
}

#[test]
fn single_mass_peeling_preserves_custom_provider_and_scaleless_precedence() {
    struct Supplied;
    impl recursive::TerminalProvider for Supplied {
        fn evaluate(
            &self,
            _: &IntegralFamily,
            _: &Integral,
            _: &Rational,
            p: Precision,
        ) -> Result<Option<ComplexFloat>> {
            Ok(Some(p.i(17)))
        }
    }
    let family = six_line();
    let backend = CountLoops::new();
    let options = FlowOptions::default();
    let context = RunContext::default();
    let p = Precision::decimal(60).unwrap();
    // epsilon=0 would make the radial prefactor singular. Neither the supplied
    // value nor a proved-scaleless sector may reach that prefactor.
    let value = recursive::RecursiveBoundary::new(&backend, &options, &context)
        .with_terminal(&Supplied)
        .evaluate(&family, &Integral(vec![1; 6]), &Rational::zero(), p)
        .unwrap();
    assert_eq!(value, p.i(17));
    let value = recursive::RecursiveBoundary::new(&backend, &options, &context)
        .evaluate(
            &family,
            &Integral(vec![1, 1, 0, 1, 0, 1]),
            &Rational::zero(),
            p,
        )
        .unwrap();
    assert_eq!(value, p.zero());
    assert_eq!(backend.calls.load(Ordering::Relaxed), 0);
}

#[test]
#[ignore = "full 20-digit Laurent acceptance; run in release mode with a 600-second external limit"]
fn single_mass_public_laurent_matches_independent_ft_and_mass_derivative() {
    let family = six_line();
    let targets = [Integral(vec![1; 6]), Integral(vec![2, 1, 1, 1, 1, 1])];
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "symbolica-amflow-single-mass-laurent-{}-{stamp}",
        std::process::id()
    ));
    std::fs::create_dir(&directory).unwrap();
    let p = Precision::decimal(100).unwrap();
    let mut previous: Option<Vec<LaurentExpansion>> = None;
    for (name, recursion, working, order) in [
        ("amf", RecursionMode::Amf, 60, 80),
        ("ft", RecursionMode::Ft, 80, 112),
    ] {
        let mut audited = CountLoops::new();
        audited.inner.max_exact_frontier = 512;
        audited.inner.native_workers = 1;
        audited.inner.checkpoints = Some(directory.join(name).join("native"));
        let backend = cache::CachedBackend {
            backend: audited,
            directory: directory.join(name).join("reductions"),
        };
        let options = FlowOptions {
            digits: 20,
            guard_digits: working - 20,
            series_order: order,
            max_steps: 1000,
            max_precision_attempts: 3,
            workers: 1,
            recursion,
            cache_directory: Some(directory.join(name).join("systems")),
            ..Default::default()
        };
        let result = solve_integrals(
            &family,
            &targets,
            &KinematicPoint::default(),
            0,
            &options,
            &backend,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(result.len(), 2);
        if name == "amf" {
            assert!(backend.backend.peeled_children.load(Ordering::Relaxed) > 0);
        }
        for fit in &result {
            assert_eq!(fit.verified_digits, Some(20));
            assert!(fit.refinements > 0 && fit.validation_samples > 0);
            assert_eq!(fit.coefficients.len(), 7);
            assert_eq!(fit.comparison_errors.len(), 7);
        }
        for power in -6..=0 {
            let expected = result[0]
                .coefficients
                .get(&(power - 1))
                .map_or_else(|| p.zero(), |c| p.neg(c));
            assert!(p.close(&result[1].coefficients[&power], &expected, 20));
            if power <= -2 {
                assert!(p.norm(&result[0].coefficients[&power]) < p.tolerance(20));
            }
            if power < 0 {
                assert!(p.norm(&result[1].coefficients[&power]) < p.tolerance(20));
            }
            if let Some(old) = &previous {
                for (a, b) in old.iter().zip(&result) {
                    assert!(p.close(&a.coefficients[&power], &b.coefficients[&power], 20));
                }
            }
        }
        previous = Some(result);
    }
    std::fs::remove_dir_all(directory).unwrap();
}
