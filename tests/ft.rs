use symbolica::prelude::*;
use symbolica_amflow::*;
fn bubble() -> IntegralFamily {
    IntegralFamily {
        name: "ft_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-2)]],
        propagators: vec![
            Propagator {
                constant: Atom::num(-1),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            Propagator {
                constant: Atom::num(-5),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}
#[test]
fn gaussian_tensor_terminal_includes_shift_and_wick_contractions() {
    let mut family = bubble();
    family.propagators[1] = family.propagators[0].clone();
    family.propagators[0] = Propagator {
        constant: Atom::num(-3),
        scalar_products: vec![Atom::num(1), Atom::num(1)],
    };
    family.physical_propagators = 1;
    let eps = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let d = Rational::from((19, 5));
    let m = p.scale(&p.i(1), 5, 2);
    let tad = |n| vacuum::tadpole(n, &m, &d, p).unwrap();
    let first = gaussian::terminal(&family, &Integral(vec![3, -1]), &eps, p).unwrap();
    assert!(p.close(&first, &p.add(&tad(2), &tad(3)), 40));
    let second = gaussian::terminal(&family, &Integral(vec![3, -2]), &eps, p).unwrap();
    let expected = p.add(
        &p.add(&tad(1), &p.add(&p.scale(&tad(2), 2, 1), &tad(3))),
        &p.scale(&p.add(&tad(2), &p.mul(&m, &tad(3))), -10, 19),
    );
    assert!(p.close(&second, &expected, 40), "{second} != {expected}");
}
#[test]
fn recursive_ft_agrees_with_direct_bubble_parameter_system() {
    let family = bubble();
    let target = Integral(vec![1, 1]);
    let eps = Rational::from((1, 10));
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let p = Precision::decimal(60).unwrap();
    let direct = ft::evaluate_bubble(&family, &target, &eps, &options, &context).unwrap();
    for (skip_reduction, refine_basis) in
        [(false, false), (true, false), (false, true), (true, true)]
    {
        let options = FlowOptions {
            skip_reduction,
            refine_basis,
            ..options.clone()
        };
        let general = ft::FtEvaluator::new(&backend, &context)
            .evaluate(&family, &target, &eps, &options)
            .unwrap();
        assert!(p.close(&general, &direct, 20), "{general} != {direct}");
    }
}

#[test]
fn stricter_ft_accuracy_recomputes_at_equal_precision_and_order() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    let mut family = bubble();
    family.propagators[0].constant = Atom::new();
    family.propagators[1].constant = Atom::num(-2);
    let target = Integral(vec![1, 1]);
    let epsilon = Rational::from((1, 10));
    let loose = FlowOptions {
        digits: 10,
        guard_digits: 30,
        series_order: 64,
        ..Default::default()
    };
    let strict = FlowOptions {
        digits: 35,
        guard_digits: 5,
        ..loose.clone()
    };
    assert_eq!(
        Precision::decimal(loose.digits + loose.guard_digits)
            .unwrap()
            .bits,
        Precision::decimal(strict.digits + strict.guard_digits)
            .unwrap()
            .bits
    );
    let steps = Arc::new(AtomicUsize::new(0));
    let observed = steps.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if matches!(event, Progress::Step { .. }) {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let backend = RustRedBackend::default();
    let evaluator = ft::FtEvaluator::new(&backend, &context);
    evaluator
        .evaluate(&family, &target, &epsilon, &loose)
        .unwrap();
    let initial_steps = steps.load(Ordering::Relaxed);
    assert!(initial_steps > 0);
    let value = evaluator
        .evaluate(&family, &target, &epsilon, &strict)
        .unwrap();
    let strict_steps = steps.load(Ordering::Relaxed);
    assert!(
        strict_steps > initial_steps,
        "a stricter requested tolerance must run numerical continuation again"
    );

    // Independent massless-bubble formula at higher arithmetic precision.
    let p = Precision::decimal(70).unwrap();
    let gamma = |argument: Rational| p.gamma_real(&p.rational(&argument).re).unwrap();
    let expected = p.mul(
        &p.div(
            &p.mul(
                &gamma(epsilon.clone()),
                &p.powi(&gamma(Rational::from(1) - &epsilon), 2),
            ),
            &gamma(Rational::from(2) - &epsilon * &Rational::from(2)),
        ),
        &p.pow(&p.i(2), &p.rational(&(-epsilon.clone()))),
    );
    assert!(p.close(&value, &expected, strict.digits));
    let repeated = evaluator
        .evaluate(&family, &target, &epsilon, &strict)
        .unwrap();
    assert_eq!(repeated, value);
    assert_eq!(steps.load(Ordering::Relaxed), strict_steps);
}

#[test]
fn recursive_ft_two_loop_sunset_matches_schwinger_parameters() {
    let family = IntegralFamily {
        name: "ft_sunset".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::num(1), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::num(1), &[]).unwrap(),
            Propagator::quadratic(&[1, 1], &[], Atom::new(), &[]).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("eps"),
        dimension: 4,
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = FlowOptions::default();
    let eps = Rational::from((3, 4));
    let p = Precision::decimal(60).unwrap();
    let answer = ft::FtEvaluator::new(&backend, &context)
        .evaluate(&family, &Integral(vec![1, 1, 1]), &eps, &options)
        .unwrap();
    let gamma = |n, d| {
        p.gamma_real(&p.rational(&Rational::from((n, d))).re)
            .unwrap()
    };
    let expected = p.neg(&p.div(
        &p.mul(&gamma(1, 2), &p.powi(&gamma(3, 4), 2)),
        &p.scale(&gamma(3, 2), 1, 4),
    ));
    assert!(p.close(&answer, &expected, 20), "{answer} != {expected}");
}

#[test]
fn high_level_ft_laurent_coefficients_match_massless_bubble() {
    let mut family = bubble();
    family.external_gram[0][0] = Atom::num(-1);
    family.propagators[0].constant = Atom::new();
    family.propagators[1].constant = Atom::num(-1);
    let options = FlowOptions {
        recursion: RecursionMode::Ft,
        ..Default::default()
    };
    let result = solve_integrals(
        &family,
        &[Integral(vec![1, 1])],
        &KinematicPoint::default(),
        0,
        &options,
        &RustRedBackend::default(),
        &RunContext::default(),
    )
    .unwrap()
    .remove(0);
    assert_eq!(result.verified_digits, Some(20));
    assert!(result.validation_samples > 0);
    assert!(result.refinements > 0);
    let p = Precision::decimal(70).unwrap();
    // Gamma(eps) Gamma(1-eps)^2 / Gamma(2-2eps) at -s=1.
    let euler_gamma = p
        .parse(
            "0.5772156649015328606065120900824024310421593359399235988057672348848677",
            "0",
        )
        .unwrap();
    for (power, expected) in [
        (-2, p.zero()),
        (-1, p.i(1)),
        (0, p.sub(&p.i(2), &euler_gamma)),
    ] {
        assert!(
            p.close(&result.coefficients[&power], &expected, 20),
            "epsilon^{power}: {} != {expected}",
            result.coefficients[&power]
        );
        assert!(result.comparison_errors.contains_key(&power));
    }
}

#[test]
fn ft_symbolic_system_restarts_without_reducer_and_separates_options() {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    struct RestartBackend {
        allowed: AtomicBool,
        calls: AtomicUsize,
    }
    impl ReductionBackend for RestartBackend {
        fn identity(&self) -> String {
            format!("ft-restart-test:{}", RustRedBackend::default().identity())
        }
        fn reduce(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            if !self.allowed.load(Ordering::Relaxed) {
                return Err(Error::Reduction("FT restart reducer disabled".into()));
            }
            RustRedBackend::default().reduce(family, targets, context)
        }
    }
    let directory = std::env::temp_dir().join(format!(
        "symbolica-amflow-ft-restart-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&directory);
    let options = FlowOptions {
        recursion: RecursionMode::Ft,
        cache_directory: Some(directory.clone()),
        ..Default::default()
    };
    let backend = RestartBackend {
        allowed: AtomicBool::new(true),
        calls: AtomicUsize::new(0),
    };
    let context = RunContext::default();
    let family = bubble();
    let target = Integral(vec![1, 1]);
    let epsilon = Rational::from((1, 10));
    let first = ft::FtEvaluator::new(&backend, &context)
        .evaluate(&family, &target, &epsilon, &options)
        .unwrap();
    let calls = backend.calls.load(Ordering::Relaxed);
    assert!(calls > 0);
    assert!(std::fs::read_dir(&directory).unwrap().any(|entry| {
        entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with("system-ft-")
    }));
    backend.allowed.store(false, Ordering::Relaxed);
    // A fresh evaluator has no in-memory systems or values.
    let restored = ft::FtEvaluator::new(&backend, &context)
        .evaluate(&family, &target, &epsilon, &options)
        .unwrap();
    assert_eq!(backend.calls.load(Ordering::Relaxed), calls);
    assert!(Precision::decimal(60).unwrap().close(&first, &restored, 35));
    let changed = FlowOptions {
        skip_reduction: true,
        ..options
    };
    assert!(matches!(
        ft::FtEvaluator::new(&backend, &context)
            .evaluate(&family, &target, &epsilon, &changed),
        Err(Error::Reduction(message)) if message == "FT restart reducer disabled"
    ));
    assert!(backend.calls.load(Ordering::Relaxed) > calls);
    std::fs::remove_dir_all(directory).unwrap();
}
