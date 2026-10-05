//! Automatic AMF must not silently change the principal single-mass sheet.
use std::sync::atomic::{AtomicBool, Ordering};
use symbolica::prelude::*;
use symbolica_amflow::*;

fn mass(real: i64, imaginary: i64) -> Atom {
    Atom::num(Complex::new(
        Rational::from(real),
        Rational::from((imaginary, 3)),
    ))
}

fn tadpole(mass: Atom) -> IntegralFamily {
    IntegralFamily {
        name: "principal_mass_sheet_tadpole".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        // 4*(k²-M²) also checks positive nonunit quadratic normalization.
        propagators: vec![Propagator::quadratic(&[2], &[], mass * 4, &[]).unwrap()],
        physical_propagators: 1,
        epsilon: symbol!("principal_mass_sheet::epsilon"),
        dimension: 4,
    }
}

fn six_line(mass: Atom) -> IntegralFamily {
    IntegralFamily {
        name: "principal_mass_sheet_six_line".into(),
        loops: vec!["q".into(), "l".into(), "k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: [
            ([1, 0, 0], mass),
            ([0, 1, 0], Atom::new()),
            ([-1, 1, 0], Atom::new()),
            ([0, 0, 1], Atom::new()),
            ([-1, 0, 1], Atom::new()),
            ([0, 1, -1], Atom::new()),
        ]
        .into_iter()
        .map(|(routing, mass)| Propagator::quadratic(&routing, &[], mass, &[]).unwrap())
        .collect(),
        physical_propagators: 6,
        epsilon: symbol!("principal_mass_sheet_six::epsilon"),
        dimension: 4,
    }
}

fn unsupported_sheet<T>(result: Result<T>) {
    match result {
        Err(Error::Unsupported(message)) => {
            assert!(message.contains("principal single-mass sheet"), "{message}")
        }
        _ => panic!("expected the typed principal-mass-sheet guard"),
    }
}

struct SwitchBackend {
    reject: AtomicBool,
    inner: RustRedBackend,
}
impl SwitchBackend {
    fn new(reject: bool) -> Self {
        Self {
            reject: AtomicBool::new(reject),
            inner: RustRedBackend::default(),
        }
    }
}
impl ReductionBackend for SwitchBackend {
    fn identity(&self) -> String {
        self.inner.identity()
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<reduction::Reduction> {
        assert!(
            !self.reject.load(Ordering::Relaxed),
            "the mass-sheet guard or prepared cache must precede reduction"
        );
        self.inner.reduce(family, targets, context)
    }
}

#[test]
fn incompatible_quadrants_and_axis_reject_before_cache_and_backend() {
    let backend = SwitchBackend::new(true);
    let context = RunContext::default();
    // This is deliberately a regular file. A cache access would fail with an
    // I/O error instead of the required prior mass-sheet error.
    let cache = std::env::temp_dir().join(format!(
        "amflow-mass-sheet-cache-file-{}",
        std::process::id()
    ));
    std::fs::write(&cache, b"not a cache directory").unwrap();
    for real in [-3, 0] {
        for (imaginary, prescription) in [(4, Prescription::PlusI0), (-4, Prescription::MinusI0)] {
            let options = FlowOptions {
                prescription,
                cache_directory: Some(cache.clone()),
                ..Default::default()
            };
            for family in [
                tadpole(mass(real, imaginary)),
                six_line(mass(real, imaginary)),
            ] {
                let target = Integral(vec![1; family.propagators.len()]);
                unsupported_sheet(PreparedFlow::new(
                    &family,
                    std::slice::from_ref(&target),
                    &KinematicPoint::default(),
                    &backend,
                    &options,
                    &context,
                ));
                unsupported_sheet(PreparedFlow::new_at_epsilon(
                    &family,
                    &[target],
                    &KinematicPoint::default(),
                    &backend,
                    &options,
                    &context,
                    &Rational::from((1, 10)),
                ));
            }
        }
    }
    std::fs::remove_file(cache).unwrap();
}

#[test]
fn exact_kinematic_substitution_is_checked_before_reduction() {
    let m = symbol!("principal_mass_sheet::m_squared");
    let family = tadpole(Atom::var(m));
    let point = KinematicPoint(std::collections::BTreeMap::from([(
        Atom::var(m),
        mass(-3, 4),
    )]));
    unsupported_sheet(PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &point,
        &SwitchBackend::new(true),
        &FlowOptions::default(),
        &RunContext::default(),
    ));
}

#[test]
fn certified_scaleless_sector_has_no_mass_sheet_to_reject() {
    let family = IntegralFamily {
        name: "principal_mass_sheet_scaleless_sector".into(),
        loops: vec!["q".into(), "l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], mass(-3, 4), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::new(), &[]).unwrap(),
            Propagator::quadratic(&[1, -1], &[], Atom::new(), &[]).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("principal_mass_sheet_scaleless::epsilon"),
        dimension: 4,
    };
    let target = Integral(vec![1, 1, 0]);
    assert!(recursive::scaleless(&family, &target).unwrap());
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let epsilon = Rational::from((1, 10));
    for prescription in [Prescription::PlusI0, Prescription::MinusI0] {
        let options = FlowOptions {
            prescription,
            ..Default::default()
        };
        for sampled in [false, true] {
            let prepared = if sampled {
                PreparedFlow::new_at_epsilon(
                    &family,
                    std::slice::from_ref(&target),
                    &KinematicPoint::default(),
                    &backend,
                    &options,
                    &context,
                    &epsilon,
                )
            } else {
                PreparedFlow::new(
                    &family,
                    std::slice::from_ref(&target),
                    &KinematicPoint::default(),
                    &backend,
                    &options,
                    &context,
                )
            }
            .unwrap();
            assert!(prepared.reduced.basis.is_empty());
            let value = prepared
                .evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context)
                .unwrap();
            assert_eq!(value, vec![Precision::decimal(60).unwrap().zero()]);
        }
    }
}

#[test]
fn cache_hit_retains_constraint_and_changed_evaluation_prescription_rejects() {
    let family = tadpole(mass(-3, 4));
    let backend = SwitchBackend::new(false);
    let cache = std::env::temp_dir().join(format!(
        "amflow-mass-sheet-cache-directory-{}",
        std::process::id()
    ));
    if cache.exists() {
        std::fs::remove_dir_all(&cache).unwrap();
    }
    let mut options = FlowOptions {
        prescription: Prescription::MinusI0,
        cache_directory: Some(cache.clone()),
        ..Default::default()
    };
    let context = RunContext::default();
    let first = PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    backend.reject.store(true, Ordering::Relaxed);
    let cached = PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    options.prescription = Prescription::PlusI0;
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    for prepared in [first, cached] {
        unsupported_sheet(prepared.evaluate(
            &Rational::from((1, 10)),
            &options,
            &provider,
            &context,
        ));
    }
    unsupported_sheet(PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    ));
    std::fs::remove_dir_all(cache).unwrap();
}

#[test]
fn compatible_tadpole_sides_match_principal_gamma_formula() {
    let epsilon = Rational::from((1, 10));
    let p = Precision::decimal(80).unwrap();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    for (imaginary, prescription) in [(-4, Prescription::PlusI0), (4, Prescription::MinusI0)] {
        let m = mass(-3, imaginary);
        let options = FlowOptions {
            guard_digits: 60,
            series_order: 112,
            prescription,
            ..Default::default()
        };
        let mut family = tadpole(m.clone());
        // The current automatic one-loop infinity provider requires unit
        // normalization; the early guard above separately tests factor four.
        family.propagators[0].scalar_products[0] = Atom::one();
        family.propagators[0].constant = -&m;
        let prepared = PreparedFlow::new(
            &family,
            &[Integral(vec![1])],
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )
        .unwrap();
        let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
        let actual = prepared
            .evaluate(&epsilon, &options, &provider, &context)
            .unwrap();
        let exponent = Rational::one() - &epsilon;
        let expected = p.neg(&p.mul(
            &p.gamma_real(&p.rational(&(-exponent.clone())).re).unwrap(),
            &p.pow(
                &p.eval(&m, &Default::default()).unwrap(),
                &p.rational(&exponent),
            ),
        ));
        assert!(
            p.close(&actual[0], &expected, 20),
            "{:?} != {expected}",
            actual
        );
    }
}

#[test]
fn compatible_six_line_sides_match_principal_homogeneity_with_refinement() {
    let epsilon = Rational::from((1, 10));
    let targets = [Integral(vec![1; 6]), Integral(vec![2, 1, 1, 1, 1, 1])];
    let backend = RustRedBackend {
        max_depth: 3,
        max_targets: 65536,
        max_sector_batch: 32,
        max_backward_frontier: 1024,
        parametric_rules: true,
        symmetry_rules: true,
        ..Default::default()
    };
    let context = RunContext::default();
    let mut previous: Option<Vec<Vec<ComplexFloat>>> = None;
    for (working, order) in [(60, 80), (80, 112)] {
        let p = Precision::decimal(working).unwrap();
        let options = FlowOptions {
            guard_digits: working - 20,
            series_order: order,
            ..Default::default()
        };
        let base = ft::FtEvaluator::new(&backend, &context)
            .evaluate(&six_line(Atom::num(3)), &targets[0], &epsilon, &options)
            .unwrap();
        let mut values = Vec::new();
        for (imaginary, prescription) in [(-4, Prescription::PlusI0), (4, Prescription::MinusI0)] {
            let options = FlowOptions {
                prescription,
                ..options.clone()
            };
            let m = mass(-3, imaginary);
            let family = six_line(m.clone());
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
            let actual = prepared
                .evaluate(&epsilon, &options, &provider, &context)
                .unwrap();
            let m = p.eval(&m, &Default::default()).unwrap();
            let degree = p.rational(&(-Rational::from(3) * &epsilon));
            let expected = p.mul(&base, &p.pow(&p.div(&m, &p.i(3)), &degree));
            assert!(
                p.close(&actual[0], &expected, 20),
                "{actual:?} != {expected}"
            );
            assert!(p.close(&actual[1], &p.mul(&p.div(&degree, &m), &expected), 20));
            values.push(actual);
        }
        if let Some(old) = previous {
            for (a, b) in old.iter().flatten().zip(values.iter().flatten()) {
                assert!(p.close(a, b, 20));
            }
        }
        previous = Some(values);
    }
}
