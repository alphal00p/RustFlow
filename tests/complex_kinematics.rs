use std::collections::BTreeMap;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};
use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn exact_complex_external_invariant_matches_the_analytic_bubble() {
    let s = Atom::var(symbol!("complex_bubble::s"));
    let gram = vec![vec![s.clone()]];
    let family = IntegralFamily {
        name: "complex_external_massless_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        propagators: vec![
            Propagator::quadratic(&[1], &[0], Atom::new(), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1], Atom::new(), &gram).unwrap(),
        ],
        external_gram: gram,
        physical_propagators: 2,
        epsilon: symbol!("complex_bubble::eps"),
        dimension: 4,
    };
    let exact_s = Atom::num(symbolica::domains::float::Complex::new(
        Rational::from(-2),
        Rational::from(1),
    ));
    let point = KinematicPoint(BTreeMap::from([(s, exact_s.clone())]));
    let epsilon = Rational::from((1, 10));
    let options = FlowOptions::default();
    let searches = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&searches);
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if let Progress::SectorReduction { integrals, .. } = event {
                observed.fetch_add(integrals, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let prepared = PreparedFlow::new_at_epsilon(
        &family,
        &[Integral(vec![1, 1])],
        &point,
        &RustRedBackend::default(),
        &options,
        &context,
        &epsilon,
    )
    .unwrap();
    assert_eq!(prepared.family.external_gram[0][0], exact_s);
    assert_eq!(prepared.family.propagators[1].constant, exact_s);
    assert!(searches.load(Ordering::Relaxed) > 0);
    let first = prepared
        .evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context)
        .unwrap();
    let refined = FlowOptions {
        guard_digits: options.guard_digits + 20,
        series_order: options.series_order + 32,
        ..options.clone()
    };
    let second = prepared
        .evaluate(&epsilon, &refined, &boundary::OneLoopBoundary, &context)
        .unwrap();

    // In the upper half-plane of s, the physical massless bubble is
    // Gamma(eps) Gamma(1-eps)^2 / Gamma(2-2eps) * (-s)^(-eps).
    // Here -s=2-i is away from the logarithm cut; the principal branch is
    // reached by continuation from negative real s with the +i0 prescription.
    let p = Precision::decimal(80).unwrap();
    let gamma = |value: Rational| p.gamma_real(&p.rational(&value).re).unwrap();
    let prefactor = p.div(
        &p.mul(
            &gamma(epsilon.clone()),
            &p.powi(&gamma(Rational::from(1) - &epsilon), 2),
        ),
        &gamma(Rational::from(2) - &epsilon * &Rational::from(2)),
    );
    let expected = p.mul(
        &prefactor,
        &p.pow(&p.complex(2, -1), &p.neg(&p.rational(&epsilon))),
    );
    for value in [&first[0], &second[0]] {
        assert!(p.close(value, &expected, 20), "{value} != {expected}");
    }
    assert!(p.close(&first[0], &second[0], 20));
}
