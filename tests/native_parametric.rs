use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};
use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn guarded_rays_preserve_massive_bubble_reduction_and_flow() {
    let family = IntegralFamily {
        name: "guarded_ray_bubble".into(),
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
        epsilon: symbol!("guarded_ray_bubble_eps"),
        dimension: 4,
    };
    let baseline = RustRedBackend::default();
    let parametric = RustRedBackend {
        parametric_rules: true,
        ..Default::default()
    };
    let applied = Arc::new(AtomicUsize::new(0));
    let observed = applied.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if let Progress::ParametricReduction { applied, .. } = event {
                observed.fetch_add(applied, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let epsilon = Rational::from((1, 10));
    let eta = symbol!("guarded_ray_bubble_eta");
    let (deformed, _) = family.deform(eta, &MassMode::All).unwrap();
    let targets = [
        Integral(vec![2, 1]),
        Integral(vec![1, 2]),
        Integral(vec![3, 1]),
    ];
    let expected = baseline
        .reduce_at_epsilon(&deformed, &targets, &epsilon, &context)
        .unwrap();
    let actual = parametric
        .reduce_at_epsilon(&deformed, &targets, &epsilon, &context)
        .unwrap();
    assert!(
        applied.load(Ordering::Relaxed) > 0,
        "test did not exercise a parametric specialization"
    );
    assert!(
        actual
            .nonzero_conditions
            .iter()
            .any(|c| !c.derivative(eta).is_zero()),
        "eta-dependent reduction guards were lost"
    );
    for target in &targets {
        let expected = expected.expand(target).unwrap();
        let actual = actual.expand(target).unwrap();
        assert_eq!(
            actual.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>()
        );
        for (integral, coefficient) in expected {
            assert!(
                (&actual[&integral] - coefficient)
                    .together()
                    .cancel()
                    .is_zero()
            );
        }
    }
    let options = FlowOptions::default();
    let target = Integral(vec![1, 1]);
    let ft = ft::evaluate_bubble(&family, &target, &epsilon, &options, &context).unwrap();
    let p = Precision::decimal(60).unwrap();
    let before_flow = applied.load(Ordering::Relaxed);
    for backend in [&baseline, &parametric] {
        let flow = PreparedFlow::new_at_epsilon(
            &family,
            std::slice::from_ref(&target),
            &KinematicPoint::default(),
            backend,
            &options,
            &context,
            &epsilon,
        )
        .unwrap();
        let values = flow
            .evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context)
            .unwrap();
        assert!(
            p.close(&values[0], &ft, 20),
            "flow {} disagrees with FT {ft}",
            values[0]
        );
    }
    assert!(
        applied.load(Ordering::Relaxed) > before_flow,
        "differential closure did not reuse guarded rules"
    );
}
