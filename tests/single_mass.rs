//! Proposed public-route regressions. Run after integrating recursive.patch.
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
        }
    }
    fn record(&self, family: &IntegralFamily) {
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
