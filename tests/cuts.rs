use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::prelude::*;
use symbolica_amflow::cache::CachedBackend;
use symbolica_amflow::cuts::*;
use symbolica_amflow::reduction::Reduction;
use symbolica_amflow::*;

fn family(reverse: bool) -> CutFamily {
    let gram = vec![vec![Atom::num(25)]];
    let ordinary = IntegralFamily {
        name: "cut_bubble_native".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[1], &[0], Atom::num(1), &gram).unwrap(),
            Propagator::quadratic(&[-1], &[1], Atom::num(4), &gram).unwrap(),
        ],
        physical_propagators: 2,
        epsilon: symbol!("cut_backend_test::eps"),
        dimension: 4,
    };
    let sign = if reverse { -1 } else { 1 };
    CutFamily::new(
        ordinary,
        CutMetadata::new(
            [
                CutLine {
                    propagator: 0,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![sign.into()],
                            external: vec![0.into()],
                        },
                    },
                },
                CutLine {
                    propagator: 1,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![(-sign).into()],
                            external: vec![sign.into()],
                        },
                    },
                },
            ],
            vec![LoopPrescription::Insensitive],
        )
        .unwrap(),
    )
    .unwrap()
}

fn targets() -> Vec<Integral> {
    vec![
        Integral(vec![1, 1]),
        Integral(vec![2, 1]),
        Integral(vec![1, 2]),
        Integral(vec![0, 1]),
        Integral(vec![1, -1]),
    ]
}
fn assert_reduction(reduction: &Reduction, dimension_minus_three: Atom) {
    assert_eq!(reduction.residuals, vec![Integral(vec![1, 1])]);
    for (target, coefficient) in [
        (Integral(vec![2, 1]), Atom::num((-7, 96))),
        (Integral(vec![1, 2]), Atom::num((-11, 192))),
    ] {
        let expansion = reduction.expand(&target).unwrap();
        assert_eq!(expansion.len(), 1);
        let expected = coefficient * &dimension_minus_three;
        assert!(
            (&expansion[&Integral(vec![1, 1])] - expected)
                .together()
                .cancel()
                .is_zero()
        );
    }
    for target in [Integral(vec![0, 1]), Integral(vec![1, -1])] {
        assert!(reduction.expand(&target).unwrap().is_empty());
    }
}

#[test]
fn cut_backend_reduces_raised_and_pinched_two_body_integrals_exactly() {
    let family = family(false);
    let context = RunContext::default();
    for factorized in [true, false] {
        let backend = RustRedBackend {
            max_depth: 3,
            max_targets: 256,
            factorized,
            ..Default::default()
        };
        let symbolic = backend.reduce_cut(&family, &targets(), &context).unwrap();
        assert_reduction(
            &symbolic,
            Atom::num(1) - Atom::num(2) * Atom::var(family.family().epsilon),
        );
        let sampled = backend
            .reduce_cut_at_epsilon(&family, &targets(), &Rational::from((1, 13)), &context)
            .unwrap();
        assert_reduction(&sampled, Atom::num((11, 13)));
    }
}

#[test]
fn cut_cache_separates_orientation_and_epsilon_and_supports_restart() {
    struct Counted {
        backend: RustRedBackend,
        calls: AtomicUsize,
    }
    impl ReductionBackend for Counted {
        fn identity(&self) -> String {
            self.backend.identity()
        }
        fn reduce(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<Reduction> {
            self.backend.reduce(family, targets, context)
        }
        fn reduce_cut(
            &self,
            family: &CutFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<Reduction> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            self.backend.reduce_cut(family, targets, context)
        }
        fn reduce_cut_at_epsilon(
            &self,
            family: &CutFamily,
            targets: &[Integral],
            epsilon: &Rational,
            context: &RunContext,
        ) -> Result<Reduction> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            self.backend
                .reduce_cut_at_epsilon(family, targets, epsilon, context)
        }
    }
    let directory =
        std::env::temp_dir().join(format!("amflow-cut-native-cache-{}", std::process::id()));
    let backend = CachedBackend {
        backend: Counted {
            backend: RustRedBackend {
                max_depth: 3,
                max_targets: 256,
                ..Default::default()
            },
            calls: AtomicUsize::new(0),
        },
        directory: directory.clone(),
    };
    let context = RunContext::default();
    let family = family(false);
    let targets = targets();
    let result = backend.reduce_cut(&family, &targets, &context).unwrap();
    let restarted = backend.reduce_cut(&family, &targets, &context).unwrap();
    assert_eq!(result.rules, restarted.rules);
    assert_eq!(backend.backend.calls.load(Ordering::Relaxed), 1);
    backend
        .reduce_cut(&self::family(true), &targets, &context)
        .unwrap();
    assert_eq!(backend.backend.calls.load(Ordering::Relaxed), 2);
    let sampled = backend
        .reduce_cut_at_epsilon(&family, &targets, &Rational::from((1, 13)), &context)
        .unwrap();
    assert_reduction(&sampled, Atom::num((11, 13)));
    backend
        .reduce_cut_at_epsilon(&family, &targets, &Rational::from((1, 13)), &context)
        .unwrap();
    assert_eq!(backend.backend.calls.load(Ordering::Relaxed), 3);
    assert_eq!(std::fs::read_dir(&directory).unwrap().count(), 3);
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn ordinary_backend_default_does_not_discard_the_cut_measure() {
    struct OrdinaryOnly;
    impl ReductionBackend for OrdinaryOnly {
        fn identity(&self) -> String {
            "ordinary-only-test".into()
        }
        fn reduce(
            &self,
            _: &IntegralFamily,
            targets: &[Integral],
            _: &RunContext,
        ) -> Result<Reduction> {
            Ok(Reduction {
                residuals: targets.to_vec(),
                rules: BTreeMap::new(),
                nonzero_conditions: vec![],
            })
        }
    }
    assert!(matches!(
        OrdinaryOnly.reduce_cut(&family(false), &targets(), &RunContext::default()),
        Err(Error::Unsupported(_))
    ));
}
