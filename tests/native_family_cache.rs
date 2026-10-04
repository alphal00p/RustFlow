use std::sync::{
    Arc,
    atomic::{AtomicBool, AtomicUsize, Ordering},
};
use symbolica::prelude::*;
use symbolica_amflow::*;

fn family() -> IntegralFamily {
    IntegralFamily {
        name: "family_bank_tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator::quadratic(&[1], &[], parse!("family_bank_m2"), &[]).unwrap()],
        physical_propagators: 1,
        epsilon: symbol!("family_bank_eps"),
        dimension: 4,
    }
}

fn counted_context() -> (RunContext, Arc<AtomicUsize>) {
    let searches = Arc::new(AtomicUsize::new(0));
    let count = Arc::clone(&searches);
    (
        RunContext {
            progress: Some(Arc::new(move |event| {
                if let Progress::SectorReduction { integrals, .. } = event {
                    count.fetch_add(integrals, Ordering::Relaxed);
                }
            })),
            ..Default::default()
        },
        searches,
    )
}

#[test]
fn different_targets_reuse_equations_but_research_leaves_and_isolate_parameters() {
    let directory = std::env::temp_dir().join(format!("amflow-family-bank-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    let backend = RustRedBackend {
        checkpoints: Some(directory.clone()),
        ..Default::default()
    };
    let fresh = RustRedBackend::default();
    let family = family();
    let epsilon = Rational::from((1, 10));
    let first = [Integral(vec![2]), Integral(vec![3])];
    backend
        .reduce_at_epsilon(&family, &first, &epsilon, &RunContext::default())
        .unwrap();
    let target = [Integral(vec![2])];
    let expected = fresh
        .reduce_at_epsilon(&family, &target, &epsilon, &RunContext::default())
        .unwrap();
    let (context, searches) = counted_context();
    let actual = backend
        .reduce_at_epsilon(&family, &target, &epsilon, &context)
        .unwrap();
    assert_eq!(
        actual.expand(&target[0]).unwrap(),
        expected.expand(&target[0]).unwrap()
    );
    assert_eq!(
        searches.load(Ordering::Relaxed),
        1,
        "the old residual must be searched again; the I(2) equation is reusable"
    );
    assert!(
        actual
            .nonzero_conditions
            .iter()
            .any(|p| !p.derivative(symbol!("family_bank_m2")).is_zero())
    );
    // The completed target-specific checkpoint is preferred over the shared
    // bank and therefore retains its own valid leaf search history.
    backend
        .reduce_at_epsilon(&family, &target, &epsilon, &context)
        .unwrap();
    assert_eq!(searches.load(Ordering::Relaxed), 1);

    let other_epsilon = Rational::from((1, 9));
    let mut other_family = family.clone();
    other_family.propagators[0].constant -= Atom::num(1);
    for (input, epsilon) in [(&family, &other_epsilon), (&other_family, &epsilon)] {
        let (context, searches) = counted_context();
        let actual = backend
            .reduce_at_epsilon(input, &target, epsilon, &context)
            .unwrap();
        let expected = fresh
            .reduce_at_epsilon(input, &target, epsilon, &RunContext::default())
            .unwrap();
        assert_eq!(
            actual.expand(&target[0]).unwrap(),
            expected.expand(&target[0]).unwrap()
        );
        assert!(
            searches.load(Ordering::Relaxed) >= 2,
            "incompatible family or epsilon reused old equations"
        );
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn raw_checkpoint_exists_before_exact_substitution_starts() {
    let directory =
        std::env::temp_dir().join(format!("amflow-before-replay-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    let backend = RustRedBackend {
        checkpoints: Some(directory.clone()),
        checkpoint_interval: std::time::Duration::from_secs(3600),
        ..Default::default()
    };
    let saved_before_replay = Arc::new(AtomicBool::new(false));
    let observed = Arc::clone(&saved_before_replay);
    let observed_directory = directory.clone();
    let cancellation = CancellationToken::default();
    let token = cancellation.clone();
    let context = RunContext {
        cancellation,
        progress: Some(Arc::new(move |event| {
            if let Progress::SubstitutionPlan { .. } = event {
                observed.store(
                    std::fs::read_dir(&observed_directory)
                        .ok()
                        .is_some_and(|entries| {
                            entries.flatten().any(|entry| {
                                entry
                                    .file_name()
                                    .to_string_lossy()
                                    .starts_with("native-stage-")
                            })
                        }),
                    Ordering::Relaxed,
                );
                token.cancel();
            }
        })),
    };
    let result = backend.reduce(&family(), &[Integral(vec![2])], &context);
    assert!(matches!(result, Err(Error::Cancelled)));
    assert!(saved_before_replay.load(Ordering::Relaxed));
    std::fs::remove_dir_all(directory).unwrap();
}
