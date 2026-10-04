use std::collections::{BTreeMap, BTreeSet};
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};
use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn verified_routing_runs_through_native_reduction_and_isolates_restarts() {
    let family = IntegralFamily {
        name: "native_symmetry_equal_mass_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-2)]],
        propagators: vec![
            Propagator {
                constant: Atom::num(-1),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            Propagator {
                constant: Atom::num(-3),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("native_symmetry_bubble_eps"),
        dimension: 4,
    };
    let eta = symbol!("native_symmetry_bubble_eta");
    let (deformed, _) = family.deform(eta, &MassMode::All).unwrap();
    let epsilon = Rational::from((1, 10));
    let targets = [Integral(vec![2, 0]), Integral(vec![0, 2])];
    let directory = std::env::temp_dir().join(format!(
        "amflow-native-symmetry-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let baseline = RustRedBackend {
        checkpoints: Some(directory.clone()),
        ..Default::default()
    };
    let symmetric = RustRedBackend {
        symmetry_rules: true,
        ..baseline.clone()
    };
    assert_eq!(
        symmetric.identity(),
        format!("{}:symmetry-v1", baseline.identity())
    );
    let expected = baseline
        .reduce_at_epsilon(&deformed, &targets, &epsilon, &RunContext::default())
        .unwrap();
    let cached_names = || {
        std::fs::read_dir(&directory)
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .filter(|name| {
                !name
                    .to_string_lossy()
                    .starts_with("native-stage-family-v1-")
            })
            .collect::<BTreeSet<_>>()
    };
    let before = cached_names();
    assert_eq!(before.len(), 1);
    let applied = Arc::new(AtomicUsize::new(0));
    let searched = Arc::new(AtomicUsize::new(0));
    let observed_applied = Arc::clone(&applied);
    let observed_searched = Arc::clone(&searched);
    let context = RunContext {
        progress: Some(Arc::new(move |event| match event {
            Progress::SymmetryReduction { applied, .. } => {
                observed_applied.fetch_add(applied, Ordering::Relaxed);
            }
            Progress::SectorReduction { integrals, .. } => {
                observed_searched.fetch_add(integrals, Ordering::Relaxed);
            }
            _ => {}
        })),
        ..Default::default()
    };
    let actual = symmetric
        .reduce_at_epsilon(&deformed, &targets, &epsilon, &context)
        .unwrap();
    assert!(applied.load(Ordering::Relaxed) > 0, "no routing rule used");
    assert!(
        searched.load(Ordering::Relaxed) > 0,
        "no concrete fallback used"
    );
    assert!(
        actual
            .nonzero_conditions
            .iter()
            .any(|c| !c.derivative(eta).is_zero())
    );
    for target in &targets {
        // The baseline keeps both translated tadpole masters. For equal
        // masses their exact shift l -> l-p identifies I(1,0) with I(0,1).
        let mut common_basis = BTreeMap::<Integral, Atom>::new();
        for (mut integral, coefficient) in expected.expand(target).unwrap() {
            if integral == Integral(vec![1, 0]) {
                integral = Integral(vec![0, 1]);
            }
            let entry = common_basis.entry(integral).or_default();
            *entry = (&*entry + coefficient).together().cancel();
        }
        let actual = actual.expand(target).unwrap();
        assert_eq!(
            actual.keys().collect::<Vec<_>>(),
            common_basis.keys().collect::<Vec<_>>()
        );
        for (integral, coefficient) in common_basis {
            assert!(
                (&actual[&integral] - coefficient)
                    .together()
                    .cancel()
                    .is_zero()
            );
        }
    }
    let after = cached_names();
    assert_eq!(
        after.len(),
        2,
        "symmetry and baseline shared a native checkpoint"
    );
    assert!(before.is_subset(&after));
    let applied_before_restart = applied.load(Ordering::Relaxed);
    let searched_before_restart = searched.load(Ordering::Relaxed);
    let resumed = symmetric
        .reduce_at_epsilon(&deformed, &targets, &epsilon, &context)
        .unwrap();
    assert_eq!(applied.load(Ordering::Relaxed), applied_before_restart);
    assert_eq!(searched.load(Ordering::Relaxed), searched_before_restart);
    for target in &targets {
        assert_eq!(
            resumed.expand(target).unwrap(),
            actual.expand(target).unwrap()
        );
    }
    std::fs::remove_dir_all(directory).unwrap();
}
