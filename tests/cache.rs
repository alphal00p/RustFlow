use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::prelude::*;
use symbolica_amflow::reduction::Reduction;
use symbolica_amflow::*;

struct CountBackend(AtomicUsize);
impl ReductionBackend for CountBackend {
    fn identity(&self) -> String {
        "cache-test-v1".into()
    }
    fn reduce(
        &self,
        _: &IntegralFamily,
        targets: &[Integral],
        _: &RunContext,
    ) -> Result<Reduction> {
        self.0.fetch_add(1, Ordering::Relaxed);
        Ok(Reduction {
            residuals: targets.to_vec(),
            nonzero_conditions: vec![parse!("m2")],
            ..Default::default()
        })
    }
}

#[test]
fn restart_invalidation_and_corruption() {
    let directory = std::env::temp_dir().join(format!(
        "symbolica-amflow-cache-test-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&directory);
    let backend = cache::CachedBackend {
        backend: CountBackend(AtomicUsize::new(0)),
        directory: directory.clone(),
    };
    let mut family = IntegralFamily {
        name: "tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: parse!("-m2"),
            scalar_products: vec![Atom::num(1)],
        }],
        physical_propagators: 1,
        epsilon: symbol!("eps"),
        dimension: 4,
    };
    let target = [Integral(vec![1])];
    let context = RunContext::default();
    backend.reduce(&family, &target, &context).unwrap();
    let restored = backend.reduce(&family, &target, &context).unwrap();
    assert_eq!(restored.nonzero_conditions, vec![parse!("m2")]);
    assert_eq!(backend.backend.0.load(Ordering::Relaxed), 1);
    family.dimension = 6;
    backend.reduce(&family, &target, &context).unwrap();
    assert_eq!(backend.backend.0.load(Ordering::Relaxed), 2);
    for file in std::fs::read_dir(&directory).unwrap() {
        std::fs::write(file.unwrap().path(), b"broken").unwrap();
    }
    assert!(matches!(
        backend.reduce(&family, &target, &context),
        Err(Error::Cache(_))
    ));
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn published_input_is_authenticated_and_preserves_targets() {
    let (family, targets) = benchmarks::paper_two_loop().unwrap();
    family.validate().unwrap();
    assert_eq!(targets.len(), 4);
    assert_eq!(targets[0].0, vec![1, 1, 1, 1, 1, 1, 1, -3, 0]);
    assert_eq!(targets[3].0, vec![1, 1, 1, 1, 1, 1, 1, 0, -3]);
    for target in &targets {
        family.validate_integral(target).unwrap();
    }
}

#[test]
fn prepared_system_restarts_without_calling_reducer() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct Counting {
        calls: AtomicUsize,
    }
    impl ReductionBackend for Counting {
        fn identity(&self) -> String {
            "prepared-cache-counting-v1".into()
        }
        fn reduce(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            RustRedBackend::default().reduce(family, targets, context)
        }
    }
    let directory = std::env::temp_dir().join(format!(
        "symbolica-amflow-system-cache-{}",
        std::process::id()
    ));
    let options = FlowOptions {
        cache_directory: Some(directory.clone()),
        ..Default::default()
    };
    let family = IntegralFamily {
        name: "restart_tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: Atom::num(-1),
            scalar_products: vec![Atom::num(1)],
        }],
        physical_propagators: 1,
        epsilon: symbol!("cache_eps"),
        dimension: 4,
    };
    let backend = Counting {
        calls: AtomicUsize::new(0),
    };
    let context = RunContext::default();
    let targets = [Integral(vec![1])];
    let first = PreparedFlow::new(
        &family,
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let calls = backend.calls.load(Ordering::Relaxed);
    let second = PreparedFlow::new(
        &family,
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    assert_eq!(backend.calls.load(Ordering::Relaxed), calls);
    assert_eq!(first.reduced.basis, second.reduced.basis);
    assert_eq!(first.reduced.matrix, second.reduced.matrix);
    assert_eq!(first.reduced.targets, second.reduced.targets);
    let changed = FlowOptions {
        skip_reduction: true,
        ..options
    };
    PreparedFlow::new(
        &family,
        &targets,
        &KinematicPoint::default(),
        &backend,
        &changed,
        &context,
    )
    .unwrap();
    assert!(backend.calls.load(Ordering::Relaxed) > calls);
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn automatic_mass_placement_uses_intrinsic_mass() {
    let (family, _) = benchmarks::paper_two_loop().unwrap();
    for i in 0..7 {
        assert_eq!(
            family.mass_squared(i).unwrap(),
            Atom::num(i64::from(i == 4))
        );
    }
    let (_, mask) = family
        .deform(symbol!("placement_eta"), &MassMode::Auto)
        .unwrap();
    assert_eq!(
        mask,
        vec![false, false, false, false, true, false, false, false, false]
    );
}

#[test]
fn native_reduction_restarts_after_increasing_the_target_budget() {
    let directory = std::env::temp_dir().join(format!(
        "symbolica-amflow-native-stage-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&directory);
    let family = IntegralFamily {
        name: "checkpoint_tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator::quadratic(&[1], &[], parse!("stage_m2"), &[]).unwrap()],
        physical_propagators: 1,
        epsilon: symbol!("stage_eps"),
        dimension: 4,
    };
    let backend = RustRedBackend {
        factorized: true,
        max_targets: 1,
        checkpoints: Some(directory.clone()),
        ..Default::default()
    };
    let target = [Integral(vec![2])];
    assert!(matches!(
        backend.reduce(&family, &target, &RunContext::default()),
        Err(Error::Limit(_))
    ));
    let calls = std::sync::Arc::new(AtomicUsize::new(0));
    let counter = calls.clone();
    let context = RunContext {
        progress: Some(std::sync::Arc::new(move |event| {
            if matches!(event, Progress::SectorReduction { .. }) {
                counter.fetch_add(1, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let resumed = RustRedBackend {
        max_targets: 16,
        ..backend
    }
    .reduce(&family, &target, &context)
    .unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), 1);
    let terms = resumed.expand(&target[0]).unwrap();
    assert!(
        (&terms[&Integral(vec![1])] - parse!("(1-stage_eps)/stage_m2"))
            .together()
            .cancel()
            .is_zero()
    );
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn mass_placement_modes_follow_intrinsic_masses_and_loop_topology() {
    let (mut family, _) = benchmarks::paper_two_loop().unwrap();
    let selected = |family: &IntegralFamily, mode| {
        family
            .deform(symbol!("placement_eta"), &mode)
            .unwrap()
            .1
            .into_iter()
            .enumerate()
            .filter_map(|(i, shifted)| shifted.then_some(i))
            .collect::<Vec<_>>()
    };
    assert_eq!(selected(&family, MassMode::Mass), vec![4]);
    assert_eq!(selected(&family, MassMode::Propagator), vec![0]);
    assert_eq!(selected(&family, MassMode::Branch), vec![6]);
    assert_eq!(selected(&family, MassMode::Loop), vec![0, 1, 2, 6]);
    family.propagators[0].constant -= Atom::num(1);
    assert_eq!(selected(&family, MassMode::Mass), vec![0, 4]);
    family.propagators[1].constant -= Atom::num(2);
    assert_eq!(selected(&family, MassMode::Mass), vec![1]);
    family.propagators[0].constant += Atom::num(1);
    family.propagators[1].constant += Atom::num(2);
    family.propagators[4].constant += Atom::num(1);
    assert!(matches!(
        family.deform(symbol!("placement_eta"), &MassMode::Mass),
        Err(Error::Unsupported(_))
    ));
}
