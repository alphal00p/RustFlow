use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::physical_transport::PhysicalTransport;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

fn problem() -> (PhysicalTransport, RustFlowCache) {
    let system = KinematicSystem {
        epsilon: symbol!("eps"),
        derivatives: BTreeMap::from([(symbol!("s"), vec![vec![parse!("eps/s")]])]),
    };
    let engine = PhysicalTransport::new(
        system,
        &[parse!("I1")],
        &Atom::num(1),
        Prescription::PlusI0,
        "positive real s, principal logarithm",
    )
    .unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut coefficients = vec![vec![p.zero()]; 5];
    coefficients[0][0] = p.i(1);
    let mut cache = RustFlowCache::default();
    cache
        .insert(CachedBoundary {
            identity: engine.identity().clone(),
            point: CachedPoint::Exact(BTreeMap::from([(symbol!("s"), Atom::num(1))])),
            kind: PointKind::Physical,
            range: EpsilonRange::new(0, 4).unwrap(),
            coefficients,
            accuracy: BoundaryAccuracy::supplied(
                80,
                p.bits,
                vec![vec![p.real(0)]; 5],
                "exact analytic boundary Y(s=1)=1 for Y=s^eps",
            )
            .unwrap(),
        })
        .unwrap();
    (engine, cache)
}

#[test]
fn repeated_physical_queries_reuse_verified_intermediates_and_persistent_boundaries() {
    let (engine, mut cache) = problem();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, target: &CachedPoint| {
            let p = Precision::decimal(40)?;
            Ok(target.evaluate(p)?[&symbol!("s")].re > p.real(0))
        },
    };
    let range = EpsilonRange::new(0, 4).unwrap();
    let options = FlowOptions::default();
    let first = engine
        .evaluate_to(
            &mut cache,
            &BTreeMap::from([(symbol!("s"), Atom::num(16))]),
            range,
            &options,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    let cold_steps = first.transport.as_ref().unwrap().diagnostics.steps;
    assert!(first.inserted_points > 2);
    assert!(first.boundary.accuracy.verified_digits() > options.digits);
    assert!(
        cache
            .entries()
            .iter()
            .any(|entry| matches!(entry.point, CachedPoint::Derived { .. }))
    );
    let directory =
        std::env::temp_dir().join(format!("amflow-physical-restart-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    cache.save(&directory).unwrap();
    let (resumed, _) = problem();
    let mut resumed_cache = RustFlowCache::load(&directory).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();
    let destination = BTreeMap::from([(symbol!("s"), Atom::num(17))]);
    let nearby = resumed
        .evaluate_to(
            &mut resumed_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    let start = nearby
        .starting_point
        .rounded_coordinates_as_exact()
        .unwrap();
    assert!(
        (&start[&symbol!("s")] - Atom::num(16))
            .together()
            .cancel()
            .is_zero()
    );
    assert!(nearby.transport.as_ref().unwrap().diagnostics.steps < cold_steps);
    let p = Precision {
        bits: nearby.boundary.accuracy.working_bits(),
    };
    let mut factorial = 1;
    for (k, row) in nearby.boundary.coefficients.iter().enumerate() {
        if k > 0 {
            factorial *= k as i64;
        }
        assert!(p.close(
            &row[0],
            &p.scale(&p.powi(&p.log(&p.i(17)), k as i64), 1, factorial),
            35
        ));
    }
    let exact_hit = resumed
        .evaluate_to(
            &mut resumed_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert!(exact_hit.transport.is_none());
    assert_eq!(exact_hit.inserted_points, 0);
    assert_eq!(
        exact_hit.boundary.coefficients,
        nearby.boundary.coefficients
    );
    // An actual intermediate value supplies another future destination; it is
    // not replaced by a nearby grid value or an interpolation of endpoints.
    let entry = resumed_cache
        .entries()
        .iter()
        .find(|entry| {
            let s = entry.point.evaluate(p).unwrap()[&symbol!("s")].clone();
            s.re > p.real(2) && s.re < p.real(8)
        })
        .unwrap();
    let intermediate = entry.point.rounded_coordinates_as_exact().unwrap();
    let hit = resumed
        .evaluate_to(
            &mut resumed_cache,
            &intermediate,
            range,
            &options,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert!(hit.transport.is_none());
}

#[test]
fn too_weak_cached_input_is_not_promoted_by_repeating_transport() {
    let (engine, mut cache) = problem();
    let mut weak = cache.entries()[0].clone();
    let p = Precision {
        bits: weak.accuracy.working_bits(),
    };
    weak.accuracy = BoundaryAccuracy::supplied(
        20,
        p.bits,
        vec![vec![p.real(0)]; 5],
        "only twenty digits established despite many working bits",
    )
    .unwrap();
    cache = BoundaryCache::default();
    cache.insert(weak).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let result = engine.evaluate_to(
        &mut cache,
        &BTreeMap::from([(symbol!("s"), Atom::num(2))]),
        EpsilonRange::new(0, 4).unwrap(),
        &FlowOptions::default(),
        &RunContext::default(),
        &policy,
    );
    assert!(matches!(result, Err(Error::Accuracy(_))));
    assert_eq!(cache.len(), 1);
}
