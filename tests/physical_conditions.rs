use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{kinematics::KinematicSystem, transport_cache::*, *};

fn symbols() -> (Symbol, Symbol) {
    symbol!(
        "physical_conditions_test::s",
        "physical_conditions_test::epsilon"
    )
}
fn flow(conditions: &[Atom]) -> RustFlow {
    let (s, epsilon) = symbols();
    RustFlow::with_conditions(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]])]),
        },
        &[parse!("physical_conditions_test::I(1)")],
        &Atom::one(),
        Prescription::PlusI0,
        "single regular sheet",
        conditions,
    )
    .unwrap()
}
fn point(value: Atom) -> CachedPoint {
    CachedPoint::Exact(BTreeMap::from([(symbols().0, value)]))
}
fn entry(flow: &RustFlow, value: Atom) -> CachedBoundary {
    let p = Precision::decimal(90).unwrap();
    CachedBoundary {
        identity: flow.identity().clone(),
        point: point(value),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 0).unwrap(),
        coefficients: vec![vec![p.i(1)]],
        accuracy: BoundaryAccuracy::supplied(
            80,
            p.bits,
            vec![vec![p.real(0)]],
            "exact constant solution one",
        )
        .unwrap(),
    }
}
#[test]
fn generic_epsilon_guards_do_not_substitute_epsilon_zero() {
    let (s, epsilon) = symbols();
    let s = Atom::var(s);
    let epsilon = Atom::var(epsilon);
    let p = Precision::decimal(60).unwrap();
    for guard in [epsilon.clone(), &s - &epsilon] {
        let engine = flow(&[guard]);
        let mut cache = RustFlowCache::default();
        cache.insert(entry(&engine, Atom::new())).unwrap();
        assert!(
            engine
                .identity()
                .conditions_admit_straight_path(&point(Atom::num(-1)), &point(Atom::num(1)), p, 20)
                .unwrap()
        );
    }
    let engine = flow(&[&epsilon * &s]);
    assert!(
        RustFlowCache::default()
            .insert(entry(&engine, Atom::new()))
            .is_err()
    );
    assert!(
        !engine
            .identity()
            .conditions_admit_straight_path(&point(Atom::num(-1)), &point(Atom::num(1)), p, 20)
            .unwrap()
    );
}
#[test]
fn original_guard_denominators_survive_later_rational_cancellation() {
    let s = Atom::var(symbols().0);
    let raw = ((s.clone() + Atom::one()).pow(2) - Atom::one()) / (&s * (&s + Atom::num(2)));
    assert!(raw.together().cancel().is_one());
    let engine = flow(&[raw]);
    assert!(!engine.identity().nonzero_conditions().is_empty());
    assert!(
        RustFlowCache::default()
            .insert(entry(&engine, Atom::new()))
            .is_err()
    );
    assert!(
        !engine
            .identity()
            .conditions_admit_straight_path(
                &point(Atom::num(-1)),
                &point(Atom::num(1)),
                Precision::decimal(60).unwrap(),
                20
            )
            .unwrap()
    );
}
#[test]
fn conditions_are_part_of_persistent_identity_and_old_schema_is_rejected() {
    let s = Atom::var(symbols().0);
    let unguarded = flow(&[]);
    let guarded = flow(&[s]);
    assert_ne!(unguarded.identity().key(), guarded.identity().key());
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&guarded, Atom::num(1))).unwrap();
    let directory =
        std::env::temp_dir().join(format!("amflow-condition-cache-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    cache.save(&directory).unwrap();
    let restored = RustFlowCache::load(&directory).unwrap();
    assert_eq!(
        restored.entries()[0].identity.key(),
        guarded.identity().key()
    );
    assert_eq!(
        restored.entries()[0].identity.nonzero_conditions(),
        guarded.identity().nonzero_conditions()
    );
    let path = directory.join("physical-boundaries.bin");
    let mut bytes = std::fs::read(&path).unwrap();
    let version = b"AMFLOW-BOUNDARIES\0".len();
    assert_eq!(bytes[version], 2);
    bytes[version] = 1;
    std::fs::write(path, bytes).unwrap();
    assert!(matches!(
        RustFlowCache::load(&directory),
        Err(Error::Cache(_))
    ));
    std::fs::remove_dir_all(directory).unwrap();
}
#[test]
fn excluded_exact_target_is_rejected_before_a_cache_hit_or_mutation() {
    let unguarded = flow(&[]);
    let guarded = flow(&[Atom::var(symbols().0)]);
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&unguarded, Atom::new())).unwrap();
    cache.insert(entry(&guarded, Atom::num(1))).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let result = guarded.evaluate_to(
        &mut cache,
        &BTreeMap::from([(symbols().0, Atom::new())]),
        EpsilonRange::new(0, 0).unwrap(),
        &FlowOptions::default(),
        &RunContext::default(),
        &policy,
    );
    assert!(matches!(result, Err(Error::InvalidInput(_))));
    assert_eq!(cache.len(), 2);
}
#[test]
fn unsafe_nearest_source_is_skipped_in_favor_of_a_safe_farther_boundary() {
    let engine = flow(&[Atom::var(symbols().0)]);
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&engine, Atom::num(-1))).unwrap();
    cache.insert(entry(&engine, Atom::num(4))).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let result = engine
        .evaluate_to(
            &mut cache,
            &BTreeMap::from([(symbols().0, Atom::num(1))]),
            EpsilonRange::new(0, 0).unwrap(),
            &FlowOptions::default(),
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert_eq!(
        result
            .starting_point
            .rounded_coordinates_as_exact()
            .unwrap()[&symbols().0],
        Atom::num(4)
    );
    assert_eq!(
        result.boundary.coefficients[0][0],
        Precision::decimal(80).unwrap().i(1)
    );
}
#[test]
fn complex_guard_coefficients_use_common_real_path_zeros() {
    let s = Atom::var(symbols().0);
    let p = Precision::decimal(60).unwrap();
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let engine = flow(&[&s + &imaginary]);
    assert!(
        engine
            .identity()
            .conditions_admit_straight_path(&point(Atom::num(-2)), &point(Atom::num(2)), p, 20)
            .unwrap()
    );
    let engine = flow(&[s]);
    assert!(
        !engine
            .identity()
            .conditions_admit_straight_path(&point(-imaginary.clone()), &point(imaginary), p, 20)
            .unwrap()
    );
}

#[test]
fn stationary_coordinates_do_not_hide_original_matrix_poles() {
    let (s, epsilon) = symbols();
    let t = symbol!("physical_conditions_test::t");
    // This scalar connection is flat: d(s/t). The restricted solution at s=0
    // is constant, but crossing t=0 leaves the domain of the supplied PDE.
    let engine = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![Atom::one() / Atom::var(t)]]),
                (t, vec![vec![-Atom::var(s) / Atom::var(t).pow(2)]]),
            ]),
        },
        &[parse!("physical_conditions_test::flat")],
        &Atom::one(),
        Prescription::PlusI0,
        "regular domain",
    )
    .unwrap();
    let point =
        |t_value| CachedPoint::Exact(BTreeMap::from([(s, Atom::new()), (t, Atom::num(t_value))]));
    assert!(
        !engine
            .identity()
            .conditions_admit_straight_path(
                &point(-1),
                &point(1),
                Precision::decimal(60).unwrap(),
                20
            )
            .unwrap()
    );
}

#[test]
fn original_matrix_epsilon_denominators_require_a_regular_laurent_chart() {
    let (s, epsilon) = symbols();
    let engine = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(
                s,
                vec![vec![Atom::one() / (Atom::var(s) - Atom::var(epsilon))]],
            )]),
        },
        &[parse!("physical_conditions_test::epsilon_denominator")],
        &Atom::one(),
        Prescription::PlusI0,
        "regular epsilon chart",
    )
    .unwrap();
    assert!(
        RustFlowCache::default()
            .insert(entry(&engine, Atom::new()))
            .is_err()
    );
    assert!(
        !engine
            .identity()
            .conditions_admit_straight_path(
                &point(Atom::num(-1)),
                &point(Atom::num(1)),
                Precision::decimal(60).unwrap(),
                20
            )
            .unwrap()
    );
    // The corresponding generic reduction guard, independent of a matrix
    // denominator, still admits this epsilon-dependent factor at s=0.
    let generic = flow(&[Atom::var(s) - Atom::var(epsilon)]);
    RustFlowCache::default()
        .insert(entry(&generic, Atom::new()))
        .unwrap();
}

#[test]
fn a_cheaper_original_matrix_pole_route_does_not_hide_a_safe_boundary() {
    let (s, epsilon) = symbols();
    let engine = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / Atom::var(s)]])]),
        },
        &[parse!("physical_conditions_test::epsilon_log")],
        &Atom::one(),
        Prescription::PlusI0,
        "regular epsilon chart",
    )
    .unwrap();
    // At the retained epsilon^0 order the solution is exactly one. Both
    // boundary values are exact; only the route from -1 crosses a matrix pole.
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&engine, Atom::num(-1))).unwrap();
    cache.insert(entry(&engine, Atom::num(4))).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let result = engine
        .evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(1))]),
            EpsilonRange::new(0, 0).unwrap(),
            &FlowOptions::default(),
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert_eq!(
        result
            .starting_point
            .rounded_coordinates_as_exact()
            .unwrap()[&s],
        Atom::num(4)
    );
}
