use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSolution, EpsilonSystem};
use symbolica_amflow::kinematics::{KinematicPath, KinematicSystem};
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{Error, FlowDiagnostics, FlowOptions, Precision, Prescription, RunContext};

fn fixture() -> (Symbol, Symbol, KinematicSystem, BoundaryIdentity) {
    let (s, epsilon) = symbol!("boundary_cache::s", "boundary_cache::eps");
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon)]])]),
    };
    let identity = BoundaryIdentity::new(
        &system,
        &[parse!("boundary_cache::I(1)")],
        &Atom::num(1),
        Prescription::PlusI0,
        "physical-sheet",
    )
    .unwrap();
    (s, epsilon, system, identity)
}

fn entry(identity: &BoundaryIdentity, s: Symbol, coordinate: Atom, digits: u32) -> CachedBoundary {
    let p = Precision::decimal(90).unwrap();
    CachedBoundary {
        identity: identity.clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, coordinate)])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 1).unwrap(),
        coefficients: vec![vec![p.i(1)], vec![p.i(2)]],
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.real(0)]; 2],
            "independent analytic boundary",
        )
        .unwrap(),
    }
}

#[test]
fn selection_filters_identity_accuracy_range_and_unsafe_routes() {
    let (s, _, system, identity) = fixture();
    let mut cache = BoundaryCache::default();
    cache.insert(entry(&identity, s, Atom::num(1), 12)).unwrap();
    cache.insert(entry(&identity, s, Atom::num(2), 35)).unwrap();
    cache.insert(entry(&identity, s, Atom::num(4), 40)).unwrap();
    for (basis, normalization, prescription, domain) in [
        (
            parse!("boundary_cache::J(1)"),
            Atom::num(1),
            Prescription::PlusI0,
            "physical-sheet",
        ),
        (
            parse!("boundary_cache::I(1)"),
            Atom::num(2),
            Prescription::PlusI0,
            "physical-sheet",
        ),
        (
            parse!("boundary_cache::I(1)"),
            Atom::num(1),
            Prescription::MinusI0,
            "physical-sheet",
        ),
        (
            parse!("boundary_cache::I(1)"),
            Atom::num(1),
            Prescription::PlusI0,
            "other-monodromy",
        ),
    ] {
        let incompatible =
            BoundaryIdentity::new(&system, &[basis], &normalization, prescription, domain).unwrap();
        assert_ne!(incompatible.key(), identity.key());
        cache
            .insert(entry(&incompatible, s, Atom::new(), 50))
            .unwrap();
    }
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::new())]));
    let query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 1).unwrap(), 20).unwrap();
    let distance = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let p = Precision::decimal(90).unwrap();
    let selected = cache.best(&query, &distance, p).unwrap().unwrap();
    assert_eq!(
        selected
            .boundary
            .point
            .rounded_coordinates_as_exact()
            .unwrap()[&s],
        Atom::num(2)
    );
    assert_eq!(selected.cost, p.real(4));
    let safe = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |boundary: &CachedBoundary, _: &CachedPoint| {
            Ok(boundary.point.rounded_coordinates_as_exact()?[&s] != Atom::num(2))
        },
    };
    assert_eq!(
        cache.best(&query, &safe, p).unwrap().unwrap().cost,
        p.real(16)
    );
    let insufficient =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 2).unwrap(), 20).unwrap();
    assert!(cache.best(&insufficient, &distance, p).unwrap().is_none());
    let omitted_pole =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(1, 1).unwrap(), 20).unwrap();
    assert!(cache.best(&omitted_pole, &distance, p).unwrap().is_none());
}

#[test]
fn scaled_distance_distinguishes_coordinates_beyond_binary64_resolution() {
    let (s, _, _, identity) = fixture();
    let huge = Atom::num(10).pow(100);
    let mut cache = BoundaryCache::default();
    cache
        .insert(entry(&identity, s, &huge + Atom::num(2), 40))
        .unwrap();
    cache
        .insert(entry(&identity, s, &huge + Atom::num(1), 40))
        .unwrap();
    let target = CachedPoint::Exact(BTreeMap::from([(s, huge)]));
    let query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 0).unwrap(), 20).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::from([(s, Atom::num(10))]),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    // The common 10^100 offset is cancelled exactly before MPFR evaluation.
    let p = Precision::decimal(40).unwrap();
    let best = cache.best(&query, &policy, p).unwrap().unwrap();
    assert!(p.close(
        &symbolica_amflow::ComplexFloat::new(best.cost, p.real(0)),
        &p.parse("0.01", "0").unwrap(),
        35
    ));
}

#[test]
fn rounded_cost_ties_preserve_nearest_sources_and_exact_hits() {
    let (s, _, _, identity) = fixture();
    let p = Precision::decimal(40).unwrap();
    let gap = Atom::num(10).pow(-80);
    let mut cache = BoundaryCache::default();
    // The slightly farther source has stronger evidence and was inserted first.
    // Both distances to 2 round to 1 at the query's working precision.
    cache
        .insert(entry(&identity, s, Atom::num(1) - &gap, 40))
        .unwrap();
    cache.insert(entry(&identity, s, Atom::num(1), 21)).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::num(2))]));
    let query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 0).unwrap(), 20).unwrap();
    let selected = cache.best(&query, &policy, p).unwrap().unwrap();
    assert_eq!(selected.cost, p.real(1));
    assert_eq!(
        selected.boundary.point.restart_coordinates().unwrap()[&s],
        Atom::num(1)
    );

    // A general cost policy keeps its own cost model while preferring a true
    // coordinate hit over an equally priced transport, independently of order.
    struct ConstantCost;
    impl TransportCost for ConstantCost {
        fn cost(
            &self,
            _: &CachedBoundary,
            _: &CachedPoint,
            p: Precision,
        ) -> symbolica_amflow::Result<Option<Float>> {
            Ok(Some(p.real(1)))
        }
    }
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::num(1))]));
    let query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 0).unwrap(), 20).unwrap();
    assert_eq!(
        cache
            .best(&query, &ConstantCost, p)
            .unwrap()
            .unwrap()
            .boundary
            .point
            .restart_coordinates()
            .unwrap()[&s],
        Atom::num(1)
    );

    // Genuine algebraic distances use consistent higher-precision comparisons.
    let root = Atom::num(2).pow(Atom::num((1, 2)));
    let mut algebraic = BoundaryCache::default();
    algebraic
        .insert(entry(&identity, s, &root - Atom::num(10).pow(-60), 40))
        .unwrap();
    algebraic
        .insert(entry(&identity, s, root.clone(), 21))
        .unwrap();
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::num(2))]));
    let query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 0).unwrap(), 20).unwrap();
    assert_eq!(
        algebraic
            .best(&query, &policy, p)
            .unwrap()
            .unwrap()
            .boundary
            .point
            .restart_coordinates()
            .unwrap()[&s],
        root
    );
}

#[test]
fn numerical_points_remain_distinct_and_detours_require_opt_in() {
    let (s, _, _, identity) = fixture();
    let p = Precision::decimal(90).unwrap();
    let mut exact = entry(&identity, s, Atom::num((1, 10)), 40);
    let numerical_point = CachedPoint::Numerical {
        coordinates: BTreeMap::from([(s, p.parse("0.1", "0").unwrap())]),
        working_bits: p.bits,
        provenance: "accepted regular physical checkpoint".into(),
    };
    let rationalized = numerical_point.rounded_coordinates_as_exact().unwrap();
    assert_ne!(rationalized[&s], Atom::num((1, 10)));
    assert_eq!(
        p.eval(&rationalized[&s], &Default::default()).unwrap(),
        numerical_point.evaluate(p).unwrap()[&s]
    );
    let mut numerical = exact.clone();
    numerical.point = numerical_point.clone();
    let mut cache = BoundaryCache::default();
    cache.insert(exact.clone()).unwrap();
    cache.insert(numerical).unwrap();
    assert_eq!(cache.len(), 2);
    exact.point = CachedPoint::Exact(BTreeMap::from([(s, Atom::new())]));
    exact.kind = PointKind::ContourDetour;
    cache.insert(exact).unwrap();
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::new())]));
    let mut query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 1).unwrap(), 20).unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    assert_eq!(
        cache
            .best(&query, &policy, p)
            .unwrap()
            .unwrap()
            .boundary
            .kind,
        PointKind::Physical
    );
    query.include_detours = true;
    assert_eq!(
        cache.best(&query, &policy, p).unwrap().unwrap().cost,
        p.real(0)
    );
    let mut numeric_only = BoundaryCache::default();
    let mut boundary = entry(&identity, s, Atom::num(1), 40);
    boundary.point = numerical_point;
    numeric_only.insert(boundary).unwrap();
    query.include_detours = false;
    assert!(numeric_only.best(&query, &policy, p).unwrap().is_some());
    query.minimum_coordinate_bits = p.bits + 1;
    assert!(numeric_only.best(&query, &policy, p).unwrap().is_none());
}

#[test]
fn boundary_evidence_caps_repeated_transport_and_rejects_false_labels() {
    let (s, _, _, identity) = fixture();
    let p = Precision::decimal(90).unwrap();
    let solution = EpsilonSolution {
        point: p.i(1),
        leading: 0,
        coefficients: vec![vec![p.i(1)]],
        diagnostics: FlowDiagnostics {
            working_bits: p.bits,
            ..Default::default()
        },
        segments: Vec::new(),
        checkpoints: Vec::new(),
        verified_digits: Some(40),
        comparison_errors: vec![vec![p.real(0)]],
    };
    let evidence = BoundaryAccuracy::from_solution(
        &solution,
        12,
        "same 12-digit boundary at both working precisions",
    )
    .unwrap();
    assert_eq!(evidence.verified_digits(), 12);
    let mut unverified = solution.clone();
    unverified.verified_digits = None;
    assert!(BoundaryAccuracy::from_solution(&unverified, 40, "unverified").is_err());
    let mut boundary = entry(&identity, s, Atom::num(1), 40);
    boundary.range = EpsilonRange::new(0, 0).unwrap();
    boundary.coefficients = solution.coefficients;
    boundary.accuracy = evidence;
    let requested = boundary
        .as_epsilon_boundary(
            &Atom::new(),
            boundary.range,
            Precision::decimal(120).unwrap(),
        )
        .unwrap();
    assert_eq!(requested.point.re, p.real(0));
    assert_eq!(boundary.accuracy.verified_digits(), 12);
    boundary.accuracy =
        BoundaryAccuracy::supplied(40, p.bits, vec![vec![p.real(1)]], "bad claimed error").unwrap();
    assert!(matches!(boundary.validate(), Err(Error::Accuracy(_))));
    assert!(BoundaryAccuracy::supplied(40, 32, vec![], "too few working bits").is_err());
}

#[test]
fn exact_input_and_regular_point_checks_reject_invalid_cache_entries() {
    let (s, epsilon, _, identity) = fixture();
    let mut invalid = entry(&identity, s, Atom::num(Float::with_val(80, 1)), 20);
    assert!(invalid.validate().is_err());
    for matrix in [Atom::var(s).pow(-1), Atom::var(s).pow(Atom::num((1, 2)))] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![matrix]])]),
        };
        invalid.identity = BoundaryIdentity::new(
            &system,
            &[parse!("boundary_cache::I(1)")],
            &Atom::num(1),
            Prescription::PlusI0,
            "physical-sheet",
        )
        .unwrap();
        invalid.point = CachedPoint::Exact(BTreeMap::from([(s, Atom::new())]));
        assert!(invalid.validate().is_err());
    }
}

#[test]
fn accepted_trajectory_endpoints_are_mapped_and_need_endpoint_evidence() {
    let (s, epsilon, _, identity) = fixture();
    let x = symbol!("cache_trajectory::x");
    let p = Precision::decimal(70).unwrap();
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::num(2) + Atom::num(3) * Atom::var(x))]),
    };
    // The physical system dY/ds = epsilon*Y pulls back to dY/dx=3*epsilon*Y.
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon)]])]),
    };
    let expanded =
        EpsilonSystem::from_differential_system(&system.pullback(&path).unwrap(), epsilon, 2)
            .unwrap();
    let boundary = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]],
    };
    let options = FlowOptions {
        digits: 30,
        guard_digits: 40,
        series_order: 16,
        ..Default::default()
    };
    let solution = expanded
        .compile(p, &Default::default())
        .unwrap()
        .transport(&boundary, &[p.i(1)], &options, &RunContext::default(), true)
        .unwrap();
    assert!(!solution.segments.is_empty());
    assert_eq!(solution.verified_digits, None);
    let mut cache = BoundaryCache::default();
    assert_eq!(
        cache
            .insert_trajectory(&identity, &path, &solution, |_, _| Ok(None))
            .unwrap(),
        0
    );
    let count = cache
        .insert_trajectory(&identity, &path, &solution, |index, segment| {
            let high = Precision::decimal(110).unwrap();
            let x = high.round(&segment.end);
            let expected = [
                high.i(1),
                high.mul(&high.i(3), &x),
                high.scale(&high.powi(&x, 2), 9, 2),
            ];
            let computed = solution.evaluate_segment(index, &segment.end)?;
            let errors = computed
                .iter()
                .zip(expected)
                .map(|(row, expected)| vec![high.norm(&high.sub(&row[0], &expected))])
                .collect();
            Ok(Some((
                PointKind::Physical,
                BoundaryAccuracy::supplied(
                    40,
                    segment.working_bits,
                    errors,
                    "independent analytic exp(3*epsilon*x) coefficients",
                )?,
            )))
        })
        .unwrap();
    assert_eq!(count, solution.segments.len());
    for entry in cache.entries() {
        assert!(matches!(entry.point, CachedPoint::Derived { .. }));
        let coordinates = entry.point.evaluate(p).unwrap();
        let x = p.scale(&p.sub(&coordinates[&s], &p.i(2)), 1, 3);
        assert!(p.close(&entry.coefficients[1][0], &p.mul(&p.i(3), &x), 40));
    }
}

#[test]
fn derived_checkpoint_preserves_multivariate_chart_despite_coordinate_rounding() {
    let (s, t, x, epsilon) = symbol!(
        "derived_chart::s",
        "derived_chart::t",
        "derived_chart::x",
        "derived_chart::eps"
    );
    let p = Precision::decimal(40).unwrap();
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]]), (t, vec![vec![Atom::new()]])]),
    };
    let identity = BoundaryIdentity::new(
        &system,
        &[parse!("derived_chart::Y")],
        &Atom::num(1),
        Prescription::PlusI0,
        "constant solution, one chart",
    )
    .unwrap();
    let large = Atom::num(10).pow(100);
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([
            (s, &large + Atom::var(x)),
            (t, &large + Atom::var(x) + Atom::num(1)),
        ]),
    };
    let boundary = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.i(1)]],
    };
    let options = FlowOptions {
        digits: 20,
        guard_digits: 20,
        series_order: 8,
        ..Default::default()
    };
    let solution =
        EpsilonSystem::from_differential_system(&system.pullback(&path).unwrap(), epsilon, 0)
            .unwrap()
            .compile(p, &Default::default())
            .unwrap()
            .transport(&boundary, &[p.i(1)], &options, &RunContext::default(), true)
            .unwrap();
    let mut cache = RustFlowCache::default();
    cache
        .insert_trajectory(&identity, &path, &solution, |_, segment| {
            Ok(Some((
                PointKind::Physical,
                BoundaryAccuracy::supplied(
                    25,
                    segment.working_bits,
                    vec![vec![p.real(0)]],
                    "analytic constant solution",
                )?,
            )))
        })
        .unwrap();
    for cached in cache.entries() {
        assert!(matches!(cached.point, CachedPoint::Derived { .. }));
        let exact = cached.point.restart_coordinates().unwrap();
        assert!(
            (&exact[&t] - &exact[&s] - Atom::num(1))
                .together()
                .cancel()
                .is_zero()
        );
        // Independent rounding loses the unit displacement at this scale. The
        // restart chart must use the retained exact image, not these display values.
        let rounded = cached.point.evaluate(p).unwrap();
        assert_eq!(rounded[&s], rounded[&t]);
    }
}

#[test]
fn batch_insertion_is_transactional_and_replaces_through_its_index() {
    let (s, _, _, identity) = fixture();
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&identity, s, Atom::num(1), 20)).unwrap();
    let valid = entry(&identity, s, Atom::num(2), 30);
    let mut invalid = entry(&identity, s, Atom::num(3), 30);
    invalid.coefficients.clear();
    assert!(cache.insert_many(vec![valid.clone(), invalid]).is_err());
    assert_eq!(cache.len(), 1);
    cache
        .insert_many(vec![valid, entry(&identity, s, Atom::num(1), 40)])
        .unwrap();
    assert_eq!(cache.len(), 2);
    assert_eq!(cache.entries()[0].accuracy.verified_digits(), 40);
    cache.insert(entry(&identity, s, Atom::num(1), 25)).unwrap();
    assert_eq!(cache.len(), 2);
    assert_eq!(cache.entries()[0].accuracy.verified_digits(), 40);
}

#[test]
fn persistent_cache_roundtrips_native_atoms_and_mpfr_and_rejects_corruption() {
    let (s, _, _, identity) = fixture();
    let mut cache = BoundaryCache::default();
    let boundary = entry(&identity, s, Atom::num((1, 7)), 40);
    cache.insert(boundary.clone()).unwrap();
    let mut derived = boundary.clone();
    derived.point = CachedPoint::Derived {
        coordinates: BTreeMap::from([(s, Atom::num((1, 7)))]),
        working_bits: boundary.accuracy.working_bits(),
        provenance: "exact image of chosen parameter".into(),
    };
    cache.insert(derived).unwrap();
    let directory = std::env::temp_dir().join(format!(
        "amflow-physical-boundary-cache-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&directory);
    assert!(BoundaryCache::load(&directory).unwrap().is_empty());
    cache.save(&directory).unwrap();
    let restored = BoundaryCache::load(&directory).unwrap();
    assert_eq!(restored.len(), 2);
    assert_eq!(restored.entries()[0].identity.key(), identity.key());
    assert_eq!(restored.entries()[0].coefficients, boundary.coefficients);
    assert_eq!(
        restored.entries()[0]
            .point
            .rounded_coordinates_as_exact()
            .unwrap(),
        boundary.point.rounded_coordinates_as_exact().unwrap()
    );
    assert!(
        matches!(&restored.entries()[1].point, CachedPoint::Derived { provenance, .. } if provenance == "exact image of chosen parameter")
    );
    assert_eq!(std::fs::read_dir(&directory).unwrap().count(), 1);
    let file = directory.join("physical-boundaries.bin");
    let pristine = std::fs::read(&file).unwrap();
    let mut corrupt = pristine.clone();
    let last = corrupt.len() - 1;
    corrupt[last] ^= 1;
    std::fs::write(&file, corrupt).unwrap();
    assert!(matches!(
        BoundaryCache::load(&directory),
        Err(Error::Cache(_))
    ));
    let mut incompatible = pristine;
    incompatible[0] ^= 1;
    std::fs::write(&file, incompatible).unwrap();
    assert!(matches!(
        BoundaryCache::load(&directory),
        Err(Error::Cache(_))
    ));
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn persistent_implementation_or_dependency_mismatch_rejects_before_exact_hit() {
    // Mirror only the envelope, leaving the encoded identities, values and
    // payload checksum unchanged; this isolates compatibility rejection.
    #[derive(serde::Serialize, serde::Deserialize)]
    struct Envelope {
        version: u32,
        implementation: String,
        dependencies: String,
        digest: String,
        payload: Vec<u8>,
    }
    let (s, _, system, identity) = fixture();
    let flow = symbolica_amflow::RustFlow::new(
        system,
        &[parse!("boundary_cache::I(1)")],
        &Atom::num(1),
        Prescription::PlusI0,
        "physical-sheet",
    )
    .unwrap();
    assert_eq!(flow.identity().key(), identity.key());
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&identity, s, Atom::num(1), 40)).unwrap();
    let directory = std::env::temp_dir().join(format!(
        "amflow-boundary-compatibility-before-hit-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&directory);
    cache.save(&directory).unwrap();
    let file = directory.join("physical-boundaries.bin");
    let pristine = std::fs::read(&file).unwrap();
    let magic = b"AMFLOW-BOUNDARIES\0\x08";
    let (envelope, consumed): (Envelope, _) = bincode::serde::decode_from_slice(
        pristine.strip_prefix(magic).unwrap(),
        bincode::config::standard(),
    )
    .unwrap();
    assert_eq!(consumed + magic.len(), pristine.len());
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let destination = BTreeMap::from([(s, Atom::num(1))]);
    let mut loaded = RustFlowCache::load(&directory).unwrap();
    let hit = flow
        .evaluate_to(
            &mut loaded,
            &destination,
            EpsilonRange::new(0, 1).unwrap(),
            &FlowOptions::default(),
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert!(
        hit.transport.is_none(),
        "pristine snapshot has an exact hit"
    );
    for implementation_mismatch in [true, false] {
        let changed = Envelope {
            version: envelope.version,
            implementation: if implementation_mismatch {
                format!("obsolete-{}", envelope.implementation)
            } else {
                envelope.implementation.clone()
            },
            dependencies: if implementation_mismatch {
                envelope.dependencies.clone()
            } else {
                format!("obsolete-{}", envelope.dependencies)
            },
            digest: envelope.digest.clone(),
            payload: envelope.payload.clone(),
        };
        let mut bytes = magic.to_vec();
        bytes.extend(bincode::serde::encode_to_vec(&changed, bincode::config::standard()).unwrap());
        std::fs::write(&file, bytes).unwrap();
        let reused = std::cell::Cell::new(false);
        let result = RustFlowCache::load(&directory).and_then(|mut restored| {
            reused.set(true);
            flow.evaluate_to(
                &mut restored,
                &destination,
                EpsilonRange::new(0, 1).unwrap(),
                &FlowOptions::default(),
                &RunContext::default(),
                &policy,
            )
        });
        assert!(matches!(result, Err(Error::Cache(_))));
        assert!(
            !reused.get(),
            "incompatible snapshot reached exact-hit reuse"
        );
    }
    for obsolete_schema in [true, false] {
        let mut payload = envelope.payload.clone();
        if !obsolete_schema {
            payload.push(0);
        }
        let changed = Envelope {
            version: if obsolete_schema { 7 } else { envelope.version },
            implementation: envelope.implementation.clone(),
            dependencies: envelope.dependencies.clone(),
            digest: blake3::hash(&payload).to_hex().to_string(),
            payload,
        };
        let mut bytes = magic.to_vec();
        bytes.extend(bincode::serde::encode_to_vec(&changed, bincode::config::standard()).unwrap());
        std::fs::write(&file, bytes).unwrap();
        let error = RustFlowCache::load(&directory).unwrap_err();
        assert!(matches!(error, Error::Cache(_)));
        if !obsolete_schema {
            assert!(
                error
                    .to_string()
                    .contains("trailing boundary-cache payload bytes")
            );
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn obsolete_residual_alias_snapshot_is_rejected_before_admitting_values() {
    // This is a real e733 solver-produced bank: its old verified API returned
    // epsilon1=0 although the exact endpoint is -999. Never import its values.
    let bytes = include_bytes!("../fixtures/regressions/residual-alias-cache-v5.bin");
    assert_eq!(bytes.len(), 1346);
    let directory = std::env::temp_dir().join(format!(
        "amflow-obsolete-residual-alias-snapshot-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join("physical-boundaries.bin"), bytes).unwrap();
    let result = RustFlowCache::load(&directory);
    std::fs::remove_dir_all(directory).unwrap();
    assert!(
        matches!(result, Err(Error::Cache(_))),
        "obsolete values were admitted: {result:?}"
    );
}
