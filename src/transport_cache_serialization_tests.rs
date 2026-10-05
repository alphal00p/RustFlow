//! Snapshot serialization checks use the uncached native exporter as an oracle.
use super::*;
use std::process::Command;

struct UncachedAtoms;
impl AtomEncoder for UncachedAtoms {
    fn atom(&mut self, atom: &Atom) -> Result<StoredAtom> {
        atom_bytes(atom)
    }
}

fn payload(cache: &RustFlowCache, atoms: &mut impl AtomEncoder) -> Vec<u8> {
    bincode::serde::encode_to_vec(
        StoredCache::encode(cache, atoms).unwrap(),
        bincode::config::standard(),
    )
    .unwrap()
}

fn fixture() -> RustFlowCache {
    let (x, epsilon, root) = symbol!(
        "snapshot_export_test::x",
        "snapshot_export_test::epsilon",
        "snapshot_export_test::root"
    );
    let p = Precision::decimal(100).unwrap();
    let one = Atom::num(1);
    let canonical = CanonicalAlgebraicSystem::new(
        epsilon,
        &[x],
        &[Atom::var(x)],
        &[vec![vec![one.clone()]]],
        vec![SquareRoot {
            symbol: root,
            radicand: Atom::var(x),
        }],
    )
    .unwrap();
    let basis = [parse!("snapshot_export_test::I(1)")];
    let algebraic = BoundaryIdentity::with_canonical_conditions(
        &canonical,
        &basis,
        &one,
        Prescription::PlusI0,
        "snapshot fixture",
        &[Atom::var(x)],
    )
    .unwrap()
    .with_prescribed_continuation(PhysicalContinuation {
        prescriptions: vec![crate::contour::PolynomialPrescription {
            polynomial: Atom::var(x),
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
        domain: "snapshot fixture affine continuation".into(),
    })
    .unwrap();
    let dense = BoundaryIdentity::new(
        &KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(x, vec![vec![Atom::var(epsilon) / Atom::var(x)]])]),
        },
        &basis,
        &one,
        Prescription::MinusI0,
        "snapshot fixture dense",
    )
    .unwrap();
    let exact = BTreeMap::from([(x, Atom::num((1, 7)))]);
    let points = [
        CachedPoint::Exact(exact.clone()),
        CachedPoint::Derived {
            coordinates: exact,
            working_bits: p.bits,
            provenance: "exact image of a rounded path parameter".into(),
        },
        CachedPoint::Numerical {
            coordinates: BTreeMap::from([(
                x,
                p.parse(
                    "0.142857142857142857142857142857142857142857142857142857",
                    "0",
                )
                .unwrap(),
            )]),
            working_bits: p.bits,
            provenance: "independently rounded numerical coordinate".into(),
        },
    ];
    let mut cache = RustFlowCache::default();
    for (identity, algebraic) in [(dense, false), (algebraic, true)] {
        for point in &points {
            let point = if algebraic {
                point
                    .clone()
                    .with_root_germ(RootGerm {
                        sheets: BTreeMap::from([(root, RootSheet::Opposite)]),
                    })
                    .unwrap()
            } else {
                point.clone()
            };
            cache
                .insert(CachedBoundary {
                    identity: identity.clone(),
                    point,
                    kind: PointKind::Physical,
                    range: EpsilonRange::new(-1, 0).unwrap(),
                    coefficients: vec![
                        vec![
                            p.parse(
                                "1.234567890123456789012345678901234567890123456789",
                                "-0.0000000000000000000000000000000000000123456789",
                            )
                            .unwrap(),
                        ],
                        vec![
                            p.parse(
                                "-0.333333333333333333333333333333333333333333333333",
                                "2.718281828459045235360287471352662497757247093699",
                            )
                            .unwrap(),
                        ],
                    ],
                    accuracy: BoundaryAccuracy::supplied(
                        40,
                        p.bits,
                        vec![vec![p.tolerance(60)], vec![p.tolerance(65)]],
                        "synthetic serialization fixture; no physics accuracy claim",
                    )
                    .unwrap(),
                })
                .unwrap();
        }
    }
    assert_eq!(cache.len(), 6);
    cache
}

fn directory(label: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("rustflow-snapshot-{label}-{}", std::process::id()))
}

// A subprocess isolates Symbolica's global polynomial/field registries from
// unrelated parallel tests, while still exercising the production native codec.
fn run_isolated(name: &str) -> bool {
    const CHILD: &str = "RUSTFLOW_SNAPSHOT_TEST_CHILD";
    if std::env::var(CHILD).as_deref() == Ok(name) {
        return false;
    }
    assert!(
        Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                &format!("transport_cache::serialization_tests::{name}"),
                "--nocapture",
                "--include-ignored"
            ])
            .env(CHILD, name)
            .status()
            .unwrap()
            .success()
    );
    true
}

#[test]
fn snapshot_payload_matches_uncached_native_export() {
    if run_isolated("snapshot_payload_matches_uncached_native_export") {
        return;
    }
    let cache = fixture();
    let _ = payload(&cache, &mut UncachedAtoms); // Finish lazy registration.
    let expected = payload(&cache, &mut UncachedAtoms);
    let mut atoms = SnapshotAtoms::default();
    let actual = payload(&cache, &mut atoms);
    assert_eq!(actual, expected);
    assert_eq!(blake3::hash(&actual), blake3::hash(&expected));
    assert!(
        atoms.0.len() < 25,
        "one memo spans both identities and every point"
    );
}

#[test]
fn snapshot_roundtrips_in_a_fresh_process_with_different_symbol_order() {
    if run_isolated("snapshot_roundtrips_in_a_fresh_process_with_different_symbol_order") {
        return;
    }
    const CHILD: &str = "RUSTFLOW_SNAPSHOT_IMPORT_DIRECTORY";
    if let Some(directory) = std::env::var_os(CHILD) {
        for i in 0..64 {
            let _ = Atom::parse(
                format!("unrelated{i}"),
                "snapshot_import_first",
                Default::default(),
            )
            .unwrap();
        }
        let loaded = RustFlowCache::load(Path::new(&directory)).unwrap();
        let expected = fixture();
        assert_eq!(loaded.len(), expected.len());
        for (a, b) in loaded.entries().iter().zip(expected.entries()) {
            assert_eq!(a.identity.key(), b.identity.key());
            assert_eq!(
                a.point.restart_coordinates().unwrap(),
                b.point.restart_coordinates().unwrap()
            );
            // Native Float serialization checks the MPFR precision as well as
            // values; the whole evidence object retains errors and provenance.
            assert_eq!(
                bincode::serde::encode_to_vec(&a.coefficients, bincode::config::standard())
                    .unwrap(),
                bincode::serde::encode_to_vec(&b.coefficients, bincode::config::standard())
                    .unwrap()
            );
            assert_eq!(
                serde_json::to_value(&a.accuracy).unwrap(),
                serde_json::to_value(&b.accuracy).unwrap()
            );
            assert_eq!(
                payload(
                    &RustFlowCache {
                        entries: vec![a.clone()],
                        ..Default::default()
                    },
                    &mut UncachedAtoms
                ),
                payload(
                    &RustFlowCache {
                        entries: vec![b.clone()],
                        ..Default::default()
                    },
                    &mut UncachedAtoms
                )
            );
        }
        return;
    }
    let path = directory("fresh-process");
    fixture().save(&path).unwrap();
    assert!(Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "transport_cache::serialization_tests::snapshot_roundtrips_in_a_fresh_process_with_different_symbol_order", "--nocapture"])
        .env(CHILD, &path)
        .status().unwrap().success());
    std::fs::remove_dir_all(path).unwrap();
}

fn polynomial_coefficient(name: &str, constant: i64) -> Atom {
    let atom = Atom::parse(name, "snapshot_growth", Default::default()).unwrap();
    let AtomView::Var(variable) = atom.as_view() else {
        panic!("test parameter must be a symbol")
    };
    (&atom + constant).set_coefficient_ring(variable.get_symbol())
}

fn finite_field_coefficient(prime: u64, value: u64) -> Atom {
    use symbolica::domains::finite_field::{FiniteFieldCore, Zp64};
    let field = Zp64::new(prime);
    let element = field.to_element(value);
    Atom::num(Coefficient::from_finite_field(field, element))
}

#[test]
fn snapshot_exports_remain_self_contained_as_symbol_state_grows() {
    if run_isolated("snapshot_exports_remain_self_contained_as_symbol_state_grows") {
        return;
    }
    const CHILD: &str = "RUSTFLOW_SNAPSHOT_GROWTH_BLOBS";
    if let Some(file) = std::env::var_os(CHILD) {
        // Deliberately occupy registry indices with different data. Construct
        // the expected imported objects only AFTER import so missing native
        // state cannot be masked by pre-registering their actual lists/fields.
        let _ = polynomial_coefficient("different_child_list", 19);
        let _ = polynomial_coefficient("another_child_list", 23);
        let _ = finite_field_coefficient(1019, 5);
        let _ = finite_field_coefficient(1021, 7);
        let data = std::fs::read(file).unwrap();
        let (blobs, consumed): (Vec<StoredAtom>, _) =
            bincode::serde::decode_from_slice(&data, bincode::config::standard()).unwrap();
        assert_eq!(consumed, data.len());
        let actual: Vec<_> = blobs
            .iter()
            .map(|bytes| atom_read(bytes).unwrap())
            .collect();
        let expected = vec![
            polynomial_coefficient("old_parameter", 3),
            finite_field_coefficient(1013, 41),
            finite_field_coefficient(1009, 37),
            polynomial_coefficient("new_parameter", 11),
            polynomial_coefficient("old_parameter", 3),
            finite_field_coefficient(1009, 37),
            parse!("snapshot_growth::a+1"),
        ];
        assert_eq!(actual, expected);
        return;
    }

    let original = parse!("snapshot_growth::a+1");
    let mut atoms = SnapshotAtoms::default();
    let before = atoms.atom(&original).unwrap();
    let old_polynomial = polynomial_coefficient("old_parameter", 3);
    let old_field = finite_field_coefficient(1009, 37);
    let old_polynomial_bytes = atoms.atom(&old_polynomial).unwrap();
    let old_field_bytes = atoms.atom(&old_field).unwrap();

    // Grow polynomial and field registries BETWEEN native export calls.
    // Concurrent changes inside a single native export are not tested here.
    let new_polynomial = polynomial_coefficient("new_parameter", 11);
    let new_field = finite_field_coefficient(1013, 41);
    let new_polynomial_bytes = atoms.atom(&new_polynomial).unwrap();
    let new_field_bytes = atoms.atom(&new_field).unwrap();
    assert_eq!(atoms.atom(&original).unwrap(), before);
    assert_eq!(atoms.atom(&old_polynomial).unwrap(), old_polynomial_bytes);
    assert_eq!(atoms.atom(&old_field).unwrap(), old_field_bytes);
    let mixed_blobs = vec![
        old_polynomial_bytes,
        new_field_bytes,
        old_field_bytes,
        new_polynomial_bytes,
        atoms.atom(&old_polynomial).unwrap(),
        atoms.atom(&old_field).unwrap(),
        before,
    ];
    let path = directory("registry-growth");
    std::fs::create_dir_all(&path).unwrap();
    let file = path.join("mixed-native-blobs.bin");
    std::fs::write(
        &file,
        bincode::serde::encode_to_vec(mixed_blobs, bincode::config::standard()).unwrap(),
    )
    .unwrap();
    assert!(Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "transport_cache::serialization_tests::snapshot_exports_remain_self_contained_as_symbol_state_grows", "--nocapture"])
        .env(CHILD, &file)
        .status().unwrap().success());
    std::fs::remove_dir_all(path).unwrap();

    // Plain unrelated symbol registration does not grow the polynomial/field
    // tables. Check that it can coexist with memo hits and new expression exports.
    let start = std::sync::Barrier::new(2);
    std::thread::scope(|scope| {
        scope.spawn(|| {
            start.wait();
            for i in 0..256 {
                let _ = Atom::parse(
                    format!("growing{i}"),
                    "snapshot_concurrent_growth",
                    Default::default(),
                )
                .unwrap();
            }
        });
        start.wait();
        for i in 0..256 {
            assert_eq!(
                atom_read(&atoms.atom(&original).unwrap()).unwrap(),
                original
            );
            let fresh = &original + i;
            assert_eq!(atom_read(&atoms.atom(&fresh).unwrap()).unwrap(), fresh);
        }
    });
    assert_eq!(
        SnapshotAtoms::default().atom(&original).unwrap(),
        atom_bytes(&original).unwrap()
    );
}

#[test]
fn rejected_snapshot_preserves_previous_file_and_failed_rename_cleans_temporary() {
    if run_isolated("rejected_snapshot_preserves_previous_file_and_failed_rename_cleans_temporary")
    {
        return;
    }
    let path = directory("atomic-failure");
    let mut cache = fixture();
    cache.save(&path).unwrap();
    let before = std::fs::read(path.join(FILE)).unwrap();
    cache.entries[0].coefficients.clear();
    assert!(cache.save(&path).is_err());
    assert_eq!(std::fs::read(path.join(FILE)).unwrap(), before);
    assert_eq!(std::fs::read_dir(&path).unwrap().count(), 1);
    std::fs::remove_file(path.join(FILE)).unwrap();
    std::fs::create_dir(path.join(FILE)).unwrap();
    assert!(fixture().save(&path).is_err());
    assert!(path.join(FILE).is_dir());
    assert_eq!(std::fs::read_dir(&path).unwrap().count(), 1);
    std::fs::remove_dir_all(path).unwrap();
}

/// Same-process microbenchmark: real 48/61 canonical matrix topology and
/// physical source coordinates, explicitly synthetic values/evidence. No saved
/// acceptance bank is read and no source-compatibility check is bypassed.
#[test]
#[ignore = "serialization profile; run alone in release mode with --nocapture"]
fn profile_canonical_snapshot_exports() {
    if run_isolated("profile_canonical_snapshot_exports") {
        return;
    }
    use crate::gg_hg::{HiggsJetIntegralSystem, PluginFamilyKind};
    use std::time::Instant;

    let p = Precision::decimal(90).unwrap();
    let mut cache = RustFlowCache::default();
    let mut matrix_entries = 0;
    let mut matrix_zeros = 0;
    for kind in [PluginFamilyKind::Planar, PluginFamilyKind::Nonplanar] {
        let namespace = format!("snapshot_profile_{}", kind.id());
        let system = HiggsJetIntegralSystem::load(kind, &namespace).unwrap();
        for atom in system
            .canonical_system()
            .constant_matrices()
            .iter()
            .flatten()
            .flatten()
        {
            matrix_entries += 1;
            matrix_zeros += usize::from(atom.is_zero());
        }
        for configuration in system.configurations(&namespace).unwrap() {
            // Three coordinate-origin variants make 48 entries without claiming
            // transported physics values or importing somebody else's evidence.
            for point in [
                CachedPoint::Exact(configuration.start.clone()),
                CachedPoint::Derived {
                    coordinates: configuration.start.clone(),
                    working_bits: p.bits,
                    provenance: "synthetic serialization profile path image".into(),
                },
                CachedPoint::Numerical {
                    coordinates: configuration
                        .start
                        .iter()
                        .map(|(&s, a)| (s, p.eval(a, &Default::default()).unwrap()))
                        .collect(),
                    working_bits: p.bits,
                    provenance: "synthetic serialization profile rounded coordinate".into(),
                },
            ] {
                cache
                    .insert(CachedBoundary {
                        identity: system.transport().identity().clone(),
                        point: point.with_root_germ(configuration.germ.clone()).unwrap(),
                        kind: PointKind::Physical,
                        range: EpsilonRange::new(0, 4).unwrap(),
                        coefficients: vec![
                            vec![
                                p.parse(
                                    "1.23456789012345678901234567890123456789",
                                    "-0.987654321098765432109876543210987654321"
                                )
                                .unwrap();
                                kind.dimension()
                            ];
                            5
                        ],
                        accuracy: BoundaryAccuracy::supplied(
                            40,
                            p.bits,
                            vec![vec![p.tolerance(60); kind.dimension()]; 5],
                            "synthetic serialization timing fixture; not a physics result",
                        )
                        .unwrap(),
                    })
                    .unwrap();
            }
        }
    }
    assert_eq!(cache.len(), 48);
    assert_eq!((matrix_entries, matrix_zeros), (362019, 357430));
    let initial_symbols = symbolica::state::State::symbol_iter().count();
    let mut profiles = Vec::new();
    for extra_symbols in [0, 1000] {
        for i in 0..extra_symbols {
            let _ = Atom::parse(
                format!("registered{i}"),
                "snapshot_profile_unrelated",
                Default::default(),
            )
            .unwrap();
        }
        // Warm both paths first. Each timed snapshot creates a fresh memo and
        // validates every entry, just like production save. The payload includes
        // bincode encoding; file I/O is measured separately below.
        let mut memo = SnapshotAtoms::default();
        let expected = payload(&cache, &mut UncachedAtoms);
        assert_eq!(expected, payload(&cache, &mut memo));
        let unique_atoms = memo.0.len();
        let mut uncached_seconds = Vec::new();
        let mut memoized_seconds = Vec::new();
        for round in 0..5 {
            for memoized in [round % 2 == 0, round % 2 != 0] {
                let start = Instant::now();
                let actual = if memoized {
                    payload(&cache, &mut SnapshotAtoms::default())
                } else {
                    payload(&cache, &mut UncachedAtoms)
                };
                let seconds = start.elapsed().as_secs_f64();
                assert_eq!(actual, expected);
                if memoized {
                    memoized_seconds.push(seconds);
                } else {
                    uncached_seconds.push(seconds);
                }
            }
        }
        profiles.push(serde_json::json!({
            "extra_unrelated_symbols": extra_symbols,
            "registered_symbols": symbolica::state::State::symbol_iter().count(),
            "distinct_exported_atoms": unique_atoms,
            "payload_bytes": expected.len(),
            "uncached_payload_seconds": uncached_seconds,
            "memoized_payload_seconds": memoized_seconds,
            "payloads_byte_identical": true,
        }));
    }
    let path = directory("canonical-profile");
    let mut save_seconds = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        cache.save(&path).unwrap();
        save_seconds.push(start.elapsed().as_secs_f64());
    }
    let restored = RustFlowCache::load(&path).unwrap();
    assert_eq!(restored.len(), 48);
    assert_eq!(
        payload(&restored, &mut SnapshotAtoms::default()),
        payload(&cache, &mut SnapshotAtoms::default())
    );
    std::fs::remove_dir_all(path).unwrap();
    println!(
        "{}",
        serde_json::json!({
            "scope": "real canonical topology, synthetic numerical values; isolated process, not notebook runtime",
            "entries": cache.len(),
            "matrix_entries": matrix_entries,
            "matrix_zeros": matrix_zeros,
            "initial_registered_symbols": initial_symbols,
            "source_digest": env!("PORT_SOURCE_DIGEST"),
            "dependency_digest": env!("DEPENDENCY_SOURCE_DIGEST"),
            "profiles": profiles,
            "memoized_atomic_save_seconds_after_registry_growth": save_seconds,
            "production_roundtrip_preserved_payload": true,
        })
    );
}
