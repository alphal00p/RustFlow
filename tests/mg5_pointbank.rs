//! Actual MG5 physical-grid data, compared with the pinned original DiffExp.
//! This exercises a closed two-master subsystem, not the complete EW amplitude.
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use symbolica::prelude::*;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::transport_cache::{
    BoundaryAccuracy, CachedBoundary, CachedPoint, EpsilonRange, PointKind, ScaledDistance,
};
use symbolica_amflow::{FlowOptions, Precision, Prescription, RunContext, RustFlow, RustFlowCache};

#[derive(Deserialize)]
struct DecimalComplex {
    real: String,
    imaginary: String,
}

#[derive(Deserialize)]
struct Point {
    coordinates: BTreeMap<String, String>,
    dependent_coordinates: BTreeMap<String, String>,
    epsilon_coefficients_by_master: Vec<Vec<DecimalComplex>>,
}

#[derive(Deserialize)]
struct Input {
    #[serde(flatten)]
    point: Point,
    input_verified_digits: u32,
    conservative_absolute_error_per_coefficient: String,
}

#[derive(Deserialize)]
struct Connection {
    epsilon_order: i64,
    matrix: Vec<Vec<String>>,
}

#[derive(Deserialize)]
struct System {
    coordinate_names: Vec<String>,
    basis: Vec<String>,
    coordinate_connections: BTreeMap<String, Connection>,
    physical_sheet: String,
}

#[derive(Deserialize)]
struct Fixture {
    system: System,
    input: Input,
    reference_points: BTreeMap<String, Point>,
}

fn exact(text: &str) -> Atom {
    Atom::parse(text, "mg5_pointbank", Default::default()).unwrap()
}

fn variable(name: &str) -> Symbol {
    match name {
        "s" => symbol!("mg5_pointbank::s"),
        "t" => symbol!("mg5_pointbank::t"),
        "b" => symbol!("mg5_pointbank::b"),
        other => panic!("unexpected independent coordinate {other}"),
    }
}

fn coordinates(point: &Point) -> BTreeMap<Symbol, Atom> {
    let values: BTreeMap<_, _> = point
        .coordinates
        .iter()
        .map(|(name, value)| (variable(name), exact(value)))
        .collect();
    assert_eq!(values.len(), 3);
    // u is dependent and must never become an extra, mismatched PDE coordinate.
    assert!(
        (&values[&variable("b")]
            - &values[&variable("s")]
            - &values[&variable("t")]
            - exact(&point.dependent_coordinates["u"]))
        .together()
        .cancel()
        .is_zero()
    );
    values
}

fn assert_point(actual: &CachedPoint, expected: &BTreeMap<Symbol, Atom>) {
    let actual = actual.rounded_coordinates_as_exact().unwrap();
    assert!(actual.keys().eq(expected.keys()));
    for (coordinate, value) in expected {
        assert!((&actual[coordinate] - value).together().cancel().is_zero());
    }
}

fn assert_reference(actual: &CachedBoundary, reference: &Point, input_digits: u32) {
    assert_point(&actual.point, &coordinates(reference));
    assert!(actual.accuracy.verified_digits() >= 20);
    assert!(actual.accuracy.verified_digits() <= input_digits);
    assert_eq!(actual.coefficients.len(), 5);
    let p = Precision {
        bits: actual.accuracy.working_bits(),
    };
    for (order, values) in actual.coefficients.iter().enumerate() {
        assert_eq!(values.len(), 2);
        for (master, value) in values.iter().enumerate() {
            let expected = &reference.epsilon_coefficients_by_master[master][order];
            let expected = p.parse(&expected.real, &expected.imaginary).unwrap();
            let error = p.norm(&p.sub(value, &expected));
            assert!(
                error < p.tolerance(20),
                "master {}, epsilon order {order}: absolute difference {error}",
                master + 1
            );
        }
    }
}

struct TemporaryBank(PathBuf);

impl TemporaryBank {
    fn new() -> Self {
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        Self(std::env::temp_dir().join(format!(
            "rustflow-mg5-pointbank-{}-{stamp}",
            std::process::id()
        )))
    }
}

impl Drop for TemporaryBank {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn actual_mg5_grid_boundary_grows_reuses_and_restarts_the_physical_point_bank() {
    let fixture: Fixture = serde_json::from_str(include_str!(
        "../fixtures/mg5-higgs-ew/np-two-master-transport.json"
    ))
    .unwrap();
    assert_eq!(fixture.system.coordinate_names, ["s", "t", "b"]);
    assert_eq!(fixture.input.input_verified_digits, 24);
    let epsilon = symbol!("mg5_pointbank::eps");
    let system = KinematicSystem {
        epsilon,
        derivatives: fixture
            .system
            .coordinate_connections
            .iter()
            .map(|(name, connection)| {
                (
                    variable(name),
                    connection
                        .matrix
                        .iter()
                        .map(|row| {
                            row.iter()
                                .map(|entry| {
                                    Atom::var(epsilon).pow(connection.epsilon_order) * exact(entry)
                                })
                                .collect()
                        })
                        .collect(),
                )
            })
            .collect(),
    };
    let engine = RustFlow::new(
        system,
        &fixture
            .system
            .basis
            .iter()
            .map(|s| exact(s))
            .collect::<Vec<_>>(),
        &Atom::num(1),
        Prescription::PlusI0,
        &fixture.system.physical_sheet,
    )
    .unwrap();
    let p = Precision::decimal(100).unwrap();
    let range = EpsilonRange::new(0, 4).unwrap();
    let input = &fixture.input.point;
    assert_eq!(input.epsilon_coefficients_by_master.len(), 2);
    assert!(
        input
            .epsilon_coefficients_by_master
            .iter()
            .all(|v| v.len() == 5)
    );
    // Fixture layout is master x epsilon; RustFlow's layout is epsilon x master.
    // Decimal strings go straight into MPFR, and coordinates stay exact Atoms.
    let coefficients = (0..5)
        .map(|order| {
            input
                .epsilon_coefficients_by_master
                .iter()
                .map(|values| {
                    p.parse(&values[order].real, &values[order].imaginary)
                        .unwrap()
                })
                .collect()
        })
        .collect();
    let inherited_error = p
        .parse(
            &fixture.input.conservative_absolute_error_per_coefficient,
            "0",
        )
        .unwrap()
        .re;
    let mut cache = RustFlowCache::default();
    cache
        .insert(CachedBoundary {
            identity: engine.identity().clone(),
            point: CachedPoint::Exact(coordinates(input)),
            kind: PointKind::Physical,
            range,
            coefficients,
            accuracy: BoundaryAccuracy::supplied(
                fixture.input.input_verified_digits,
                p.bits,
                vec![vec![inherited_error; 2]; 5],
                "Original MG5 NP_EW1 grid estimate rounded upward to 1e-27; evidence capped at 24 digits, independently of stored mantissa length",
            )
            .unwrap(),
        })
        .unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |source: &CachedBoundary, target: &CachedPoint| {
            let from = source.point.rounded_coordinates_as_exact()?;
            let to = target.rounded_coordinates_as_exact()?;
            for name in ["t", "b"] {
                if !(&from[&variable(name)] - &to[&variable(name)])
                    .together()
                    .cancel()
                    .is_zero()
                {
                    return Ok(false);
                }
            }
            for values in [&from, &to] {
                let s = p.eval(&values[&variable("s")], &Default::default())?;
                if s.im != p.real(0) || s.re <= p.real(1) {
                    return Ok(false);
                }
            }
            Ok(true)
        },
    };
    let options = FlowOptions {
        digits: 20,
        ..Default::default()
    };
    let context = RunContext::default();
    let checkpoint = &fixture.reference_points["checkpoint"];
    let destination = &fixture.reference_points["destination_direct"];
    let first = engine
        .evaluate_to(
            &mut cache,
            &coordinates(checkpoint),
            range,
            &options,
            &context,
            &policy,
        )
        .unwrap();
    assert_point(&first.starting_point, &coordinates(input));
    assert!(first.transport.is_some());
    assert!(first.inserted_points > 0);
    assert!(cache.len() > 1);
    assert_reference(
        &first.boundary,
        checkpoint,
        fixture.input.input_verified_digits,
    );

    // Save before the second query: a resumed bank must select the newly computed
    // checkpoint too, rather than quietly restart from the original grid seed.
    let directory = TemporaryBank::new();
    cache.save(&directory.0).unwrap();
    let second = engine
        .evaluate_to(
            &mut cache,
            &coordinates(destination),
            range,
            &options,
            &context,
            &policy,
        )
        .unwrap();
    assert_point(&second.starting_point, &coordinates(checkpoint));
    assert!(second.transport.is_some());
    assert_reference(
        &second.boundary,
        destination,
        fixture.input.input_verified_digits,
    );
    assert_reference(
        &second.boundary,
        &fixture.reference_points["destination_from_checkpoint"],
        fixture.input.input_verified_digits,
    );
    let repeated = engine
        .evaluate_to(
            &mut cache,
            &coordinates(destination),
            range,
            &options,
            &context,
            &policy,
        )
        .unwrap();
    assert!(repeated.transport.is_none());
    assert_eq!(repeated.inserted_points, 0);
    assert_eq!(repeated.boundary.coefficients, second.boundary.coefficients);
    assert_point(&repeated.boundary.point, &coordinates(destination));

    let mut restarted_cache = RustFlowCache::load(&directory.0).unwrap();
    let resumed = engine
        .evaluate_to(
            &mut restarted_cache,
            &coordinates(destination),
            range,
            &options,
            &context,
            &policy,
        )
        .unwrap();
    assert_point(&resumed.starting_point, &coordinates(checkpoint));
    assert!(resumed.transport.is_some());
    assert_reference(
        &resumed.boundary,
        destination,
        fixture.input.input_verified_digits,
    );
    assert_eq!(resumed.boundary.coefficients, second.boundary.coefficients);
}
