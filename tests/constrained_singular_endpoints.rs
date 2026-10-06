use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, SquareRoot};
use symbolica_amflow::kinematics::{KinematicPath, KinematicSystem};
use symbolica_amflow::singular_endpoint::*;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

fn a(text: &str) -> Atom {
    Atom::parse(text, "constrained_endpoint_tests", Default::default()).unwrap()
}
fn s() -> Symbol {
    symbol!("constrained_endpoint_tests::s")
}
fn z() -> Symbol {
    symbol!("constrained_endpoint_tests::z")
}
fn eps() -> Symbol {
    symbol!("constrained_endpoint_tests::eps")
}
fn root() -> Symbol {
    symbol!("constrained_endpoint_tests::r")
}
fn ordinary(matrix: &[&[&str]]) -> RustFlow {
    let system = KinematicSystem {
        epsilon: eps(),
        derivatives: BTreeMap::from([(
            s(),
            matrix
                .iter()
                .map(|row| row.iter().map(|value| a(value)).collect())
                .collect(),
        )]),
    };
    let basis = (0..matrix.len())
        .map(|i| a(&format!("I{i}")))
        .collect::<Vec<_>>();
    RustFlow::new(
        system,
        &basis,
        &Atom::num(1),
        Prescription::PlusI0,
        "declared endpoint domain",
    )
    .unwrap()
}
fn rooted(matrix: &[&[&str]]) -> RustFlow<AlgebraicKinematicSystem> {
    let system = AlgebraicKinematicSystem {
        epsilon: eps(),
        derivatives: BTreeMap::from([(
            s(),
            matrix
                .iter()
                .map(|row| row.iter().map(|value| a(value)).collect())
                .collect(),
        )]),
        roots: vec![SquareRoot {
            symbol: root(),
            radicand: a("s"),
        }],
    };
    let basis = (0..matrix.len())
        .map(|i| a(&format!("I{i}")))
        .collect::<Vec<_>>();
    RustFlow::new_algebraic(
        system,
        &basis,
        &Atom::num(1),
        Prescription::PlusI0,
        "declared endpoint domain",
    )
    .unwrap()
}
fn request(range: EpsilonRange, sheet: Option<RootSheet>) -> EndpointRequest {
    EndpointRequest {
        chart: EndpointChart {
            path: KinematicPath {
                parameter: z(),
                coordinates: BTreeMap::from([(s(), a("z"))]),
            },
            matching_parameter: a("1/16"),
            matching_germ: sheet.map(|sheet| RootGerm {
                sheets: BTreeMap::from([(root(), sheet)]),
            }),
            winding: 0,
            homotopy: "radial approach within the declared branch".into(),
        },
        range,
        options: EndpointOptions {
            max_lift_dimension: 32,
            series_order: 32,
        },
    }
}
fn bank(
    identity: &BoundaryIdentity,
    range: EpsilonRange,
    coordinate: &str,
    sheet: Option<RootSheet>,
    values: Vec<Vec<ComplexFloat>>,
    digits: u32,
) -> RustFlowCache {
    let p = Precision::decimal(90).unwrap();
    let point = CachedPoint::Exact(BTreeMap::from([(s(), a(coordinate))]));
    let point = sheet
        .map(|sheet| {
            point
                .clone()
                .with_root_germ(RootGerm {
                    sheets: BTreeMap::from([(root(), sheet)]),
                })
                .unwrap()
        })
        .unwrap_or(point);
    let errors = values
        .iter()
        .map(|row| vec![p.real(0); row.len()])
        .collect();
    let mut bank = RustFlowCache::default();
    bank.insert(CachedBoundary {
        identity: identity.clone(),
        point,
        kind: PointKind::Physical,
        range,
        coefficients: values,
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            errors,
            "independent analytic supplied boundary",
        )
        .unwrap(),
    })
    .unwrap();
    bank
}
fn cost() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}
fn admit(_: &CachedBoundary, _: &EndpointChart) -> Result<bool> {
    Ok(true)
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 35,
        series_order: 40,
        ..FlowOptions::default()
    }
}

use symbolica_amflow::asymptotic::{
    AsymptoticSelector, ExactAsymptoticConstraints, ExactAsymptoticRelation,
};
fn constraints(selectors: &[(usize, &str, usize, &str)]) -> EndpointConstraints {
    EndpointConstraints {
        asymptotic: ExactAsymptoticConstraints {
            relations: selectors
                .iter()
                .map(|&(component, power, log_power, value)| {
                    ExactAsymptoticRelation::coefficient(
                        AsymptoticSelector {
                            component,
                            power: a(power),
                            log_power,
                        },
                        a(value),
                    )
                })
                .collect(),
            provenance: "explicit analytic physical-sector declaration".into(),
        },
        limits: Default::default(),
    }
}
fn run(
    flow: &RustFlow,
    cache: &mut RustFlowCache,
    request: &EndpointRequest,
    constraints: &EndpointConstraints,
) -> Result<EndpointResult> {
    flow.evaluate_constrained_endpoint(
        cache,
        request,
        constraints,
        &options(),
        &RunContext::default(),
        &cost(),
        &admit,
    )
}
#[test]
fn logarithmic_constraint_keeps_uncertain_free_amplitude_and_regular_anchor() {
    let flow = ordinary(&[&["0", "1/s"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let req = request(range, None);
    let constraints = constraints(&[(0, "0", 1, "0")]);
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(3), p.zero()]],
        65,
    );
    assert!(
        flow.evaluate_endpoint(
            &mut cache,
            &req,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit
        )
        .is_err()
    );
    let result = run(&flow, &mut cache, &req, &constraints).unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(3), 40));
    assert!(result.boundary.coefficients[0][1].is_zero());
    assert!(result.boundary.accuracy.comparison_errors()[0][0] > p.real(0));
    assert!(result.boundary.accuracy.verified_digits() <= 65);
    assert_eq!(
        result.boundary.matching_boundary.coefficients[0],
        vec![p.i(3), p.zero()]
    );
    assert_eq!(cache.len(), 1);
    assert_eq!(cache.endpoint_entries().len(), 1);
}
#[test]
fn exact_correlated_poles_and_affine_finite_value_are_supported() {
    let flow = ordinary(&[&["-1/(2*s)", "1/(2*s)"], &["1/(2*s)", "-1/(2*s)"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let req = request(range, None);
    let constraints = constraints(&[(0, "-1", 0, "0"), (0, "0", 0, "5")]);
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(5), p.i(5)]],
        60,
    );
    let result = run(&flow, &mut cache, &req, &constraints).unwrap();
    assert!(
        result.boundary.coefficients[0]
            .iter()
            .all(|v| p.close(v, &p.i(5), 40))
    );
    assert!(result.boundary.accuracy.verified_digits() <= 60);
}
#[test]
fn every_regular_component_must_match_the_asserted_affine_space() {
    let flow = ordinary(&[&["0", "1/s"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let constraints = constraints(&[(0, "0", 1, "0")]);
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(3), p.i(1)]],
        60,
    );
    let error = run(&flow, &mut cache, &request(range, None), &constraints)
        .err()
        .unwrap();
    assert!(matches!(error, Error::Accuracy(_)), "{error}");
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
}
#[test]
fn tiny_exact_divergent_rhs_and_missing_constraints_are_never_numerical_zero() {
    let flow = ordinary(&[&["0", "1/s"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(3), p.zero()]],
        60,
    );
    for constraint in [
        constraints(&[(0, "0", 1, "1/10^300")]),
        constraints(&[(0, "0", 0, "3")]),
    ] {
        assert!(run(&flow, &mut cache, &request(range, None), &constraint).is_err());
        assert!(cache.endpoint_entries().is_empty());
    }
}
#[test]
fn coupled_epsilon_constraints_include_negative_leading_coefficients() {
    let flow = ordinary(&[&["0", "eps/s"], &["0", "0"]]);
    let range = EpsilonRange::new(-1, 1).unwrap();
    let p = Precision::decimal(90).unwrap();
    let values = vec![
        vec![p.i(2), p.zero()],
        vec![p.i(3), p.zero()],
        vec![p.i(4), p.i(5)],
    ];
    let mut cache = bank(flow.identity(), range, "1/16", None, values.clone(), 65);
    let constraints = constraints(&[(2, "0", 1, "0"), (4, "0", 1, "0")]);
    let result = run(&flow, &mut cache, &request(range, None), &constraints).unwrap();
    for (actual, expected) in result
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(values.iter().flatten())
    {
        assert!(p.close(actual, expected, 40));
    }
    let shorter = request(EpsilonRange::new(-1, 0).unwrap(), None);
    assert!(run(&flow, &mut cache, &shorter, &constraints).is_err());
}
#[test]
fn native_atom_codec_preserves_constraints_and_exact_range_identity() {
    let flow = ordinary(&[&["0", "1/s"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(3), p.zero()]],
        65,
    );
    let constraints = constraints(&[(0, "0", 1, "0")]);
    let req = request(range, None);
    let result = run(&flow, &mut cache, &req, &constraints).unwrap();
    let directory = std::env::temp_dir().join(format!(
        "amflow-constrained-endpoint-{}",
        std::process::id()
    ));
    cache.save(&directory).unwrap();
    let mut loaded = RustFlowCache::load(&directory).unwrap();
    std::fs::remove_dir_all(directory).unwrap();
    let hit = run(&flow, &mut loaded, &req, &constraints).unwrap();
    assert!(hit.cache_hit);
    assert_eq!(hit.boundary.coefficients, result.boundary.coefficients);
    let mut other = constraints.clone();
    other.asymptotic.provenance = "different explicit declaration".into();
    assert!(!run(&flow, &mut loaded, &req, &other).unwrap().cache_hit);
    assert_eq!(loaded.endpoint_entries().len(), 2);
    assert!(
        flow.evaluate_endpoint(
            &mut loaded,
            &req,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit
        )
        .is_err()
    );
}
#[test]
fn constrained_root_lift_uses_proved_sheet_space_and_checks_matching_data() {
    let flow = rooted(&[&["1/r"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        vec![vec![p.exp(&p.rational(&Rational::from((1, 2))))]],
        65,
    );
    let result = flow.evaluate_constrained_endpoint(
        &mut cache,
        &request(range, Some(RootSheet::Principal)),
        &constraints(&[(0, "0", 0, "1")]),
        &options(),
        &RunContext::default(),
        &cost(),
        &admit,
    );
    let result = result.unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 40));
    assert_eq!(cache.endpoint_entries().len(), 1);
}
#[test]
fn short_normalized_prefix_refines_instead_of_guessing_missing_coefficients() {
    let flow = ordinary(&[&["0", "1/s^21"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(7), p.zero()]],
        65,
    );
    let mut req = request(range, None);
    req.options.series_order = 8;
    let result = run(&flow, &mut cache, &req, &constraints(&[(0, "-20", 0, "0")])).unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(7), 35));
    assert!(result.boundary.coefficients[0][1].is_zero());
}
#[test]
fn constrained_matching_uses_the_declared_common_log_winding() {
    let flow = ordinary(&[&["1/s", "1/s"], &["0", "1/s"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let z = p.rational(&Rational::from((1, 16)));
    let constraints = constraints(&[(0, "1", 0, "3"), (1, "1", 0, "2")]);
    for winding in [-1, 1] {
        let mut req = request(range, None);
        req.chart.winding = winding;
        let log = p.add(
            &p.log(&z),
            &p.scale(&p.log(&p.i(-1)), 2 * i64::from(winding), 1),
        );
        let values = vec![vec![
            p.mul(&z, &p.add(&p.i(3), &p.mul(&p.i(2), &log))),
            p.mul(&z, &p.i(2)),
        ]];
        let mut cache = bank(flow.identity(), range, "1/16", None, values, 65);
        let result = run(&flow, &mut cache, &req, &constraints).unwrap();
        assert!(
            result.boundary.coefficients[0]
                .iter()
                .all(ComplexFloat::is_zero)
        );
        req.chart.winding = 0;
        assert!(run(&flow, &mut cache, &req, &constraints).is_err());
        assert_eq!(cache.endpoint_entries().len(), 1);
    }
}
#[test]
fn independent_precision_and_order_profiles_agree_and_keep_input_cap() {
    let flow = ordinary(&[&["1/(1-s)", "0"], &["0", "-1/s"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut outcomes = Vec::new();
    for (guard, order) in [(35, 8), (55, 8), (35, 64)] {
        let mut cache = bank(
            flow.identity(),
            range,
            "1/16",
            None,
            vec![vec![p.rational(&Rational::from((16, 15))), p.zero()]],
            60,
        );
        let mut req = request(range, None);
        req.options.series_order = order;
        let mut opts = options();
        opts.guard_digits = guard;
        let result = flow
            .evaluate_constrained_endpoint(
                &mut cache,
                &req,
                &constraints(&[(1, "-1", 0, "0")]),
                &opts,
                &RunContext::default(),
                &cost(),
                &admit,
            )
            .unwrap();
        assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 30));
        assert!(result.boundary.accuracy.verified_digits() <= 60);
        outcomes.push(result.boundary.coefficients);
    }
    for value in &outcomes[1..] {
        assert!(p.close(&value[0][0], &outcomes[0][0][0], 30));
    }
}
#[test]
fn constrained_admission_cancellation_and_resource_refusal_are_transactional() {
    let flow = ordinary(&[&["0", "1/s"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let original = bank(
        flow.identity(),
        range,
        "1/4",
        None,
        vec![vec![p.i(3), p.zero()]],
        65,
    );
    let req = request(range, None);
    let declared = constraints(&[(0, "0", 1, "0")]);
    for cancel in [false, true] {
        let mut cache = original.clone();
        let context = RunContext::default();
        let admission = |_: &CachedBoundary, _: &EndpointChart| {
            if cancel {
                context.cancellation.cancel();
                Ok(true)
            } else {
                Ok(false)
            }
        };
        assert!(
            flow.evaluate_constrained_endpoint(
                &mut cache,
                &req,
                &declared,
                &options(),
                &context,
                &cost(),
                &admission
            )
            .is_err()
        );
        assert_eq!(cache.len(), 1);
        assert!(cache.endpoint_entries().is_empty());
        assert_eq!(
            cache.entries()[0].point.evaluate(p).unwrap()[&s()],
            p.rational(&Rational::from((1, 4)))
        );
    }
    let mut cache = original;
    let mut excessive = req.clone();
    assert!(EpsilonRange::new(i32::MIN, i32::MAX).is_err());
    excessive.range = EpsilonRange::new(0, 1023).unwrap();
    assert!(matches!(
        run(&flow, &mut cache, &excessive, &declared),
        Err(Error::Limit(_))
    ));
    let mut limited = declared.clone();
    limited.limits.max_dimension = 1;
    assert!(matches!(
        run(&flow, &mut cache, &req, &limited),
        Err(Error::Limit(_))
    ));
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
}
#[test]
fn previous_schema_six_is_explicitly_rejected() {
    let flow = ordinary(&[&["0"]]);
    let p = Precision::decimal(90).unwrap();
    let cache = bank(
        flow.identity(),
        EpsilonRange::new(0, 0).unwrap(),
        "1/16",
        None,
        vec![vec![p.i(1)]],
        60,
    );
    let directory = std::env::temp_dir().join(format!(
        "amflow-constrained-old-schema-{}",
        std::process::id()
    ));
    cache.save(&directory).unwrap();
    let file = directory.join("physical-boundaries.bin");
    let mut bytes = std::fs::read(&file).unwrap();
    bytes[b"AMFLOW-BOUNDARIES\0".len()] = 6;
    std::fs::write(&file, bytes).unwrap();
    assert!(matches!(
        RustFlowCache::load(&directory),
        Err(Error::Cache(_))
    ));
    std::fs::remove_dir_all(directory).unwrap();
}
#[test]
fn pure_epsilon_constraints_agree_with_the_independent_analytic_origin_owner() {
    use symbolica_amflow::algebraic::{AlgebraicSystem, AnalyticOriginOptions};
    use symbolica_amflow::diffexp::EpsilonSystem;
    let flow = ordinary(&[&["eps/s", "eps"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 3).unwrap();
    let mut matrices = vec![vec![vec![Atom::new(); 2]; 2]; 4];
    matrices[1] = vec![vec![a("1/s"), a("1")], vec![a("0"), a("0")]];
    let owner = AlgebraicSystem {
        system: EpsilonSystem {
            variable: s(),
            matrices,
        },
        roots: Vec::new(),
        nonzero_conditions: vec![a("s")],
    };
    let mut boundary = vec![vec![Atom::new(); 2]; 4];
    boundary[0][1] = Atom::one();
    let p = Precision::decimal(90).unwrap();
    let seed = owner
        .analytic_origin(
            &boundary,
            &Default::default(),
            &Default::default(),
            p,
            &p.rational(&Rational::from((1, 16))),
            &AnalyticOriginOptions {
                order: 40,
                check_digits: 40,
                ..Default::default()
            },
        )
        .unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        seed.boundary.coefficients,
        65,
    );
    let result = run(
        &flow,
        &mut cache,
        &request(range, None),
        &constraints(&[(2, "0", 1, "0"), (4, "0", 1, "0"), (6, "0", 1, "0")]),
    )
    .unwrap();
    for (actual, expected) in result
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(boundary.iter().flatten())
    {
        assert!(p.close(actual, &p.eval(expected, &Default::default()).unwrap(), 35));
    }
}
#[test]
fn constrained_prescribed_transport_retains_independent_route_and_terminal_admission() {
    let flow = ordinary(&[&["0", "1/s"], &["0", "0"]])
        .with_prescribed_continuation(PhysicalContinuation {
            prescriptions: vec![symbolica_amflow::contour::PolynomialPrescription {
                polynomial: a("s"),
                prescription: Prescription::PlusI0,
            }],
            unprescribed_side: Prescription::PlusI0,
            domain: "upper detour to constrained endpoint".into(),
        })
        .unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "-1",
        None,
        vec![vec![p.i(3), p.zero()]],
        65,
    );
    let req = request(range, None);
    let declared = constraints(&[(0, "0", 1, "0")]);
    assert!(matches!(
        run(&flow, &mut cache, &req, &declared),
        Err(Error::InvalidInput(_))
    ));
    let deny = |_: &CachedBoundary,
                _: &CachedPoint,
                _: &symbolica_amflow::physical_transport::PhysicalRoute| Ok(false);
    assert!(
        flow.evaluate_prescribed_constrained_endpoint(
            &mut cache,
            &req,
            &declared,
            &options(),
            &RunContext::default(),
            &cost(),
            &deny,
            &admit
        )
        .is_err()
    );
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
    let calls = std::cell::Cell::new(0);
    let allow = |_: &CachedBoundary,
                 _: &CachedPoint,
                 path: &symbolica_amflow::physical_transport::PhysicalRoute| {
        assert!(!path.crossings.is_empty());
        calls.set(calls.get() + 1);
        Ok(true)
    };
    let result = flow
        .evaluate_prescribed_constrained_endpoint(
            &mut cache,
            &req,
            &declared,
            &options(),
            &RunContext::default(),
            &cost(),
            &allow,
            &admit,
        )
        .unwrap();
    assert!(calls.get() > 0);
    assert!(result.matching_transport.is_some());
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(3), 30));
    let deny_endpoint = |_: &CachedBoundary, _: &EndpointChart| Ok(false);
    assert!(
        flow.evaluate_prescribed_constrained_endpoint(
            &mut cache,
            &req,
            &declared,
            &options(),
            &RunContext::default(),
            &cost(),
            &allow,
            &deny_endpoint
        )
        .is_err()
    );
}
#[test]
fn finite_epsilon_constraints_do_not_select_a_symbolic_dimensional_sector() {
    let flow = ordinary(&[&["eps/s"]]);
    let range = EpsilonRange::new(0, 1).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![
            vec![p.i(1)],
            vec![p.log(&p.rational(&Rational::from((1, 16))))],
        ],
        65,
    );
    for declared in [
        constraints(&[(0, "0", 0, "1")]),
        constraints(&[(0, "0", 0, "1"), (1, "0", 1, "0")]),
    ] {
        assert!(run(&flow, &mut cache, &request(range, None), &declared).is_err());
        assert!(cache.endpoint_entries().is_empty());
    }
}
#[test]
fn exact_constraints_follow_shifted_and_infinite_endpoint_charts() {
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    for (entry, map, matching) in [("1/(s-1)", "1+z", "17/16"), ("-1/s", "1/z", "16")] {
        let flow = ordinary(&[&["0", entry], &["0", "0"]]);
        let mut req = request(range, None);
        req.chart.path.coordinates.insert(s(), a(map));
        let mut cache = bank(
            flow.identity(),
            range,
            matching,
            None,
            vec![vec![p.i(3), p.zero()]],
            65,
        );
        let result = run(&flow, &mut cache, &req, &constraints(&[(0, "0", 1, "0")])).unwrap();
        assert!(p.close(&result.boundary.coefficients[0][0], &p.i(3), 35));
        assert!(result.boundary.coefficients[0][1].is_zero());
    }
}
#[test]
fn constraints_do_not_erase_retained_source_holes() {
    let flow = ordinary(&[&["1/(s-1/32)", "1/s"], &["0", "0"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.zero(), p.zero()]],
        65,
    );
    assert!(matches!(
        run(
            &flow,
            &mut cache,
            &request(range, None),
            &constraints(&[(0, "0", 1, "0")])
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(cache.endpoint_entries().is_empty());
}

#[test]
fn root_correlations_cancel_divergence_on_both_sheets_and_windings() {
    let flow = rooted(&[
        &["0", "1/(s*r)", "-1/(s*r)"],
        &["0", "0", "0"],
        &["0", "0", "0"],
    ]);
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    for sheet in [RootSheet::Principal, RootSheet::Opposite] {
        for winding in [0, 1] {
            let mut req = request(range, Some(sheet));
            req.chart.winding = winding;
            let mut cache = bank(
                flow.identity(),
                range,
                "1/16",
                Some(sheet),
                vec![vec![p.i(3), p.i(5), p.i(5)]],
                60,
            );
            let result = flow
                .evaluate_constrained_endpoint(
                    &mut cache,
                    &req,
                    &constraints(&[(0, "-1/2", 0, "0")]),
                    &options(),
                    &RunContext::default(),
                    &cost(),
                    &admit,
                )
                .unwrap();
            for (value, expected) in result.boundary.coefficients[0].iter().zip([3, 5, 5]) {
                assert!(p.close(value, &p.i(expected), 35));
            }
            assert!(result.boundary.accuracy.verified_digits() <= 60);
            assert_eq!(cache.len(), 1);
            assert_eq!(cache.endpoint_entries().len(), 1);
        }
    }
}

#[test]
fn root_log_constraint_rejects_nonzero_divergent_component() {
    let flow = rooted(&[&["0", "1/(s*r)"], &["0", "1/(2*s)"]]);
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    for (second, accepted) in [
        (p.zero(), true),
        (p.rational(&Rational::from((1, 4))), false),
    ] {
        let mut cache = bank(
            flow.identity(),
            range,
            "1/16",
            Some(RootSheet::Principal),
            vec![vec![p.i(3), second]],
            60,
        );
        let result = flow.evaluate_constrained_endpoint(
            &mut cache,
            &request(range, Some(RootSheet::Principal)),
            &constraints(&[(0, "0", 1, "0")]),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        );
        assert_eq!(result.is_ok(), accepted, "{:?}", result.as_ref().err());
        if let Ok(result) = result {
            assert!(p.close(&result.boundary.coefficients[0][0], &p.i(3), 35));
        } else {
            assert!(cache.endpoint_entries().is_empty());
        }
    }
}

#[test]
fn coupled_epsilon_root_constraints_survive_cache_restart_and_exact_range() {
    let flow = rooted(&[&["eps/r"]]);
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(-1, 1).unwrap();
    let req = request(range, Some(RootSheet::Opposite));
    let declared = constraints(&[(0, "-1", 0, "0")]);
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Opposite),
        vec![
            vec![p.i(1)],
            vec![p.rational(&Rational::from((-1, 2)))],
            vec![p.rational(&Rational::from((1, 8)))],
        ],
        60,
    );
    let result = flow
        .evaluate_constrained_endpoint(
            &mut cache,
            &req,
            &declared,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 35));
    for row in &result.boundary.coefficients[1..] {
        assert!(p.close(&row[0], &p.zero(), 35));
    }
    let directory = std::env::temp_dir().join(format!(
        "amflow-root-constrained-endpoint-{}",
        std::process::id()
    ));
    cache.save(&directory).unwrap();
    let mut loaded = RustFlowCache::load(&directory).unwrap();
    std::fs::remove_dir_all(directory).unwrap();
    let hit = flow
        .evaluate_constrained_endpoint(
            &mut loaded,
            &req,
            &declared,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(hit.cache_hit);
    assert_eq!(hit.boundary.coefficients, result.boundary.coefficients);
    let mut narrow = req.clone();
    narrow.range = EpsilonRange::new(-1, 0).unwrap();
    let recomputed = flow
        .evaluate_constrained_endpoint(
            &mut loaded,
            &narrow,
            &declared,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(!recomputed.cache_hit);
}

#[test]
fn root_refinement_profiles_agree_with_analytic_exponential() {
    let flow = rooted(&[&["1/r"]]);
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let mut outcomes = Vec::new();
    for (guard, order, sheet, sign) in [
        (35, 8, RootSheet::Principal, 1),
        (55, 8, RootSheet::Principal, 1),
        (35, 32, RootSheet::Opposite, -1),
    ] {
        let mut req = request(range, Some(sheet));
        req.options.series_order = order;
        let mut cache = bank(
            flow.identity(),
            range,
            "1/16",
            Some(sheet),
            vec![vec![p.scale(
                &p.exp(&p.rational(&Rational::from((sign, 2)))),
                3,
                1,
            )]],
            60,
        );
        let mut opts = options();
        opts.guard_digits = guard;
        let result = flow
            .evaluate_constrained_endpoint(
                &mut cache,
                &req,
                &constraints(&[(0, "-1", 0, "0")]),
                &opts,
                &RunContext::default(),
                &cost(),
                &admit,
            )
            .unwrap();
        assert!(p.close(&result.boundary.coefficients[0][0], &p.i(3), 30));
        assert!(result.boundary.accuracy.comparison_errors()[0][0] > p.real(0));
        outcomes.push(result.boundary.coefficients[0][0].clone());
    }
    for value in &outcomes[1..] {
        assert!(p.close(value, &outcomes[0], 30));
    }
}

#[test]
fn exact_half_power_declaration_uses_matching_germ_and_common_winding_on_complex_rays() {
    let flow = rooted(&[&["1/r"]]);
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    for coordinate in ["1/16", "𝑖/16"] {
        let square_root = p
            .eval(&a(&format!("({coordinate})^(1/2)")), &Default::default())
            .unwrap();
        for (sheet, sign) in [(RootSheet::Principal, 1), (RootSheet::Opposite, -1)] {
            for winding in [0, 1] {
                let amplitude = if winding == 0 { sign } else { -sign };
                let boundary = p.scale(&p.exp(&p.scale(&square_root, 2 * sign, 1)), amplitude, 1);
                let mut req = request(range, Some(sheet));
                req.chart.matching_parameter = a(coordinate);
                req.chart.winding = winding;
                let mut cache = bank(
                    flow.identity(),
                    range,
                    coordinate,
                    Some(sheet),
                    vec![vec![boundary]],
                    60,
                );
                let result = flow
                    .evaluate_constrained_endpoint(
                        &mut cache,
                        &req,
                        &constraints(&[(0, "1/2", 0, "2")]),
                        &options(),
                        &RunContext::default(),
                        &cost(),
                        &admit,
                    )
                    .unwrap();
                assert!(p.close(&result.boundary.coefficients[0][0], &p.i(amplitude), 35));
            }
        }
    }
}

#[test]
fn constrained_root_lift_dimension_and_cancellation_preserve_the_bank() {
    let flow = rooted(&[&["1/r"]]);
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let original = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        vec![vec![p.exp(&p.rational(&Rational::from((1, 2))))]],
        60,
    );
    let req = request(range, Some(RootSheet::Principal));
    let mut limited = constraints(&[(0, "-1", 0, "0")]);
    limited.limits.max_dimension = 1;
    let mut cache = original.clone();
    assert!(matches!(
        flow.evaluate_constrained_endpoint(
            &mut cache,
            &req,
            &limited,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit
        ),
        Err(Error::Limit(_))
    ));
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
    let context = RunContext::default();
    let cancel = |_: &CachedBoundary, _: &EndpointChart| {
        context.cancellation.cancel();
        Ok(true)
    };
    let mut cache = original;
    assert!(matches!(
        flow.evaluate_constrained_endpoint(
            &mut cache,
            &req,
            &constraints(&[(0, "-1", 0, "0")]),
            &options(),
            &context,
            &cost(),
            &cancel
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
}

#[test]
fn rational_coefficient_height_limit_does_not_limit_working_precision() {
    let flow = ordinary(&[&["0"]]);
    let p = Precision::decimal(120).unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let mut cache = bank(flow.identity(), range, "1/16", None, vec![vec![p.i(3)]], 60);
    let mut declared = constraints(&[(0, "-1", 0, "0")]);
    declared.limits.max_coefficient_bits = 32;
    let mut opts = options();
    opts.guard_digits = 80;
    let result = flow
        .evaluate_constrained_endpoint(
            &mut cache,
            &request(range, None),
            &declared,
            &opts,
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(3), 50));
}
