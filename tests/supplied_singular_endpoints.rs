use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, SquareRoot};
use symbolica_amflow::kinematics::{KinematicPath, KinematicSystem};
use symbolica_amflow::singular_endpoint::*;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

fn a(text: &str) -> Atom {
    Atom::parse(text, "supplied_endpoint_tests", Default::default()).unwrap()
}
fn s() -> Symbol {
    symbol!("supplied_endpoint_tests::s")
}
fn z() -> Symbol {
    symbol!("supplied_endpoint_tests::z")
}
fn eps() -> Symbol {
    symbol!("supplied_endpoint_tests::eps")
}
fn root() -> Symbol {
    symbol!("supplied_endpoint_tests::r")
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

#[test]
fn rooted_public_endpoint_restarts_and_keeps_regular_anchors() {
    let flow = rooted(&[&["1/r"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let request = request(range, Some(RootSheet::Principal));
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/4",
        Some(RootSheet::Principal),
        vec![vec![p.exp(&p.i(1))]],
        70,
    );
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(!result.cache_hit);
    assert!(result.matching_transport.is_some());
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 30));
    assert_eq!(cache.endpoint_entries().len(), 1);
    assert!(
        cache
            .entries()
            .iter()
            .all(|entry| entry.point.evaluate(p).unwrap()[&s()].re > p.real(0))
    );
    let directory =
        std::env::temp_dir().join(format!("amflow-endpoint-restart-{}", std::process::id()));
    cache.save(&directory).unwrap();
    let mut loaded = RustFlowCache::load(&directory).unwrap();
    std::fs::remove_dir_all(directory).unwrap();
    let hit = flow
        .evaluate_endpoint(
            &mut loaded,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(hit.cache_hit);
    assert_eq!(hit.boundary.coefficients, result.boundary.coefficients);
    let regular = flow
        .evaluate_to(
            &mut loaded,
            &BTreeMap::from([(s(), a("3/32"))]),
            request.chart.matching_germ.as_ref().unwrap(),
            range,
            &options(),
            &RunContext::default(),
            &cost(),
        )
        .unwrap();
    assert!(regular.starting_point.evaluate(p).unwrap()[&s()].re > p.real(0));
    assert!(p.close(
        &regular.boundary.coefficients[0][0],
        &p.exp(&p.scale(
            &p.pow(
                &p.rational(&Rational::from((3, 32))),
                &p.scale(&p.i(1), 1, 2)
            ),
            2,
            1
        )),
        25
    ));
}

#[test]
fn opposite_roots_and_log_winding_are_distinct_terminal_identities() {
    let flow = rooted(&[&["1/r"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    for (sheet, sign) in [(RootSheet::Principal, 1), (RootSheet::Opposite, -1)] {
        let mut cache = bank(
            flow.identity(),
            range,
            "1/16",
            Some(sheet),
            vec![vec![p.exp(&p.scale(&p.i(sign), 1, 2))]],
            70,
        );
        for winding in [0, 1, -1] {
            let mut request = request(range, Some(sheet));
            request.chart.winding = winding;
            let result = flow
                .evaluate_endpoint(
                    &mut cache,
                    &request,
                    &options(),
                    &RunContext::default(),
                    &cost(),
                    &admit,
                )
                .unwrap();
            assert!(!result.cache_hit);
            assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 30));
        }
        assert_eq!(cache.endpoint_entries().len(), 3);
        let mut changed = request(range, Some(sheet));
        changed.chart.homotopy = "another explicitly admitted terminal homotopy".into();
        let result = flow
            .evaluate_endpoint(
                &mut cache,
                &changed,
                &options(),
                &RunContext::default(),
                &cost(),
                &admit,
            )
            .unwrap();
        assert!(!result.cache_hit);
        assert_eq!(cache.endpoint_entries().len(), 4);
    }
}

#[test]
fn coupled_half_power_logs_have_finite_stable_limits() {
    let flow = rooted(&[
        &["1/(2*s)", "1/r", "0"],
        &["0", "0", "0"],
        &["1/s", "0", "0"],
    ]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let log = p.log(&p.rational(&Rational::from((1, 16))));
    let values = vec![vec![
        p.scale(&log, 1, 4),
        p.i(1),
        p.scale(&p.sub(&log, &p.i(2)), 1, 2),
    ]];
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        values,
        70,
    );
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request(range, Some(RootSheet::Principal)),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    for (value, expected) in result.boundary.coefficients[0].iter().zip([0, 1, 0]) {
        assert!(p.close(value, &p.i(expected), 30));
    }
}

#[test]
fn coupled_epsilon_hierarchy_and_negative_leading_power_are_preserved() {
    let flow = rooted(&[&["eps/r"]]);
    let range = EpsilonRange::new(-2, 1).unwrap();
    let p = Precision::decimal(90).unwrap();
    let values = vec![
        vec![p.i(1)],
        vec![p.scale(&p.i(1), 1, 2)],
        vec![p.scale(&p.i(1), 1, 8)],
        vec![p.scale(&p.i(1), 1, 48)],
    ];
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        values,
        70,
    );
    let request = request(range, Some(RootSheet::Principal));
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert_eq!(result.boundary.range, range);
    for (i, row) in result.boundary.coefficients.iter().enumerate() {
        assert!(p.close(&row[0], &p.i(i64::from(i == 0)), 30));
    }
    let mut lower = request;
    lower.range = EpsilonRange::new(-2, -1).unwrap();
    let hit = flow
        .evaluate_endpoint(
            &mut cache,
            &lower,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(hit.cache_hit);
    assert_eq!(hit.boundary.coefficients.len(), 2);
    assert_eq!(hit.boundary.matching_boundary.coefficients.len(), 2);
}

#[test]
fn infinity_is_an_exact_chart_not_a_fake_finite_cache_point() {
    let flow = ordinary(&[&["1/s^2"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut request = request(range, None);
    request.chart.path.coordinates.insert(s(), a("1/z"));
    let mut cache = bank(
        flow.identity(),
        range,
        "16",
        None,
        vec![vec![p.exp(&p.rational(&Rational::from((-1, 16))))]],
        70,
    );
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 30));
    assert_eq!(cache.len(), 1);
    assert_eq!(cache.endpoint_entries().len(), 1);
}

#[test]
fn divergent_finite_epsilon_coefficients_and_unstable_zero_sectors_are_refused_atomically() {
    let p = Precision::decimal(90).unwrap();
    for (matrix, range, values) in [
        (
            vec![vec!["eps/s"]],
            EpsilonRange::new(0, 1).unwrap(),
            vec![vec![p.i(1)], vec![p.zero()]],
        ),
        (
            vec![vec!["-1/s", "0"], vec!["0", "0"]],
            EpsilonRange::new(0, 0).unwrap(),
            vec![vec![p.zero(), p.i(1)]],
        ),
    ] {
        let rows = matrix.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let flow = ordinary(&rows);
        let mut cache = bank(flow.identity(), range, "1", None, values, 70);
        let error = flow
            .evaluate_endpoint(
                &mut cache,
                &request(range, None),
                &options(),
                &RunContext::default(),
                &cost(),
                &admit,
            )
            .err()
            .unwrap();
        assert!(matches!(error, Error::Accuracy(_)), "{error}");
        assert_eq!(cache.len(), 1);
        assert!(cache.endpoint_entries().is_empty());
    }
}

#[test]
fn endpoint_admission_and_cancellation_roll_back_regular_transport() {
    let flow = ordinary(&[&["1/(2*s)"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(flow.identity(), range, "1", None, vec![vec![p.i(1)]], 70);
    let request = request(range, None);
    let deny = |_: &CachedBoundary, _: &EndpointChart| Ok(false);
    assert!(
        flow.evaluate_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &deny
        )
        .is_err()
    );
    assert_eq!(cache.len(), 1);
    let context = RunContext::default();
    let cancel = |_: &CachedBoundary, _: &EndpointChart| {
        context.cancellation.cancel();
        Ok(true)
    };
    assert!(matches!(
        flow.evaluate_endpoint(&mut cache, &request, &options(), &context, &cost(), &cancel),
        Err(Error::Cancelled)
    ));
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
}

#[test]
fn source_holes_inside_matching_disk_and_old_schema_are_rejected() {
    let flow = ordinary(&[&["1/(s-1/32)"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(flow.identity(), range, "1/16", None, vec![vec![p.i(1)]], 70);
    let error = flow
        .evaluate_endpoint(
            &mut cache,
            &request(range, None),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .err()
        .unwrap();
    assert!(matches!(error, Error::InvalidInput(_)), "{error}");
    let directory =
        std::env::temp_dir().join(format!("amflow-endpoint-old-schema-{}", std::process::id()));
    cache.save(&directory).unwrap();
    let file = directory.join("physical-boundaries.bin");
    let mut bytes = std::fs::read(&file).unwrap();
    bytes[b"AMFLOW-BOUNDARIES\0".len()] = 5;
    std::fs::write(file, bytes).unwrap();
    assert!(matches!(
        RustFlowCache::load(&directory),
        Err(Error::Cache(_))
    ));
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn independently_raised_precision_and_order_agree_without_upgrading_input_evidence() {
    let flow = rooted(&[&["1/r"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut outcomes = Vec::new();
    for (guard, order) in [(35, 8), (55, 8), (35, 64)] {
        let mut cache = bank(
            flow.identity(),
            range,
            "1/16",
            Some(RootSheet::Principal),
            vec![vec![p.exp(&p.scale(&p.i(1), 1, 2))]],
            45,
        );
        let mut request = request(range, Some(RootSheet::Principal));
        request.options.series_order = order;
        let options = FlowOptions {
            digits: 30,
            guard_digits: guard,
            ..options()
        };
        let result = flow
            .evaluate_endpoint(
                &mut cache,
                &request,
                &options,
                &RunContext::default(),
                &cost(),
                &admit,
            )
            .unwrap();
        assert!(result.boundary.accuracy.verified_digits() <= 45);
        assert_eq!(result.boundary.accuracy.input_verified_digits(), 45);
        assert!(
            result
                .boundary
                .accuracy
                .provenance()
                .contains("independent (bits, order) profiles")
        );
        assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 40));
        outcomes.push(result.boundary.coefficients[0][0].clone());
    }
    assert!(
        outcomes
            .windows(2)
            .all(|pair| p.close(&pair[0], &pair[1], 40))
    );
}

#[test]
fn amplified_input_uncertainty_is_not_repaired_by_more_working_precision() {
    let flow = ordinary(&[&["-100"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.exp(&p.rational(&Rational::from((-25, 4))))]],
        20,
    );
    let error = flow
        .evaluate_endpoint(
            &mut cache,
            &request(range, None),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .err()
        .unwrap();
    assert!(matches!(error, Error::Accuracy(_)), "{error}");
    assert_eq!(cache.len(), 1);
    assert!(cache.endpoint_entries().is_empty());
}

#[test]
fn prescribed_matching_route_requires_its_own_admission() {
    use std::cell::Cell;
    let flow = ordinary(&[&["1/s"]])
        .with_prescribed_continuation(PhysicalContinuation {
            prescriptions: vec![symbolica_amflow::contour::PolynomialPrescription {
                polynomial: a("s"),
                prescription: Prescription::PlusI0,
            }],
            unprescribed_side: Prescription::PlusI0,
            domain: "upper detour then positive radial approach".into(),
        })
        .unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(flow.identity(), range, "-1", None, vec![vec![p.i(-1)]], 70);
    let request = request(range, None);
    assert!(matches!(
        flow.evaluate_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit
        ),
        Err(Error::InvalidInput(_))
    ));
    let deny = |_: &CachedBoundary,
                _: &CachedPoint,
                _: &symbolica_amflow::physical_transport::PhysicalRoute| Ok(false);
    assert!(
        flow.evaluate_prescribed_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &deny,
            &admit
        )
        .is_err()
    );
    assert_eq!(cache.len(), 1);
    let called = Cell::new(0);
    let route = |_: &CachedBoundary,
                 _: &CachedPoint,
                 path: &symbolica_amflow::physical_transport::PhysicalRoute| {
        called.set(called.get() + 1);
        assert!(!path.crossings.is_empty());
        Ok(true)
    };
    let result = flow
        .evaluate_prescribed_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &route,
            &admit,
        )
        .unwrap();
    assert!(called.get() > 0);
    assert!(p.close(&result.boundary.coefficients[0][0], &p.zero(), 20));
    assert!(result.matching_transport.is_some());
}

#[test]
fn wrong_root_germ_and_lift_limit_cannot_reuse_an_endpoint() {
    let flow = rooted(&[&["1/r"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        vec![vec![p.exp(&p.scale(&p.i(1), 1, 2))]],
        70,
    );
    let good = request(range, Some(RootSheet::Principal));
    flow.evaluate_endpoint(
        &mut cache,
        &good,
        &options(),
        &RunContext::default(),
        &cost(),
        &admit,
    )
    .unwrap();
    let wrong = request(range, Some(RootSheet::Opposite));
    assert!(
        flow.evaluate_endpoint(
            &mut cache,
            &wrong,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit
        )
        .is_err()
    );
    assert_eq!(cache.endpoint_entries().len(), 1);
    let mut limited = good;
    limited.chart.homotopy = "distinct request requiring preparation".into();
    limited.options.max_lift_dimension = 1;
    assert!(matches!(
        flow.evaluate_endpoint(
            &mut cache,
            &limited,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit
        ),
        Err(Error::Limit(_))
    ));
    assert_eq!(cache.endpoint_entries().len(), 1);
}

#[test]
fn tiny_later_resonant_logarithms_are_not_discarded_before_admission() {
    // r=1/s on the positive sheet. In the rational lift the auxiliary
    // exponent -1 resonates at k=1, creating a physical log with coefficient
    // 10^-300, below the ordinary recurrence's numerical trimming threshold.
    let system = AlgebraicKinematicSystem {
        epsilon: eps(),
        derivatives: BTreeMap::from([(
            s(),
            vec![vec![a("0"), a("r/10^300")], vec![a("0"), a("0")]],
        )]),
        roots: vec![SquareRoot {
            symbol: root(),
            radicand: a("1/s^2"),
        }],
    };
    let flow = RustFlow::new_algebraic(
        system,
        &[a("I0"), a("I1")],
        &a("1"),
        Prescription::PlusI0,
        "positive sheet",
    )
    .unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        vec![vec![p.zero(), p.i(1)]],
        70,
    );
    let error = flow
        .evaluate_endpoint(
            &mut cache,
            &request(range, Some(RootSheet::Principal)),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .err()
        .unwrap();
    assert!(matches!(error, Error::Accuracy(_)), "{error}");
    assert!(cache.endpoint_entries().is_empty());
}

#[test]
fn rounded_root_cancellation_does_not_prove_absence_of_a_divergent_log() {
    let system = AlgebraicKinematicSystem {
        epsilon: eps(),
        derivatives: BTreeMap::from([(
            s(),
            vec![vec![a("0"), a("r-(1+1/10^300)/s")], vec![a("0"), a("0")]],
        )]),
        roots: vec![SquareRoot {
            symbol: root(),
            radicand: a("1/s^2"),
        }],
    };
    let flow = RustFlow::new_algebraic(
        system,
        &[a("I0"), a("I1")],
        &a("1"),
        Prescription::PlusI0,
        "fixed positive ray",
    )
    .unwrap();
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    for sheet in [RootSheet::Principal, RootSheet::Opposite] {
        let mut cache = bank(
            flow.identity(),
            range,
            "1/16",
            Some(sheet),
            vec![vec![p.zero(), p.i(1)]],
            70,
        );
        let error = flow
            .evaluate_endpoint(
                &mut cache,
                &request(range, Some(sheet)),
                &options(),
                &RunContext::default(),
                &cost(),
                &admit,
            )
            .err()
            .unwrap();
        assert!(matches!(error, Error::Accuracy(_)), "{error}");
        assert_eq!(cache.len(), 1);
        assert!(cache.endpoint_entries().is_empty());
    }
}

#[test]
fn incompatible_normalization_cannot_reuse_terminal_or_regular_evidence() {
    let flow = ordinary(&[&["1/(2*s)"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.scale(&p.i(1), 1, 4)]],
        70,
    );
    let request = request(range, None);
    flow.evaluate_endpoint(
        &mut cache,
        &request,
        &options(),
        &RunContext::default(),
        &cost(),
        &admit,
    )
    .unwrap();
    let other = RustFlow::new(
        KinematicSystem {
            epsilon: eps(),
            derivatives: BTreeMap::from([(s(), vec![vec![a("1/(2*s)")]])]),
        },
        &[a("I0")],
        &a("2"),
        Prescription::PlusI0,
        "declared endpoint domain",
    )
    .unwrap();
    assert!(
        other
            .evaluate_endpoint(
                &mut cache,
                &request,
                &options(),
                &RunContext::default(),
                &cost(),
                &admit
            )
            .is_err()
    );
    assert_eq!(cache.endpoint_entries().len(), 1);
}

#[test]
fn tiny_exact_positive_exponents_are_not_false_sampled_resonances() {
    let flow = ordinary(&[&["0", "0"], &["0", "1/(10^300*s)"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        None,
        vec![vec![p.i(1), p.i(1)]],
        70,
    );
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request(range, None),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 30));
    assert!(p.close(&result.boundary.coefficients[0][1], &p.zero(), 30));
}

#[test]
fn canonical_registered_root_endpoint_reuses_the_same_evidence_path() {
    let system = symbolica_amflow::algebraic::CanonicalAlgebraicSystem::new(
        eps(),
        &[s()],
        &[a("1+r")],
        &[vec![vec![a("1")]]],
        vec![SquareRoot {
            symbol: root(),
            radicand: a("s"),
        }],
    )
    .unwrap();
    let flow = RustFlow::new_canonical(
        system,
        &[a("I0")],
        &a("1"),
        Prescription::PlusI0,
        "positive real root",
    )
    .unwrap();
    let range = EpsilonRange::new(0, 2).unwrap();
    let p = Precision::decimal(90).unwrap();
    let log = p.log(&p.rational(&Rational::from((5, 4))));
    let mut cache = bank(
        flow.identity(),
        range,
        "1/16",
        Some(RootSheet::Principal),
        vec![
            vec![p.i(1)],
            vec![log.clone()],
            vec![p.scale(&p.mul(&log, &log), 1, 2)],
        ],
        70,
    );
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request(range, Some(RootSheet::Principal)),
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    for (i, row) in result.boundary.coefficients.iter().enumerate() {
        assert!(p.close(&row[0], &p.i(i64::from(i == 0)), 30));
    }
}

#[test]
fn endpoint_hierarchy_limits_are_checked_before_source_preparation() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    let flow = ordinary(&[&["1/s"]]);
    let mut cache = RustFlowCache::default();
    let mut request = request(EpsilonRange::new(0, 1023).unwrap(), None);
    let phases = Arc::new(AtomicUsize::new(0));
    let seen = phases.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |_| {
            seen.fetch_add(1, Ordering::SeqCst);
        })),
        ..RunContext::default()
    };
    let error = flow
        .evaluate_endpoint(&mut cache, &request, &options(), &context, &cost(), &admit)
        .err()
        .unwrap();
    assert!(matches!(error,Error::Limit(message) if message.contains("before source pullback")));
    assert_eq!(phases.load(Ordering::SeqCst), 0);
    request.range = EpsilonRange {
        leading: i32::MIN,
        last: i32::MAX,
    };
    assert!(matches!(
        flow.evaluate_endpoint(&mut cache, &request, &options(), &context, &cost(), &admit),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(phases.load(Ordering::SeqCst), 0);
    assert!(cache.is_empty());
}

#[test]
fn exact_complex_matching_parameter_preserves_the_registered_sheet() {
    let flow = rooted(&[&["1/r"]]);
    let range = EpsilonRange::new(0, 0).unwrap();
    let p = Precision::decimal(90).unwrap();
    let mut request = request(range, Some(RootSheet::Principal));
    request.chart.matching_parameter =
        Atom::num(Complex::new(Rational::zero(), Rational::from((1, 16))));
    let matching = p.complex(0, 1);
    let matching = p.scale(&matching, 1, 16);
    let value = p.exp(&p.scale(&p.pow(&matching, &p.scale(&p.i(1), 1, 2)), 2, 1));
    let mut cache = RustFlowCache::default();
    cache
        .insert(CachedBoundary {
            identity: flow.identity().clone(),
            point: request.chart.matching_point().unwrap(),
            kind: PointKind::Physical,
            range,
            coefficients: vec![vec![value]],
            accuracy: BoundaryAccuracy::supplied(
                70,
                p.bits,
                vec![vec![p.real(0)]],
                "analytic exp(2 sqrt(s)) on principal complex sheet",
            )
            .unwrap(),
        })
        .unwrap();
    let result = flow
        .evaluate_endpoint(
            &mut cache,
            &request,
            &options(),
            &RunContext::default(),
            &cost(),
            &admit,
        )
        .unwrap();
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 30));
    assert_eq!(
        result
            .boundary
            .matching_boundary
            .point
            .restart_coordinates()
            .unwrap()[&s()],
        request.chart.matching_parameter
    );
}
