use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{reduction::Reduction, transport_cache::*, *};

fn prepared(weight: Atom) -> PreparedPhysicalFamily {
    prepared_with_constant(weight, Atom::num(-1))
}
fn prepared_with_constant(weight: Atom, constant: Atom) -> PreparedPhysicalFamily {
    let epsilon = symbol!("target_projection::eps");
    let family = IntegralFamily {
        name: "supplied_target_reduction".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant,
            scalar_products: vec![Atom::one()],
        }],
        physical_propagators: 1,
        epsilon,
        dimension: 4,
    };
    let backend = TableBackend {
        name: "explicit projection test relation".into(),
        reduction: Reduction {
            rules: BTreeMap::from([(
                Integral(vec![2]),
                BTreeMap::from([(Integral(vec![1]), weight)]),
            )]),
            residuals: vec![Integral(vec![1])],
            ..Default::default()
        },
    };
    PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![2])],
        &[symbol!("target_projection::s")],
        &backend,
        &FlowOptions::default(),
        "regular supplied table",
        &RunContext::default(),
    )
    .unwrap()
}
fn boundary(prepared: &PreparedPhysicalFamily, last: i32, verified: u32) -> CachedBoundary {
    let p = Precision::decimal(100).unwrap();
    CachedBoundary {
        identity: prepared.flow().identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(
            symbol!("target_projection::s"),
            Atom::num(1),
        )])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(-2, last).unwrap(),
        coefficients: (-2..=last)
            .map(|k| vec![p.i(if k < 0 { 0 } else { i64::from(k) + 1 })])
            .collect(),
        accuracy: BoundaryAccuracy::supplied(
            verified,
            p.bits,
            vec![vec![p.real(0)]; (last + 3) as usize],
            "exact polynomial test boundary for the supplied constant differential system",
        )
        .unwrap(),
    }
}
#[test]
fn epsilon_pole_weights_require_extra_master_orders_before_projection() {
    let family = prepared(parse!("1/target_projection::eps^2"));
    let source = boundary(&family, 2, 80);
    let range = EpsilonRange::new(-1, 0).unwrap();
    assert_eq!(
        family
            .required_master_range(&source.point.restart_coordinates().unwrap(), range)
            .unwrap(),
        EpsilonRange::new(-2, 2).unwrap()
    );
    assert!(matches!(
        family.project_targets(&boundary(&family, 1, 80), range, 20),
        Err(Error::InvalidInput(_))
    ));
    let result = family.project_targets(&source, range, 20).unwrap();
    let p = Precision::decimal(100).unwrap();
    assert_eq!(result[0].coefficients[&-1], p.i(2));
    assert_eq!(result[0].coefficients[&0], p.i(3));
    assert_eq!(result[0].verified_digits, 20);
}
#[test]
fn positive_weights_preserve_unknown_series_orders_and_share_master_data() {
    let family = prepared(parse!("target_projection::eps^2"));
    let source = boundary(&family, 0, 80);
    let range = EpsilonRange::new(2, 2).unwrap();
    assert_eq!(
        family
            .required_master_range(&source.point.restart_coordinates().unwrap(), range)
            .unwrap(),
        EpsilonRange::new(-2, 0).unwrap()
    );
    let result = family.project_targets(&source, range, 20).unwrap();
    assert_eq!(
        result[0].coefficients[&2],
        Precision::decimal(100).unwrap().i(1)
    );
    assert!(matches!(
        family.project_targets(&source, EpsilonRange::new(3, 3).unwrap(), 20),
        Err(Error::InvalidInput(_))
    ));
    let other = prepared(parse!("2*target_projection::eps^2"));
    // The master basis/system identity is shared; the retained target map is
    // applied by each prepared object rather than cached as a master value.
    let doubled = other.project_targets(&source, range, 20).unwrap();
    assert_eq!(
        doubled[0].coefficients[&2],
        Precision::decimal(100).unwrap().i(2)
    );
}
#[test]
fn propagated_uncertainty_cannot_be_replaced_by_more_working_bits() {
    let family = prepared(Atom::one());
    let mut source = boundary(&family, 0, 30);
    let p = Precision::decimal(100).unwrap();
    source.coefficients[2][0] = p.i(1);
    assert!(matches!(
        family.project_targets(&source, EpsilonRange::new(0, 0).unwrap(), 40),
        Err(Error::Accuracy(_))
    ));
    let mut wrong = boundary(&family, 0, 80);
    wrong.identity = BoundaryIdentity::new(
        family.flow().system(),
        &[parse!("target_projection::different_master")],
        &Atom::one(),
        Prescription::PlusI0,
        "regular supplied table",
    )
    .unwrap();
    assert!(matches!(
        family.project_targets(&wrong, EpsilonRange::new(0, 0).unwrap(), 20),
        Err(Error::InvalidInput(_))
    ));
    assert!(
        family
            .project_targets(
                &source,
                EpsilonRange {
                    leading: 1,
                    last: 0
                },
                20
            )
            .is_err()
    );
    let mut omitted_poles = boundary(&family, 0, 80);
    omitted_poles.range = EpsilonRange::new(0, 0).unwrap();
    omitted_poles.coefficients.drain(..2);
    omitted_poles.accuracy = BoundaryAccuracy::supplied(
        80,
        p.bits,
        vec![vec![p.real(0)]],
        "finite coefficient alone is insufficient for target projection",
    )
    .unwrap();
    assert!(matches!(
        family.project_targets(&omitted_poles, EpsilonRange::new(0, 0).unwrap(), 20),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn large_cancellation_propagates_master_errors_through_native_series_products() {
    let family = prepared(parse!(
        "10^50*(target_projection::eps^2-1)+target_projection::eps^2"
    ));
    let p = Precision::decimal(100).unwrap();
    let mut source = boundary(&family, 2, 30);
    source.coefficients[4][0] = p.i(1);
    let range = EpsilonRange::new(2, 2).unwrap();
    assert!(matches!(
        family.project_targets(&source, range, 20),
        Err(Error::Accuracy(_))
    ));
    source.accuracy = BoundaryAccuracy::supplied(
        90,
        p.bits,
        vec![vec![p.real(0)]; 5],
        "independent higher-accuracy polynomial boundary",
    )
    .unwrap();
    let result = family.project_targets(&source, range, 20).unwrap();
    assert!(p.close(&result[0].coefficients[&2], &p.i(1), 20));
    assert!(result[0].absolute_errors[&2] < p.tolerance(20));
}

#[test]
fn epsilon_dependent_families_require_a_separate_master_pole_bound() {
    let family = prepared_with_constant(Atom::one(), parse!("-target_projection::eps^10"));
    let source = boundary(&family, 0, 80);
    let range = EpsilonRange::new(-2, 0).unwrap();
    assert!(matches!(
        family.required_master_range(&source.point.restart_coordinates().unwrap(), range),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        family.project_targets(&source, range, 20),
        Err(Error::Unsupported(_))
    ));
}

mod projection_diagnostics {
    use super::*;
    use serde_json::Value;

    const REFUSAL: &str =
        "target reduction amplifies cached master uncertainty beyond the requested accuracy";
    const MARKER: &str = "\nprojection-rejection-v1=";

    fn refusal(result: Result<Vec<ProjectedLaurentExpansion>>) -> Value {
        let message = match result {
            Err(Error::Accuracy(message)) => message,
            other => panic!("expected the public projection accuracy refusal, got {other:?}"),
        };
        let payload = message.strip_prefix(REFUSAL).unwrap();
        assert_eq!(payload.matches(MARKER).count(), 1);
        let diagnostic: Value =
            serde_json::from_str(payload.strip_prefix(MARKER).unwrap()).unwrap();
        assert_eq!(diagnostic["schema"], "projection-rejection-v1");
        assert_eq!(diagnostic["checked_error_is_finite"], true);
        assert_eq!(diagnostic["term_products_are_diagnostic"], true);
        assert!(diagnostic.get("term_detail_error").is_none());
        assert!(
            diagnostic["aggregate_authority"]
                .as_str()
                .unwrap()
                .contains("original fixed-series convolution")
        );
        diagnostic
    }

    fn real_field(object: &Value, field: &str, bits: u32) -> Float {
        Float::parse(object[field].as_str().unwrap(), Some(bits)).unwrap()
    }

    fn term(diagnostic: &Value, master_order: i32) -> &Value {
        let terms = diagnostic["terms"].as_array().unwrap();
        let matches = terms
            .iter()
            .filter(|entry| {
                entry["basis_column_zero_based"] == 0
                    && entry["master_epsilon_order"] == master_order
            })
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1);
        matches[0]
    }

    fn prepared_rows(weights: Vec<Atom>) -> PreparedPhysicalFamily {
        // Reuse the existing public supplied-table fixture's family and chart.
        // These test relations are declared data, not claimed IBP identities.
        let template = prepared(Atom::one());
        let targets = (0..weights.len())
            .map(|index| Integral(vec![i16::try_from(index).unwrap() + 2]))
            .collect::<Vec<_>>();
        let backend = TableBackend {
            name: "ordered diagnostic projection test relations".into(),
            reduction: Reduction {
                rules: targets
                    .iter()
                    .cloned()
                    .zip(
                        weights
                            .into_iter()
                            .map(|weight| BTreeMap::from([(Integral(vec![1]), weight)])),
                    )
                    .collect(),
                residuals: vec![Integral(vec![1])],
                ..Default::default()
            },
        };
        PreparedPhysicalFamily::new(
            template.family(),
            &targets,
            &[symbol!("target_projection::s")],
            &backend,
            &FlowOptions::default(),
            "regular supplied table",
            &RunContext::default(),
        )
        .unwrap()
    }

    fn cancellation_refuses_zero_comparison_errors_and_higher_evidence_succeeds() {
        // The exact coefficient at epsilon^2 is (M+1)*1 - M*1 = 1.
        // Independent zero stored comparisons do not erase the verified-digit
        // floor. More working bits alone do not resolve the inherited error.
        let family = prepared(parse!(
            "10^50*(target_projection::eps^2-1)+target_projection::eps^2"
        ));
        let p = Precision::decimal(100).unwrap();
        let refined = Precision { bits: p.bits + 64 };
        let mut source = boundary(&family, 2, 30);
        source.coefficients[4][0] = p.i(1);
        source.validate().unwrap();
        let before_coefficients = source.coefficients.clone();
        let before_comparisons = source.accuracy.comparison_errors().to_vec();
        let range = EpsilonRange::new(2, 2).unwrap();
        let diagnostic = refusal(family.project_targets(&source, range, 20));
        assert_eq!(source.coefficients, before_coefficients);
        assert_eq!(
            source.accuracy.comparison_errors(),
            before_comparisons.as_slice()
        );
        assert_eq!(source.accuracy.verified_digits(), 30);
        assert_eq!(diagnostic["target_index_zero_based"], 0);
        assert_eq!(diagnostic["epsilon_order"], 2);
        assert_eq!(diagnostic["requested_digits"], 20);
        assert_eq!(diagnostic["inherited_verified_digits"], 30);
        assert_eq!(diagnostic["working_bits"], p.bits);
        assert_eq!(diagnostic["refined_bits"], refined.bits);
        assert_eq!(
            diagnostic["central_low"],
            serde_json::json!([p.real(1).to_string(), p.real(0).to_string()])
        );
        assert_eq!(
            diagnostic["central_refined"],
            serde_json::json!([refined.real(1).to_string(), refined.real(0).to_string()])
        );
        let inherited = real_field(&diagnostic, "propagated_inherited_error", refined.bits);
        let difference = real_field(&diagnostic, "refinement_difference", refined.bits);
        let checked = real_field(&diagnostic, "checked_total_absolute_error", p.bits);
        let allowance = real_field(&diagnostic, "allowance", p.bits);
        assert!(inherited > refined.parse("1e19", "0").unwrap().re);
        assert!(inherited < refined.parse("1e21", "0").unwrap().re);
        assert_eq!(difference, refined.real(0));
        assert!(checked > allowance);
        assert_eq!(allowance, p.tolerance(20));
        for (master_order, weight_order, expected_weight) in [
            (0, 2, "100000000000000000000000000000000000000000000000001"),
            (2, 0, "-100000000000000000000000000000000000000000000000000"),
        ] {
            let entry = term(&diagnostic, master_order);
            assert_eq!(entry["weight_epsilon_order"], weight_order);
            let exact_weight = Atom::parse(
                entry["exact_specialized_weight_coefficient"]
                    .as_str()
                    .unwrap(),
                "projection_diagnostic_test",
                Default::default(),
            )
            .unwrap();
            let expected = refined.parse(expected_weight, "0").unwrap();
            assert_eq!(
                refined.eval(&exact_weight, &Default::default()).unwrap(),
                expected
            );
            assert_eq!(
                real_field(entry, "weight_norm_refined", refined.bits),
                refined.norm(&expected)
            );
            assert_eq!(
                real_field(entry, "stored_comparison_error", refined.bits),
                refined.real(0)
            );
            let floor = real_field(entry, "mandatory_verified_digit_floor", refined.bits);
            assert_eq!(floor, refined.tolerance(30));
            assert_eq!(
                real_field(entry, "effective_master_uncertainty", refined.bits),
                floor
            );
            assert!(
                real_field(entry, "individual_rounded_product", refined.bits)
                    > refined.parse("1e19", "0").unwrap().re
            );
            assert_eq!(
                entry["master"],
                serde_json::json!([p.real(1).to_string(), p.real(0).to_string()])
            );
        }
        // Independent stronger fixture evidence is an existing success control,
        // not permission to upgrade evidence on a production boundary.
        source.accuracy = BoundaryAccuracy::supplied(
            90,
            p.bits,
            vec![vec![p.real(0)]; 5],
            "analytic polynomial fixture with independently stronger supplied evidence",
        )
        .unwrap();
        let result = family.project_targets(&source, range, 20).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coefficients[&2], p.i(1));
        assert!(result[0].absolute_errors[&2] > p.real(0));
        assert!(result[0].absolute_errors[&2] < p.tolerance(20));
        assert_eq!(result[0].verified_digits, 20);
        assert_eq!(result[0].working_bits, p.bits);
    }

    fn first_failed_target_and_order_keep_the_comparison_floor_breakdown() {
        let family = prepared_rows(vec![
            parse!("target_projection::eps^2"),
            parse!("10^50"),
            parse!("10^100"),
        ]);
        let p = Precision::decimal(100).unwrap();
        let refined = Precision { bits: p.bits + 64 };
        let mut source = boundary(&family, 1, 30);
        source.coefficients[3][0] = p.zero();
        let mut comparisons = vec![vec![p.real(0)]; 4];
        comparisons[3][0] = p.tolerance(35);
        source.accuracy = BoundaryAccuracy::supplied(
            30,
            p.bits,
            comparisons,
            "analytic polynomial with nonzero comparison smaller than mandatory floor",
        )
        .unwrap();
        source.validate().unwrap();
        // Every target passes at order zero. The first target also passes both
        // orders. Target 1 fails at order 1; target 2 would fail there as well.
        let earlier = family
            .project_targets(&source, EpsilonRange::new(0, 0).unwrap(), 20)
            .unwrap();
        assert_eq!(earlier.len(), 3);
        assert_eq!(earlier[0].coefficients[&0], p.zero());
        assert_eq!(earlier[1].coefficients[&0], p.parse("1e50", "0").unwrap());
        assert_eq!(earlier[2].coefficients[&0], p.parse("1e100", "0").unwrap());
        let first = prepared(parse!("target_projection::eps^2"))
            .project_targets(&source, EpsilonRange::new(0, 1).unwrap(), 20)
            .unwrap();
        assert_eq!(first[0].coefficients[&0], p.zero());
        assert_eq!(first[0].coefficients[&1], p.zero());
        let diagnostic =
            refusal(family.project_targets(&source, EpsilonRange::new(0, 1).unwrap(), 20));
        assert_eq!(diagnostic["target_index_zero_based"], 1);
        assert_eq!(diagnostic["epsilon_order"], 1);
        let entry = term(&diagnostic, 1);
        assert_eq!(entry["weight_epsilon_order"], 0);
        let comparison = real_field(entry, "stored_comparison_error", refined.bits);
        let floor = real_field(entry, "mandatory_verified_digit_floor", refined.bits);
        assert_eq!(
            entry["stored_comparison_error"].as_str().unwrap(),
            source.accuracy.comparison_errors()[3][0].to_string()
        );
        assert!(comparison > refined.real(0));
        assert!(floor > comparison);
        assert_eq!(floor, refined.tolerance(30));
        assert_eq!(
            real_field(entry, "master_norm_refined", refined.bits),
            refined.real(0)
        );
        assert_eq!(
            real_field(entry, "effective_master_uncertainty", refined.bits),
            floor
        );
        assert!(
            real_field(&diagnostic, "checked_total_absolute_error", p.bits)
                > real_field(&diagnostic, "allowance", p.bits)
        );
    }

    fn earlier_failures_do_not_fabricate_a_target_order_diagnostic() {
        let family = prepared(Atom::one());
        let mut source = boundary(&family, 0, 30);
        let range = EpsilonRange::new(0, 0).unwrap();
        match family.project_targets(&source, range, 40) {
            Err(Error::Accuracy(message)) => {
                assert_eq!(
                    message,
                    "target projection requires sufficient inherited accuracy"
                );
                assert!(!message.contains(MARKER));
            }
            other => panic!("expected inherited-accuracy validation, got {other:?}"),
        }
        let p = Precision::decimal(100).unwrap();
        source.accuracy = BoundaryAccuracy::supplied(
            30,
            p.bits,
            vec![vec![p.tolerance(10)]; 3],
            "deliberately invalid claimed accuracy for the validation regression",
        )
        .unwrap();
        match family.project_targets(&source, range, 20) {
            Err(Error::Accuracy(message)) => {
                assert_eq!(message, "cached comparison error exceeds claimed accuracy");
                assert!(!message.contains(MARKER));
            }
            other => panic!("expected invalid cached comparison evidence, got {other:?}"),
        }
        let pole = prepared(parse!("1/target_projection::eps^2"));
        let truncated = boundary(&pole, 1, 80);
        match pole.project_targets(&truncated, EpsilonRange::new(-1, 0).unwrap(), 20) {
            Err(Error::InvalidInput(message)) => {
                assert!(message.contains("requires master epsilon order 2"));
                assert!(!message.contains(MARKER));
            }
            other => panic!("expected missing master-order refusal, got {other:?}"),
        }
    }

    #[test]
    fn public_projection_rejection_contract() {
        // One libtest invocation, one original OS thread. All fixtures, Atoms,
        // native-capable errors, numeric values and their drops stay in scope.
        cancellation_refuses_zero_comparison_errors_and_higher_evidence_succeeds();
        first_failed_target_and_order_keep_the_comparison_floor_breakdown();
        earlier_failures_do_not_fabricate_a_target_order_diagnostic();
    }
}
