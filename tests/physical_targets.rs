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
    assert!(matches!(
        family.prepare_target_projection(&source.point.restart_coordinates().unwrap(), range),
        Err(Error::Unsupported(_))
    ));
}

fn assert_same_projection(
    found: &[ProjectedLaurentExpansion],
    expected: &[ProjectedLaurentExpansion],
) {
    assert_eq!(found.len(), expected.len());
    for (found, expected) in found.iter().zip(expected) {
        assert_eq!(found.coefficients, expected.coefficients);
        assert_eq!(found.absolute_errors, expected.absolute_errors);
        assert_eq!(found.verified_digits, expected.verified_digits);
        assert_eq!(found.working_bits, expected.working_bits);
    }
}

#[test]
fn prepared_projection_matches_existing_apis_at_80_and_201_bits() {
    // Substitution at s=1 removes the double pole; s=2 retains it. The
    // remaining rational series exercises nontrivial depth and p/p+64 rounding.
    let family = prepared(parse!(
        "(target_projection::s-1)/target_projection::eps^2+1/(3*target_projection::eps*(1-target_projection::eps))"
    ));
    let s = symbol!("target_projection::s");
    let range = EpsilonRange::new(-2, 1).unwrap();
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let p = Precision { bits };
        for coordinate in [1, 2] {
            let point = BTreeMap::from([(s, Atom::num(coordinate))]);
            let plan: PreparedTargetProjection<'_> =
                family.prepare_target_projection(&point, range).unwrap();
            let required = plan.required_master_range();
            assert_eq!(
                required,
                family.required_master_range(&point, range).unwrap()
            );
            assert_eq!(required, EpsilonRange::new(-2, coordinate + 1).unwrap());
            let mut source = boundary(&family, required.last, verified);
            source.point = CachedPoint::Exact(point.clone());
            source.coefficients = (-2..=required.last)
                .map(|k| vec![p.i(if k < 0 { 0 } else { i64::from(k) + 1 })])
                .collect();
            source.accuracy = BoundaryAccuracy::supplied(
                verified,
                bits,
                vec![vec![p.real(0)]; source.coefficients.len()],
                "supplied polynomial master coefficients at the requested working bits",
            )
            .unwrap();
            let expected = family.project_targets(&source, range, digits).unwrap();
            let found = plan.project_targets(&source, digits).unwrap();
            assert_same_projection(&found, &expected);
            assert_eq!(found[0].working_bits, bits);
            assert_eq!(found[0].coefficients[&-2], p.i(i64::from(coordinate - 1)));
            // A prepared plan can be reused without changing its exact rows.
            assert_same_projection(&plan.project_targets(&source, digits).unwrap(), &expected);
            // Coordinate origin is not upgraded by preparation; an identical
            // retained path image can use the plan, with its own accuracy cap.
            source.point = CachedPoint::Derived {
                coordinates: point,
                working_bits: bits,
                provenance: "same exact path image for projection regression".into(),
            };
            assert_same_projection(&plan.project_targets(&source, digits).unwrap(), &expected);
        }
    }
}

#[test]
fn prepared_projection_rejects_mismatches_and_preserves_depth_and_accuracy_guards() {
    let family = prepared(parse!("1/target_projection::eps^2"));
    let source = boundary(&family, 2, 80);
    let point = source.point.restart_coordinates().unwrap();
    let range = EpsilonRange::new(-1, 0).unwrap();
    let plan = family.prepare_target_projection(&point, range).unwrap();
    assert_eq!(
        plan.required_master_range(),
        EpsilonRange::new(-2, 2).unwrap()
    );

    let mut wrong_point = source.clone();
    wrong_point.point = CachedPoint::Exact(BTreeMap::from([(
        symbol!("target_projection::s"),
        Atom::num(2),
    )]));
    wrong_point.validate().unwrap();
    assert!(family.project_targets(&wrong_point, range, 20).is_ok());
    assert!(matches!(
        plan.project_targets(&wrong_point, 20),
        Err(Error::InvalidInput(message)) if message.contains("different physical point coordinates")
    ));
    let mut wrong_identity = source.clone();
    let other_family = prepared_with_constant(Atom::one(), Atom::num(-2));
    wrong_identity.identity = other_family.flow().identity().clone();
    wrong_identity.validate().unwrap();
    assert!(matches!(
        plan.project_targets(&wrong_identity, 20),
        Err(Error::InvalidInput(message)) if message.contains("different physical basis identity")
    ));

    let short = boundary(&family, 1, 80);
    assert!(matches!(
        plan.project_targets(&short, 20),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(
        plan.project_targets(&short, 20).unwrap_err().to_string(),
        family
            .project_targets(&short, range, 20)
            .unwrap_err()
            .to_string()
    );
    let mut omitted_poles = source.clone();
    omitted_poles.range = EpsilonRange::new(0, 2).unwrap();
    omitted_poles.coefficients.drain(..2);
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    omitted_poles.accuracy = BoundaryAccuracy::supplied(
        80,
        p.bits,
        vec![vec![p.real(0)]; 3],
        "deliberately omitted master poles",
    )
    .unwrap();
    assert!(matches!(
        plan.project_targets(&omitted_poles, 20),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(
        plan.project_targets(&omitted_poles, 20)
            .unwrap_err()
            .to_string(),
        family
            .project_targets(&omitted_poles, range, 20)
            .unwrap_err()
            .to_string()
    );
    for digits in [0, 81] {
        assert!(matches!(
            plan.project_targets(&source, digits),
            Err(Error::Accuracy(_))
        ));
        assert_eq!(
            plan.project_targets(&source, digits)
                .unwrap_err()
                .to_string(),
            family
                .project_targets(&source, range, digits)
                .unwrap_err()
                .to_string()
        );
    }
    assert!(matches!(
        family.prepare_target_projection(&BTreeMap::new(), range),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        family.prepare_target_projection(
            &point,
            EpsilonRange {
                leading: 1,
                last: 0
            }
        ),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn prepared_projection_preserves_master_orders_below_the_ordinary_bound() {
    // The second weight is beyond the prepared absolute truncation and only
    // contributes when a boundary supplies a nonzero earlier master pole.
    for weight in [
        parse!("1/(1-target_projection::eps)"),
        parse!("target_projection::eps^3"),
    ] {
        let family = prepared(weight);
        let mut source = boundary(&family, 0, 80);
        let point = source.point.restart_coordinates().unwrap();
        let range = EpsilonRange::new(0, 0).unwrap();
        let plan = family.prepare_target_projection(&point, range).unwrap();
        source.range = EpsilonRange::new(-3, 0).unwrap();
        let p = Precision {
            bits: source.accuracy.working_bits(),
        };
        source.coefficients = vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]];
        source.accuracy = BoundaryAccuracy::supplied(
            80,
            p.bits,
            vec![vec![p.real(0)]; 4],
            "nonzero master coefficient below the ordinary pole bound",
        )
        .unwrap();
        let expected = family.project_targets(&source, range, 20).unwrap();
        let found = plan.project_targets(&source, 20).unwrap();
        assert_same_projection(&found, &expected);
        assert_eq!(found[0].coefficients[&0], p.i(1));
    }
}

#[test]
fn prepared_projection_retains_the_source_floor_and_cancellation_guard() {
    let family = prepared(parse!(
        "10^50*(target_projection::eps^2-1)+target_projection::eps^2"
    ));
    let mut source = boundary(&family, 2, 30);
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    source.coefficients[4][0] = p.i(1);
    let point = source.point.restart_coordinates().unwrap();
    let range = EpsilonRange::new(2, 2).unwrap();
    let plan = family.prepare_target_projection(&point, range).unwrap();
    // Comparison errors are zero, so only inherited source floors reveal this
    // excessive cancellation. Working precision alone cannot raise source credit.
    assert!(matches!(
        plan.project_targets(&source, 20),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(
        plan.project_targets(&source, 20).unwrap_err().to_string(),
        family
            .project_targets(&source, range, 20)
            .unwrap_err()
            .to_string()
    );
    source.accuracy = BoundaryAccuracy::supplied(
        90,
        p.bits,
        vec![vec![p.real(0)]; 5],
        "independently higher-accuracy polynomial test input",
    )
    .unwrap();
    assert_same_projection(
        &plan.project_targets(&source, 20).unwrap(),
        &family.project_targets(&source, range, 20).unwrap(),
    );
}

fn complex_projection_boundary(
    family: &PreparedPhysicalFamily,
    point: &BTreeMap<Symbol, Atom>,
    last: i32,
    bits: u32,
    verified: u32,
) -> CachedBoundary {
    let p = Precision { bits };
    let mut source = boundary(family, last, verified);
    source.point = CachedPoint::Exact(point.clone());
    source.coefficients = (-2..=last)
        .map(|k| {
            (0..family.basis().len())
                .map(|column| p.complex(i64::from(k) + 3 + column as i64, 2 * i64::from(k) + 1))
                .collect()
        })
        .collect();
    source.accuracy = BoundaryAccuracy::supplied(
        verified,
        bits,
        vec![vec![p.real(0); family.basis().len()]; source.coefficients.len()],
        "exact complex polynomial master fixture",
    )
    .unwrap();
    source
}

#[test]
fn kinematic_projection_matches_exact_complex_values_and_errors_at_80_and_201_bits() {
    let family = prepared(parse!(
        "(target_projection::s-1)/target_projection::eps^2+1/(3*target_projection::eps*(1-target_projection::eps))"
    ));
    let range = EpsilonRange::new(-4, 1).unwrap();
    // This public type annotation also exercises the exported API and lifetime.
    let kinematic: PreparedKinematicTargetProjection =
        family.prepare_kinematic_target_projection(range).unwrap();
    assert_eq!(
        kinematic.required_master_range(),
        EpsilonRange::new(-2, 3).unwrap()
    );
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let p = Precision { bits };
        for coordinate in [1, 2] {
            let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::num(coordinate))]);
            let plan: PreparedTargetProjection<'_> = kinematic.at_point(&family, &point).unwrap();
            let legacy = family.prepare_target_projection(&point, range).unwrap();
            let required = plan.required_master_range();
            assert_eq!(required, legacy.required_master_range());
            assert_eq!(required.last, coordinate + 1);
            let mut source =
                complex_projection_boundary(&family, &point, required.last, bits, verified);
            let expected = family.project_targets(&source, range, digits).unwrap();
            assert_same_projection(&legacy.project_targets(&source, digits).unwrap(), &expected);
            let found = plan.project_targets(&source, digits).unwrap();
            assert_same_projection(&found, &expected);
            assert_eq!(
                found[0].coefficients[&-4],
                p.complex(i64::from(coordinate - 1), -3 * i64::from(coordinate - 1))
            );
            // Reuse does not bind again or promote derived coordinates.
            source.point = CachedPoint::Derived {
                coordinates: point,
                working_bits: bits,
                provenance: "same exact coordinates for retained kinematic projection".into(),
            };
            assert_same_projection(&plan.project_targets(&source, digits).unwrap(), &expected);
        }
    }
}

#[test]
fn kinematic_projection_is_owned_and_rejects_other_family_target_maps() {
    struct Session {
        family: PreparedPhysicalFamily,
        projection: PreparedKinematicTargetProjection,
    }
    let range = EpsilonRange::new(0, 0).unwrap();
    let family = prepared(Atom::one());
    let projection = family.prepare_kinematic_target_projection(range).unwrap();
    // Both values can move into one session, with no borrowed sibling field.
    let session = Session { family, projection };
    let source = boundary(&session.family, 0, 80);
    let point = source.point.restart_coordinates().unwrap();
    let plan = session
        .projection
        .at_point(&session.family, &point)
        .unwrap();
    // The point plan's lifetime is tied only to the supplied family.
    drop(session.projection);
    assert_same_projection(
        &plan.project_targets(&source, 20).unwrap(),
        &session.family.project_targets(&source, range, 20).unwrap(),
    );
    let projection = session
        .family
        .prepare_kinematic_target_projection(range)
        .unwrap();
    let other_targets = prepared(Atom::num(2));
    assert_eq!(
        other_targets.flow().identity().key(),
        session.family.flow().identity().key()
    );
    assert!(matches!(
        projection.at_point(&other_targets, &point),
        Err(Error::InvalidInput(_))
    ));
    let other_basis = prepared_with_constant(Atom::one(), Atom::num(-2));
    assert!(matches!(
        projection.at_point(&other_basis, &point),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn kinematic_projection_specializes_exceptional_denominators_before_using_generic_depth() {
    let range = EpsilonRange::new(0, 0).unwrap();
    let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::num(0))]);
    for (weight, generic_last, point_last) in [
        ("1/(target_projection::s+target_projection::eps)", 0, 1),
        // Initially outside the retained series: checking only retained rows
        // would silently miss its specialized epsilon pole.
        (
            "target_projection::eps^3/(target_projection::s+target_projection::eps^4)",
            -2,
            1,
        ),
        (
            "(target_projection::s+target_projection::eps)/(target_projection::s+target_projection::eps^2)",
            0,
            1,
        ),
        // Coefficients of the generic series are singular but the specialized
        // exact rational weight is the regular constant 1/2.
        (
            "(target_projection::s+target_projection::eps)/(target_projection::s+2*target_projection::eps)",
            0,
            0,
        ),
    ] {
        let family = prepared(Atom::parse(weight, "", Default::default()).unwrap());
        let kinematic = family.prepare_kinematic_target_projection(range).unwrap();
        assert_eq!(kinematic.required_master_range().last, generic_last);
        let plan = kinematic.at_point(&family, &point).unwrap();
        assert_eq!(plan.required_master_range().last, point_last);
        assert_eq!(
            plan.required_master_range(),
            family.required_master_range(&point, range).unwrap()
        );
        let source = complex_projection_boundary(&family, &point, point_last, 201, 48);
        assert_same_projection(
            &plan.project_targets(&source, 32).unwrap(),
            &family.project_targets(&source, range, 32).unwrap(),
        );
        if point_last > 0 {
            let short = complex_projection_boundary(&family, &point, point_last - 1, 201, 48);
            assert!(matches!(
                plan.project_targets(&short, 32),
                Err(Error::InvalidInput(_))
            ));
        }
    }
    // A genuine singularity is not interpreted as a vanishing weight.
    let family = prepared(parse!("1/target_projection::s"));
    let kinematic = family.prepare_kinematic_target_projection(range).unwrap();
    assert!(kinematic.at_point(&family, &point).is_err());
    assert!(family.prepare_target_projection(&point, range).is_err());

    // The denominator guard must specialize in Q(i), including cancellations
    // that would be invisible if the imaginary unit stayed a free parameter.
    let family = prepared(parse!(
        "1/(target_projection::s^2+1+target_projection::eps)"
    ));
    let point = BTreeMap::from([(
        symbol!("target_projection::s"),
        Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
    )]);
    let kinematic = family.prepare_kinematic_target_projection(range).unwrap();
    assert_eq!(kinematic.required_master_range().last, 0);
    let plan = kinematic.at_point(&family, &point).unwrap();
    assert_eq!(plan.required_master_range().last, 1);
    let source = complex_projection_boundary(&family, &point, 1, 201, 48);
    assert_same_projection(
        &plan.project_targets(&source, 32).unwrap(),
        &family.project_targets(&source, range, 32).unwrap(),
    );
}

#[test]
fn kinematic_projection_preserves_point_identity_depth_accuracy_and_domain_guards() {
    let range = EpsilonRange::new(-1, 0).unwrap();
    let family = prepared(parse!("1/target_projection::eps^2"));
    let mut source = boundary(&family, 2, 80);
    let point = source.point.restart_coordinates().unwrap();
    let kinematic = family.prepare_kinematic_target_projection(range).unwrap();
    let plan = kinematic.at_point(&family, &point).unwrap();
    let mut wrong_point = source.clone();
    wrong_point.point = CachedPoint::Exact(BTreeMap::from([(
        symbol!("target_projection::s"),
        Atom::num(2),
    )]));
    assert!(
        matches!(plan.project_targets(&wrong_point, 20), Err(Error::InvalidInput(message)) if message.contains("different physical point"))
    );
    let other = prepared_with_constant(Atom::one(), Atom::num(-2));
    let mut wrong_identity = source.clone();
    wrong_identity.identity = other.flow().identity().clone();
    for (candidate, digits) in [
        (boundary(&family, 1, 80), 20),
        (wrong_identity, 20),
        (source.clone(), 0),
        (source.clone(), 81),
    ] {
        assert_eq!(
            plan.project_targets(&candidate, digits)
                .unwrap_err()
                .to_string(),
            family
                .project_targets(&candidate, range, digits)
                .unwrap_err()
                .to_string(),
        );
    }
    // Missing ordinary pole coefficients retain the old error, too.
    source.range = EpsilonRange::new(0, 2).unwrap();
    source.coefficients.drain(..2);
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    source.accuracy =
        BoundaryAccuracy::supplied(80, p.bits, vec![vec![p.real(0)]; 3], "omitted poles").unwrap();
    assert_eq!(
        plan.project_targets(&source, 20).unwrap_err().to_string(),
        family
            .project_targets(&source, range, 20)
            .unwrap_err()
            .to_string(),
    );
    assert!(matches!(
        kinematic.at_point(&family, &BTreeMap::new()),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        family.prepare_kinematic_target_projection(EpsilonRange {
            leading: 1,
            last: 0
        }),
        Err(Error::InvalidInput(_))
    ));

    let domain = prepared_with_constant(
        parse!("1/target_projection::s"),
        parse!("-target_projection::s"),
    );
    let zero = BTreeMap::from([(symbol!("target_projection::s"), Atom::num(0))]);
    assert!(matches!(
        domain
            .prepare_kinematic_target_projection(range)
            .unwrap()
            .at_point(&domain, &zero),
        Err(Error::InvalidInput(_))
    ));
    let epsilon_dependent =
        prepared_with_constant(Atom::one(), parse!("-target_projection::eps^10"));
    assert!(matches!(
        epsilon_dependent.prepare_kinematic_target_projection(range),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn kinematic_projection_keeps_source_floors_and_earlier_master_poles() {
    let range = EpsilonRange::new(2, 2).unwrap();
    let family = prepared(parse!(
        "10^50*(target_projection::eps^2-1)+target_projection::eps^2"
    ));
    let mut source = boundary(&family, 2, 30);
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    source.coefficients[4][0] = p.i(1);
    let point = source.point.restart_coordinates().unwrap();
    let kinematic = family.prepare_kinematic_target_projection(range).unwrap();
    let plan = kinematic.at_point(&family, &point).unwrap();
    assert!(matches!(
        plan.project_targets(&source, 20),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(
        plan.project_targets(&source, 20).unwrap_err().to_string(),
        family
            .project_targets(&source, range, 20)
            .unwrap_err()
            .to_string(),
    );
    source.accuracy = BoundaryAccuracy::supplied(
        90,
        p.bits,
        vec![vec![p.real(0)]; 5],
        "refined polynomial source",
    )
    .unwrap();
    assert_same_projection(
        &plan.project_targets(&source, 20).unwrap(),
        &family.project_targets(&source, range, 20).unwrap(),
    );

    for weight in [
        parse!("1/(1-target_projection::eps)"),
        parse!("target_projection::eps^3"),
        parse!("(target_projection::s-1)*target_projection::eps^3"),
    ] {
        let family = prepared(weight);
        let range = EpsilonRange::new(0, 0).unwrap();
        let kinematic = family.prepare_kinematic_target_projection(range).unwrap();
        let plan = kinematic.at_point(&family, &point).unwrap();
        let mut source = boundary(&family, 0, 80);
        source.range = EpsilonRange::new(-3, 0).unwrap();
        source.coefficients = vec![
            vec![p.complex(1, 2)],
            vec![p.zero()],
            vec![p.zero()],
            vec![p.zero()],
        ];
        source.accuracy =
            BoundaryAccuracy::supplied(80, p.bits, vec![vec![p.real(0)]; 4], "earlier master pole")
                .unwrap();
        assert_same_projection(
            &plan.project_targets(&source, 20).unwrap(),
            &family.project_targets(&source, range, 20).unwrap(),
        );
    }
}

#[test]
fn kinematic_projection_uses_owner_epsilon_shearing_and_master_column_order() {
    // The same two-master bubble relations used by epsilon_physical_family,
    // with a supplied boundary: no AMF calculation is needed for projection.
    let (s, eps) = symbol!("target_projection::s", "target_projection::eps");
    let mass = Atom::var(s);
    let epsilon = Atom::var(eps);
    let family = IntegralFamily {
        name: "projection_evanescent_bubble".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-1)]],
        propagators: vec![
            Propagator {
                constant: -&mass,
                scalar_products: vec![Atom::one(), Atom::new()],
            },
            Propagator {
                constant: -&mass - Atom::one(),
                scalar_products: vec![Atom::one(), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: eps,
        dimension: 4,
    };
    let (bubble, t3) = (Integral(vec![1, 1]), Integral(vec![3, 0]));
    let denominator = 4 * &mass + Atom::one();
    let mixed = BTreeMap::from([
        (bubble.clone(), (Atom::one() - 2 * &epsilon) / &denominator),
        (t3.clone(), 2 * &mass / (&epsilon * &denominator)),
    ]);
    let backend = TableBackend {
        name: "supplied bubble shearing fixture".into(),
        reduction: Reduction {
            rules: BTreeMap::from([
                (Integral(vec![2, 1]), mixed.clone()),
                (Integral(vec![1, 2]), mixed),
                (
                    Integral(vec![4, 0]),
                    BTreeMap::from([(t3.clone(), -(Atom::one() + &epsilon) / (3 * &mass))]),
                ),
            ]),
            residuals: vec![bubble.clone(), t3.clone()],
            ..Default::default()
        },
    };
    let targets = [bubble, t3, Integral(vec![2, 1])];
    let context = RunContext::default();
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &targets,
        &[s],
        &backend,
        &FlowOptions::default(),
        "Euclidean projection fixture",
        &context,
    )
    .unwrap()
    .with_epsilon_shearing(&context)
    .unwrap();
    assert_eq!(prepared.epsilon_shearing().unwrap().weights(), &[-1, 0]);
    let range = EpsilonRange::new(-2, 0).unwrap();
    let kinematic = prepared.prepare_kinematic_target_projection(range).unwrap();
    assert_eq!(
        kinematic.required_master_range(),
        EpsilonRange::new(-2, 1).unwrap()
    );
    for coordinate in [1, 2] {
        let point = BTreeMap::from([(s, Atom::num(coordinate))]);
        let plan = kinematic.at_point(&prepared, &point).unwrap();
        let source = complex_projection_boundary(&prepared, &point, 1, 201, 48);
        assert_same_projection(
            &plan.project_targets(&source, 32).unwrap(),
            &prepared.project_targets(&source, range, 32).unwrap(),
        );
        // Original B=epsilon^-1 J_B; t3 is the unshifted second column.
        let projected = plan.project_targets(&source, 32).unwrap();
        assert_eq!(projected[0].coefficients[&0], source.coefficients[3][0]);
        assert_eq!(projected[1].coefficients[&0], source.coefficients[2][1]);
    }
}
