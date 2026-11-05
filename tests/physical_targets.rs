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
fn retained_exact_binding_preserves_long_rational_and_complex_coordinates() {
    let (s, eps) = symbol!("target_projection::s", "target_projection::eps");
    let s_atom = Atom::var(s);
    let epsilon = Atom::var(eps);
    let complex = Atom::num(Complex::new(
        Rational::from((2, 7)),
        Rational::from((3, 11)),
    ));
    let weight = (&s_atom - Atom::one()) / epsilon.pow(2)
        + (s_atom.pow(3) + complex)
            / ((s_atom.pow(2) + Atom::num(2))
                * (Atom::one() - &epsilon)
                * (&s_atom + Atom::num(3)));
    let family = prepared(weight);
    let range = EpsilonRange::new(-4, 1).unwrap();
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let almost_one = parse!("1+1/10^100");
    let points = [
        almost_one.clone(),
        Atom::num(Complex::new(
            Rational::from((5, 7)),
            Rational::from((2, 11)),
        )),
        almost_one + Atom::num(Complex::new(Rational::from(0), Rational::from((2, 11)))),
    ];
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        for coordinate in &points {
            let point = BTreeMap::from([(s, coordinate.clone())]);
            let plan = retained.at_point(&family, &point).unwrap();
            let legacy = family.prepare_target_projection(&point, range).unwrap();
            // Rounding almost_one to either working precision would erase its
            // double pole. Exact binding must retain it and its master demand.
            assert_eq!(
                plan.required_master_range(),
                EpsilonRange::new(-2, 3).unwrap()
            );
            assert_eq!(plan.required_master_range(), legacy.required_master_range());
            let source = complex_projection_boundary(&family, &point, 3, bits, verified);
            let expected = family.project_targets(&source, range, digits).unwrap();
            assert_same_projection(&legacy.project_targets(&source, digits).unwrap(), &expected);
            let found = plan.project_targets(&source, digits).unwrap();
            assert_same_projection(&found, &expected);
            assert_ne!(found[0].coefficients[&-4], Precision { bits }.zero());
            assert_same_projection(
                &retained
                    .at_point(&family, &point)
                    .unwrap()
                    .project_targets(&source, digits)
                    .unwrap(),
                &expected,
            );
            let short = complex_projection_boundary(&family, &point, 2, bits, verified);
            assert_eq!(
                plan.project_targets(&short, digits)
                    .unwrap_err()
                    .to_string(),
                family
                    .project_targets(&short, range, digits)
                    .unwrap_err()
                    .to_string(),
            );
        }
    }
    for invalid in [parse!("1.25`80"), Atom::var(s)] {
        assert!(matches!(
            retained.at_point(&family, &BTreeMap::from([(s, invalid)])),
            Err(Error::InvalidInput(_))
        ));
    }
}

#[test]
fn retained_exact_binding_checks_complex_denominators_before_inversion() {
    let s = Atom::var(symbol!("target_projection::s"));
    let epsilon = Atom::var(symbol!("target_projection::eps"));
    let root = Atom::num(Complex::new(Rational::from((1, 3)), Rational::from((2, 7))));
    let vanishing = (&s - &root).pow(2);
    let point = BTreeMap::from([(symbol!("target_projection::s"), root)]);
    let range = EpsilonRange::new(0, 0).unwrap();
    // The generic denominator leading coefficient vanishes in Q(i). Include a
    // retained pole weight, an omitted weight acquiring a pole, and a removable
    // coefficient singularity whose original weight specializes to 1/2.
    for (weight, generic_last, point_last) in [
        (Atom::one() / (&vanishing + &epsilon), 0, 1),
        (epsilon.pow(3) / (&vanishing + epsilon.pow(4)), -2, 1),
        ((&vanishing + &epsilon) / (&vanishing + 2 * &epsilon), 0, 0),
    ] {
        let family = prepared(weight);
        let retained = family.prepare_kinematic_target_projection(range).unwrap();
        assert_eq!(retained.required_master_range().last, generic_last);
        for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
            let plan = retained.at_point(&family, &point).unwrap();
            assert_eq!(plan.required_master_range().last, point_last);
            let source = complex_projection_boundary(&family, &point, point_last, bits, verified);
            assert_same_projection(
                &plan.project_targets(&source, digits).unwrap(),
                &family.project_targets(&source, range, digits).unwrap(),
            );
        }
    }
    // A true denominator zero must retain the legacy rejection, including when
    // the singularity can only be seen after exact complex specialization.
    let singular = prepared(Atom::one() / vanishing);
    assert!(singular.prepare_target_projection(&point, range).is_err());
    assert!(
        singular
            .prepare_kinematic_target_projection(range)
            .unwrap()
            .at_point(&singular, &point)
            .is_err()
    );
}

#[test]
fn retained_exact_binding_uses_legacy_for_algebraic_coordinates() {
    let family = prepared(parse!("target_projection::s"));
    let range = EpsilonRange::new(0, 0).unwrap();
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let point = BTreeMap::from([(symbol!("target_projection::s"), parse!("2^(1/2)"))]);
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let plan = retained.at_point(&family, &point).unwrap();
        let source = complex_projection_boundary(&family, &point, 0, bits, verified);
        assert_same_projection(
            &plan.project_targets(&source, digits).unwrap(),
            &family.project_targets(&source, range, digits).unwrap(),
        );
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

fn assert_numeric_projection_covered(
    found: &[ProjectedLaurentExpansion],
    exact: &[ProjectedLaurentExpansion],
) {
    assert_eq!(found.len(), exact.len());
    for (found, exact) in found.iter().zip(exact) {
        assert_eq!(found.verified_digits, exact.verified_digits);
        assert_eq!(found.working_bits, exact.working_bits);
        assert_eq!(
            found.coefficients.keys().collect::<Vec<_>>(),
            exact.coefficients.keys().collect::<Vec<_>>()
        );
        let p = Precision {
            bits: found.working_bits + 64,
        };
        for (power, value) in &found.coefficients {
            let difference = p.norm(&p.sub(value, &exact.coefficients[power]));
            let allowance = found.absolute_errors[power].add_round(
                &exact.absolute_errors[power],
                p.bits,
                symbolica::domains::float::RoundingDirection::Up,
            );
            assert!(
                difference <= allowance,
                "epsilon^{power}: difference {difference} exceeds reported errors {allowance}"
            );
        }
    }
}

#[test]
fn numeric_binding_covers_exact_long_rational_and_complex_results_at_80_and_201_bits() {
    let family = prepared(
        parse!(
            "(target_projection::s-1)/target_projection::eps^2+(target_projection::s^3+2/7)/(3*(target_projection::s^2+2)*(1-target_projection::eps))"
        ) + Atom::num(Complex::new(
            Rational::from((2, 7)),
            Rational::from((3, 11)),
        )),
    );
    let range = EpsilonRange::new(-4, 1).unwrap();
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let s = symbol!("target_projection::s");
    let almost_one = parse!("1+1/10^100");
    let points = [
        almost_one.clone(),
        Atom::num(Complex::new(
            Rational::from((5, 7)),
            Rational::from((2, 11)),
        )),
        almost_one + Atom::num(Complex::new(Rational::from(0), Rational::from((2, 11)))),
    ];
    for (bits, verified, digits) in [(80, 18, 16), (201, 48, 32)] {
        for coordinate in &points {
            let point = BTreeMap::from([(s, coordinate.clone())]);
            let plan: PreparedTargetProjection<'_> =
                retained.at_point_numeric(&family, &point).unwrap();
            assert_eq!(
                plan.required_master_range(),
                retained.required_master_range()
            );
            let mut source = complex_projection_boundary(&family, &point, 3, bits, verified);
            let exact = retained
                .at_point(&family, &point)
                .unwrap()
                .project_targets(&source, digits)
                .unwrap();
            let found = plan.project_targets(&source, digits).unwrap();
            assert_numeric_projection_covered(&found, &exact);
            assert_numeric_projection_covered(
                &plan.project_targets(&source, digits).unwrap(),
                &exact,
            );
            source.point = CachedPoint::Derived {
                coordinates: point.clone(),
                working_bits: bits,
                provenance: "same exact numeric-plan coordinates".into(),
            };
            assert_numeric_projection_covered(
                &plan.project_targets(&source, digits).unwrap(),
                &exact,
            );
            // A coordinate which rounds identically at both working precisions
            // is still a different exact point, even if its effect is tiny.
            source.point =
                CachedPoint::Exact(BTreeMap::from([(s, coordinate + parse!("1/10^120"))]));
            assert!(matches!(
                plan.project_targets(&source, digits),
                Err(Error::InvalidInput(_))
            ));
        }
    }
    let point = BTreeMap::from([(s, Atom::num(2))]);
    let plan = retained.at_point_numeric(&family, &point).unwrap();
    drop(retained);
    let source = complex_projection_boundary(&family, &point, 3, 201, 48);
    assert_numeric_projection_covered(
        &plan.project_targets(&source, 32).unwrap(),
        &family.project_targets(&source, range, 32).unwrap(),
    );
}

#[test]
fn numeric_binding_accounts_for_cancellation_lost_at_both_precisions() {
    let family = prepared(parse!("10^100*(target_projection::s-1)"));
    let range = EpsilonRange::new(0, 0).unwrap();
    let point = BTreeMap::from([(symbol!("target_projection::s"), parse!("1+1/10^100"))]);
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let numeric = retained.at_point_numeric(&family, &point).unwrap();
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let source = complex_projection_boundary(&family, &point, 0, bits, verified);
        // The true weight is one, while center-only evaluation can see zero at
        // both p and p+64. Coefficient radii must prevent false source credit.
        assert!(family.project_targets(&source, range, digits).is_ok());
        assert!(matches!(
            numeric.project_targets(&source, digits),
            Err(Error::Accuracy(_))
        ));
    }
    // Resolved cancellation is usable; it need not be rejected categorically.
    let family = prepared(parse!("10^10*(target_projection::s-1)"));
    let point = BTreeMap::from([(symbol!("target_projection::s"), parse!("1+1/10^15"))]);
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let numeric = retained.at_point_numeric(&family, &point).unwrap();
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let source = complex_projection_boundary(&family, &point, 0, bits, verified);
        assert_numeric_projection_covered(
            &numeric.project_targets(&source, digits).unwrap(),
            &family.project_targets(&source, range, digits).unwrap(),
        );
    }
}

#[test]
fn numeric_binding_retains_source_floors_comparison_errors_and_credit() {
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
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let plan = retained.at_point_numeric(&family, &point).unwrap();
    assert!(matches!(
        plan.project_targets(&source, 20),
        Err(Error::Accuracy(_))
    ));
    source.accuracy = BoundaryAccuracy::supplied(
        90,
        p.bits,
        vec![vec![p.real(0)]; 5],
        "higher accuracy cancellation fixture",
    )
    .unwrap();
    assert_numeric_projection_covered(
        &plan.project_targets(&source, 20).unwrap(),
        &family.project_targets(&source, range, 20).unwrap(),
    );
    source.accuracy = BoundaryAccuracy::supplied(
        90,
        p.bits,
        vec![vec![p.real(1)]; 5],
        "explicit nonzero master errors",
    )
    .unwrap();
    assert!(matches!(
        plan.project_targets(&source, 20),
        Err(Error::Accuracy(_))
    ));
    for digits in [0, 91] {
        assert!(matches!(
            plan.project_targets(&source, digits),
            Err(Error::Accuracy(_))
        ));
    }
}

#[test]
fn numeric_binding_uses_exact_fallback_for_short_or_earlier_master_series() {
    let range = EpsilonRange::new(-2, 1).unwrap();
    let family = prepared(parse!(
        "(target_projection::s-1)/target_projection::eps^2+1/(3*target_projection::eps*(1-target_projection::eps))"
    ));
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    for coordinate in [1, 2] {
        let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::num(coordinate))]);
        let numeric = retained.at_point_numeric(&family, &point).unwrap();
        assert_eq!(numeric.required_master_range().last, 3);
        for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
            let source = complex_projection_boundary(&family, &point, 2, bits, verified);
            if coordinate == 1 {
                // Only exact specialization may establish that order 3 is unnecessary.
                assert_same_projection(
                    &numeric.project_targets(&source, digits).unwrap(),
                    &family.project_targets(&source, range, digits).unwrap(),
                );
            } else {
                assert!(matches!(
                    numeric.project_targets(&source, digits),
                    Err(Error::InvalidInput(_))
                ));
            }
        }
    }
    for weight in [
        parse!("1/(1-target_projection::eps)"),
        parse!("target_projection::eps^3"),
        parse!("(target_projection::s-1)*target_projection::eps^3"),
    ] {
        let family = prepared(weight);
        let range = EpsilonRange::new(0, 0).unwrap();
        let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::num(1))]);
        let retained = family.prepare_kinematic_target_projection(range).unwrap();
        let numeric = retained.at_point_numeric(&family, &point).unwrap();
        for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
            let p = Precision { bits };
            let mut source = complex_projection_boundary(&family, &point, 0, bits, verified);
            source.range = EpsilonRange::new(-3, 0).unwrap();
            source.coefficients = vec![
                vec![p.complex(1, 2)],
                vec![p.zero()],
                vec![p.zero()],
                vec![p.zero()],
            ];
            source.accuracy = BoundaryAccuracy::supplied(
                verified,
                bits,
                vec![vec![p.real(0)]; 4],
                "earlier master pole",
            )
            .unwrap();
            assert_same_projection(
                &numeric.project_targets(&source, digits).unwrap(),
                &family.project_targets(&source, range, digits).unwrap(),
            );
        }
    }
}

#[test]
fn numeric_binding_checks_exact_exceptional_and_unresolved_denominators() {
    let (s, eps) = symbol!("target_projection::s", "target_projection::eps");
    let root = Atom::num(Complex::new(Rational::from((1, 3)), Rational::from((2, 7))));
    let vanishing = (Atom::var(s) - &root).pow(2);
    let epsilon = Atom::var(eps);
    let point = BTreeMap::from([(s, root)]);
    let range = EpsilonRange::new(0, 0).unwrap();
    for (weight, point_last) in [
        (Atom::one() / (&vanishing + &epsilon), 1),
        (epsilon.pow(3) / (&vanishing + epsilon.pow(4)), 1),
        ((&vanishing + &epsilon) / (&vanishing + 2 * &epsilon), 0),
    ] {
        let family = prepared(weight);
        let retained = family.prepare_kinematic_target_projection(range).unwrap();
        let numeric = retained.at_point_numeric(&family, &point).unwrap();
        assert_eq!(numeric.required_master_range().last, point_last);
        for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
            let source = complex_projection_boundary(&family, &point, point_last, bits, verified);
            assert_same_projection(
                &numeric.project_targets(&source, digits).unwrap(),
                &family.project_targets(&source, range, digits).unwrap(),
            );
            if point_last > 0 {
                let short =
                    complex_projection_boundary(&family, &point, point_last - 1, bits, verified);
                assert!(matches!(
                    numeric.project_targets(&short, digits),
                    Err(Error::InvalidInput(_))
                ));
            }
        }
    }
    let singular = prepared(Atom::one() / vanishing);
    assert!(
        singular
            .prepare_kinematic_target_projection(range)
            .unwrap()
            .at_point_numeric(&singular, &point)
            .is_err()
    );

    // Guard vanishing requires exact Q(i), with i^2=-1.
    let family = prepared(parse!(
        "1/(target_projection::s^2+1+target_projection::eps)"
    ));
    let point = BTreeMap::from([(
        s,
        Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
    )]);
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let numeric = retained.at_point_numeric(&family, &point).unwrap();
    assert_eq!(numeric.required_master_range().last, 1);
    let source = complex_projection_boundary(&family, &point, 1, 201, 48);
    assert_same_projection(
        &numeric.project_targets(&source, 32).unwrap(),
        &family.project_targets(&source, range, 32).unwrap(),
    );

    // Nonzero exact guard, but its numerical enclosure contains zero at both
    // working precisions. Projection falls back instead of using an inverse.
    let family = prepared(parse!("1/(target_projection::s-1)"));
    let point = BTreeMap::from([(s, parse!("1+1/10^100"))]);
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let numeric = retained.at_point_numeric(&family, &point).unwrap();
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let source = complex_projection_boundary(&family, &point, 0, bits, verified);
        assert_same_projection(
            &numeric.project_targets(&source, digits).unwrap(),
            &family.project_targets(&source, range, digits).unwrap(),
        );
    }
}

#[test]
fn numeric_binding_preserves_family_point_basis_and_domain_guards() {
    let range = EpsilonRange::new(0, 0).unwrap();
    let family = prepared(parse!("target_projection::s"));
    let retained = family.prepare_kinematic_target_projection(range).unwrap();
    let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::num(1))]);
    for other in [
        prepared(Atom::num(2)),
        prepared_with_constant(Atom::one(), Atom::num(-2)),
    ] {
        assert!(matches!(
            retained.at_point_numeric(&other, &point),
            Err(Error::InvalidInput(_))
        ));
    }
    for invalid in [
        BTreeMap::new(),
        BTreeMap::from([(symbol!("target_projection::s"), parse!("1.25`80"))]),
    ] {
        assert!(matches!(
            retained.at_point_numeric(&family, &invalid),
            Err(Error::InvalidInput(_))
        ));
    }
    let plan = retained.at_point_numeric(&family, &point).unwrap();
    let mut source = boundary(&family, 0, 80);
    let other = prepared_with_constant(Atom::one(), Atom::num(-2));
    source.identity = other.flow().identity().clone();
    assert!(matches!(
        plan.project_targets(&source, 20),
        Err(Error::InvalidInput(_))
    ));
    let mut source = boundary(&family, 0, 80);
    source.range = EpsilonRange::new(0, 0).unwrap();
    source.coefficients.drain(..2);
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    source.accuracy =
        BoundaryAccuracy::supplied(80, p.bits, vec![vec![p.real(0)]], "missing ordinary poles")
            .unwrap();
    assert!(matches!(
        plan.project_targets(&source, 20),
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
            .at_point_numeric(&domain, &zero),
        Err(Error::InvalidInput(_))
    ));
    // Algebraic coordinates preserve the legacy evaluator and error semantics.
    let point = BTreeMap::from([(symbol!("target_projection::s"), parse!("2^(1/2)"))]);
    let plan = retained.at_point_numeric(&family, &point).unwrap();
    let source = complex_projection_boundary(&family, &point, 0, 201, 48);
    assert_same_projection(
        &plan.project_targets(&source, 32).unwrap(),
        &family.project_targets(&source, range, 32).unwrap(),
    );
}

#[test]
fn numeric_projection_uses_owner_epsilon_shearing_and_master_column_order() {
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
        let plan = kinematic.at_point_numeric(&prepared, &point).unwrap();
        let source = complex_projection_boundary(&prepared, &point, 1, 201, 48);
        assert_numeric_projection_covered(
            &plan.project_targets(&source, 32).unwrap(),
            &prepared.project_targets(&source, range, 32).unwrap(),
        );
        // Original B=epsilon^-1 J_B; t3 is the unshifted second column.
        let projected = plan.project_targets(&source, 32).unwrap();
        assert_eq!(projected[0].coefficients[&0], source.coefficients[3][0]);
        assert_eq!(projected[1].coefficients[&0], source.coefficients[2][1]);
    }
}

fn combination_family(
    first: BTreeMap<Integral, Atom>,
    second: BTreeMap<Integral, Atom>,
) -> PreparedPhysicalFamily {
    // A supplied constant system isolates projection algebra, not integration.
    let family = IntegralFamily {
        name: "shared_master_combination_fixture".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-1)]],
        propagators: [-1, -2]
            .into_iter()
            .enumerate()
            .map(|(slot, constant)| Propagator {
                constant: Atom::num(constant),
                scalar_products: vec![Atom::one(), Atom::num(2 * slot as i64)],
            })
            .collect(),
        physical_propagators: 2,
        epsilon: symbol!("target_projection::eps"),
        dimension: 4,
    };
    let targets = [Integral(vec![2, 0]), Integral(vec![3, 0])];
    let backend = TableBackend {
        name: "exact shared master test relations".into(),
        reduction: Reduction {
            rules: BTreeMap::from([(targets[0].clone(), first), (targets[1].clone(), second)]),
            residuals: vec![Integral(vec![0, 1]), Integral(vec![1, 0])],
            ..Default::default()
        },
    };
    PreparedPhysicalFamily::new(
        &family,
        &targets,
        &[symbol!("target_projection::s")],
        &backend,
        &FlowOptions::default(),
        "regular combination test table",
        &RunContext::default(),
    )
    .unwrap()
}

#[test]
fn target_combinations_cancel_shared_master_before_errors_at_80_and_201_bits() {
    let (j1, j2) = (Integral(vec![0, 1]), Integral(vec![1, 0]));
    let family = combination_family(
        BTreeMap::from([(j1.clone(), parse!("10^20")), (j2.clone(), Atom::one())]),
        BTreeMap::from([(j1.clone(), parse!("-10^20"))]),
    );
    assert_eq!(family.basis(), &[j1, j2]);
    let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::one())]);
    let range = EpsilonRange::new(-2, 0).unwrap();
    let weights = [vec![Atom::one(), Atom::one()]];
    let plan: PreparedTargetProjection<'_> = family
        .prepare_target_combinations(&point, &weights, range)
        .unwrap();
    assert_eq!(plan.required_master_range(), range);
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let p = Precision { bits };
        let mut source = complex_projection_boundary(&family, &point, 0, bits, verified);
        let reordered = family
            .prepare_target_combinations(
                &point,
                &[
                    vec![Atom::new(), Atom::one()],
                    vec![Atom::one(), Atom::new()],
                ],
                range,
            )
            .unwrap();
        let individual = family.project_targets(&source, range, digits).unwrap();
        let mut reverse = individual.clone();
        reverse.reverse();
        assert_same_projection(
            &reordered.project_targets(&source, digits).unwrap(),
            &reverse,
        );
        let found = plan.project_targets(&source, digits).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].working_bits, bits);
        assert_eq!(found[0].verified_digits, digits);
        for (index, power) in (-2..=0).enumerate() {
            assert_eq!(found[0].coefficients[&power], source.coefficients[index][1]);
            assert!(found[0].absolute_errors[&power] > p.real(0));
            assert!(
                found[0].absolute_errors[&power]
                    < individual[0].absolute_errors[&power].clone()
                        + &individual[1].absolute_errors[&power]
            );
        }
        // Compare to the surviving J2 identity projection, including its
        // inherited source floor. Cancellation removes J1's dependence only.
        let reference = prepared(Atom::one());
        let mut reference_source = source.clone();
        reference_source.identity = reference.flow().identity().clone();
        reference_source.coefficients = source
            .coefficients
            .iter()
            .map(|row| vec![row[1].clone()])
            .collect();
        reference_source.accuracy = BoundaryAccuracy::supplied(
            verified,
            bits,
            vec![vec![p.real(0)]; 3],
            "same surviving J2 coefficients",
        )
        .unwrap();
        assert_same_projection(
            &found,
            &reference
                .project_targets(&reference_source, range, digits)
                .unwrap(),
        );
        // An inconsistent claimed boundary grade is rejected even when its
        // affected column cancels: the input contract is not weakened.
        source.accuracy = BoundaryAccuracy::supplied(
            verified,
            bits,
            vec![vec![p.real(1), p.real(0)]; 3],
            "invalid comparison error claim",
        )
        .unwrap();
        assert!(matches!(
            plan.project_targets(&source, digits),
            Err(Error::Accuracy(_))
        ));
    }
}

#[test]
fn target_combinations_bind_exact_gaussian_weights_and_cancel_poles_before_depth() {
    let (j1, j2) = (Integral(vec![0, 1]), Integral(vec![1, 0]));
    let family = combination_family(
        BTreeMap::from([
            (j1.clone(), parse!("10^20/target_projection::eps^5")),
            (j2, Atom::one()),
        ]),
        BTreeMap::from([(j1, parse!("-10^20/target_projection::eps^5"))]),
    );
    let (s, eps) = symbol!("target_projection::s", "target_projection::eps");
    let point = BTreeMap::from([(s, parse!("1+1/10^100"))]);
    let range = EpsilonRange::new(-3, 0).unwrap();
    let gaussian = Atom::num(Complex::new(
        Rational::from((2, 7)),
        Rational::from((3, 11)),
    ));
    let weight: Atom = (&gaussian * (Atom::var(s) - Atom::one()) * parse!("10^100"))
        / (3 * Atom::var(eps) * (Atom::one() - Atom::var(eps)));
    let plan = family
        .prepare_target_combinations(&point, &[vec![weight.clone(), weight]], range)
        .unwrap();
    // Both fifth-order poles cancel; only the weight's simple pole remains.
    assert_eq!(
        plan.required_master_range(),
        EpsilonRange::new(-2, 1).unwrap()
    );
    assert_eq!(family.required_master_range(&point, range).unwrap().last, 5);
    let reference =
        prepared(gaussian / parse!("3*target_projection::eps*(1-target_projection::eps)"));
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let source = complex_projection_boundary(&family, &point, 1, bits, verified);
        let mut reference_source = source.clone();
        reference_source.identity = reference.flow().identity().clone();
        reference_source.coefficients = source
            .coefficients
            .iter()
            .map(|row| vec![row[1].clone()])
            .collect();
        let p = Precision { bits };
        reference_source.accuracy = BoundaryAccuracy::supplied(
            verified,
            bits,
            vec![vec![p.real(0)]; 4],
            "same surviving J2 coefficients",
        )
        .unwrap();
        assert_same_projection(
            &plan.project_targets(&source, digits).unwrap(),
            &reference
                .project_targets(&reference_source, range, digits)
                .unwrap(),
        );
        let short = complex_projection_boundary(&family, &point, 0, bits, verified);
        assert!(matches!(
            plan.project_targets(&short, digits),
            Err(Error::InvalidInput(_))
        ));
    }
    // At s=1 the exact coefficient vanishes; no rounded-zero inference is used.
    let zero = BTreeMap::from([(s, Atom::one())]);
    let plan = family
        .prepare_target_combinations(
            &zero,
            &[vec![Atom::var(s) - Atom::one(), Atom::var(s) - Atom::one()]],
            range,
        )
        .unwrap();
    assert_eq!(
        plan.required_master_range(),
        EpsilonRange::new(-2, -2).unwrap()
    );
    let source = complex_projection_boundary(&family, &zero, -2, 80, 16);
    let found = plan.project_targets(&source, 10).unwrap();
    assert!(
        found[0]
            .coefficients
            .values()
            .all(|v| *v == Precision { bits: 80 }.zero())
    );
    assert!(found[0].absolute_errors.values().all(|e| e.is_zero()));
}

#[test]
fn target_combinations_validate_matrix_point_identity_depth_and_source_credit() {
    let family = prepared(parse!("1/target_projection::eps^2"));
    let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::one())]);
    let range = EpsilonRange::new(-1, 0).unwrap();
    for weights in [
        vec![],
        vec![vec![]],
        vec![vec![Atom::one(), Atom::one()]],
        vec![vec![Atom::one()], vec![]],
    ] {
        assert!(matches!(
            family.prepare_target_combinations(&point, &weights, range),
            Err(Error::InvalidInput(_))
        ));
    }
    for invalid in [
        BTreeMap::new(),
        BTreeMap::from([(family.family().epsilon, Atom::one())]),
    ] {
        assert!(matches!(
            family.prepare_target_combinations(&invalid, &[vec![Atom::one()]], range),
            Err(Error::InvalidInput(_))
        ));
    }
    assert!(matches!(
        family.prepare_target_combinations(
            &point,
            &[vec![Atom::one()]],
            EpsilonRange {
                leading: 1,
                last: 0
            },
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        family.prepare_target_combinations(
            &point,
            &[vec![parse!("sin(target_projection::eps)")]],
            range,
        ),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        family.prepare_target_combinations(
            &point,
            &[vec![parse!("target_projection::undeclared")]],
            range,
        ),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        family.prepare_target_combinations(
            &point,
            &[vec![parse!("1/(target_projection::s-1)")]],
            range,
        ),
        Err(Error::InvalidInput(_))
    ));
    let plan = family
        .prepare_target_combinations(&point, &[vec![Atom::one()]], range)
        .unwrap();
    let source = boundary(&family, 2, 80);
    assert_same_projection(
        &plan.project_targets(&source, 20).unwrap(),
        &family.project_targets(&source, range, 20).unwrap(),
    );
    assert!(matches!(
        plan.project_targets(&boundary(&family, 1, 80), 20),
        Err(Error::InvalidInput(_))
    ));
    let mut wrong_point = source.clone();
    wrong_point.point = CachedPoint::Exact(BTreeMap::from([(
        symbol!("target_projection::s"),
        Atom::num(2),
    )]));
    assert!(
        matches!(plan.project_targets(&wrong_point, 20), Err(Error::InvalidInput(message)) if message.contains("different physical point coordinates"))
    );
    let other = prepared_with_constant(Atom::one(), Atom::num(-2));
    let mut wrong_family = source.clone();
    wrong_family.identity = other.flow().identity().clone();
    assert!(
        matches!(plan.project_targets(&wrong_family, 20), Err(Error::InvalidInput(message)) if message.contains("different physical basis identity"))
    );
    for digits in [0, 81] {
        assert!(matches!(
            plan.project_targets(&source, digits),
            Err(Error::Accuracy(_))
        ));
    }
    let mut omitted = source.clone();
    omitted.range = EpsilonRange::new(0, 2).unwrap();
    omitted.coefficients.drain(..2);
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    omitted.accuracy =
        BoundaryAccuracy::supplied(80, p.bits, vec![vec![p.real(0)]; 3], "omitted poles").unwrap();
    assert!(matches!(
        plan.project_targets(&omitted, 20),
        Err(Error::InvalidInput(_))
    ));
    let epsilon_family = prepared_with_constant(Atom::one(), parse!("-target_projection::eps^10"));
    assert!(matches!(
        epsilon_family.prepare_target_combinations(&point, &[vec![Atom::one()]], range),
        Err(Error::Unsupported(_))
    ));
    let domain = prepared(parse!("1/target_projection::s"));
    let zero = BTreeMap::from([(symbol!("target_projection::s"), Atom::new())]);
    assert!(matches!(
        domain.prepare_target_combinations(&zero, &[vec![Atom::one()]], range),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn target_combinations_preserve_earlier_poles_and_uncancelled_source_floors() {
    let point = BTreeMap::from([(symbol!("target_projection::s"), Atom::one())]);
    for weight in [
        parse!("1/(1-target_projection::eps)"),
        parse!("target_projection::eps^3"),
    ] {
        let family = prepared(Atom::one());
        let reference = prepared(weight.clone());
        let range = EpsilonRange::new(0, 0).unwrap();
        let plan = family
            .prepare_target_combinations(&point, &[vec![weight]], range)
            .unwrap();
        let mut source = boundary(&family, 0, 80);
        let p = Precision {
            bits: source.accuracy.working_bits(),
        };
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
        let found = plan.project_targets(&source, 20).unwrap();
        assert_eq!(found[0].coefficients[&0], p.complex(1, 2));
        assert_same_projection(
            &found,
            &reference.project_targets(&source, range, 20).unwrap(),
        );
    }
    let (j1, j2) = (Integral(vec![0, 1]), Integral(vec![1, 0]));
    let family = combination_family(
        BTreeMap::from([(j1, parse!("10^20"))]),
        BTreeMap::from([(j2, parse!("-10^20"))]),
    );
    let range = EpsilonRange::new(0, 0).unwrap();
    let plan = family
        .prepare_target_combinations(&point, &[vec![Atom::one(), Atom::one()]], range)
        .unwrap();
    for (bits, verified, digits, refined_verified) in [(80, 16, 10, 40), (201, 48, 32, 60)] {
        let p = Precision { bits };
        let mut source = complex_projection_boundary(&family, &point, 0, bits, verified);
        source.coefficients = vec![vec![p.i(1), p.i(1)]; 3];
        // Equal central values of different masters do not imply covariance.
        // Zero comparison errors still leave the inherited source floors.
        assert!(matches!(
            plan.project_targets(&source, digits),
            Err(Error::Accuracy(_))
        ));
        let refined = Precision { bits: bits + 128 };
        source.coefficients = vec![vec![refined.i(1), refined.i(1)]; 3];
        source.accuracy = BoundaryAccuracy::supplied(
            refined_verified,
            refined.bits,
            vec![vec![refined.real(0); 2]; 3],
            "independently refined source",
        )
        .unwrap();
        let found = plan.project_targets(&source, digits).unwrap();
        assert_eq!(found[0].coefficients[&0], refined.zero());
        assert!(found[0].absolute_errors[&0] > refined.real(0));
    }
}

#[test]
fn target_combinations_use_cached_column_order_and_epsilon_shearing() {
    // Same supplied evanescent bubble as the existing shearing regressions.
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
    let context = RunContext::default();
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[bubble.clone(), t3.clone(), Integral(vec![2, 1])],
        &[s],
        &backend,
        &FlowOptions::default(),
        "Euclidean projection fixture",
        &context,
    )
    .unwrap()
    .with_epsilon_shearing(&context)
    .unwrap();
    assert_eq!(prepared.basis(), &[bubble, t3]);
    assert_eq!(prepared.epsilon_shearing().unwrap().weights(), &[-1, 0]);
    let weights = [
        vec![epsilon.clone(), Atom::new(), Atom::new()],
        vec![Atom::new(), Atom::one(), Atom::new()],
        vec![
            (Atom::one() - 2 * &epsilon) / &denominator,
            2 * &mass / (&epsilon * &denominator),
            Atom::new(),
        ],
    ];
    let range = EpsilonRange::new(-3, 0).unwrap();
    for coordinate in [1, 2] {
        let point = BTreeMap::from([(s, Atom::num(coordinate))]);
        let plan = prepared
            .prepare_target_combinations(&point, &weights, range)
            .unwrap();
        assert_eq!(
            plan.required_master_range(),
            EpsilonRange::new(-2, 1).unwrap()
        );
        for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
            let source = complex_projection_boundary(&prepared, &point, 1, bits, verified);
            let found = plan.project_targets(&source, digits).unwrap();
            let legacy = prepared.project_targets(&source, range, digits).unwrap();
            assert_eq!(found.len(), 3);
            assert_same_projection(&found[1..], &legacy[1..]);
            for (index, power) in (-2..=0).enumerate() {
                // epsilon * original bubble is the cached first master.
                assert_eq!(found[0].coefficients[&power], source.coefficients[index][0]);
            }
            // The mixed target's natural bound is -3 here. No physical
            // cancellation or stronger pole bound is inferred from the range.
            assert_ne!(found[2].coefficients[&-3], Precision { bits }.zero());
        }
    }
}
