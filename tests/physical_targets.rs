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
