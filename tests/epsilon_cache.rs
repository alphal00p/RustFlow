use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{
    contour::PolynomialPrescription, kinematics::KinematicSystem, transport_cache::*, *,
};

fn names() -> (Symbol, Symbol) {
    (symbol!("epsilon_cache::s"), symbol!("epsilon_cache::eps"))
}
fn atom(s: &str) -> Atom {
    Atom::parse(s, "epsilon_cache", Default::default()).unwrap()
}
fn system(entry: Atom) -> KinematicSystem {
    let (s, epsilon) = names();
    KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(
            s,
            vec![
                vec![Atom::new(), entry],
                vec![Atom::var(epsilon), Atom::new()],
            ],
        )]),
    }
}
fn flow(entry: Atom) -> Result<RustFlow> {
    RustFlow::new(
        system(entry),
        &[atom("I(1)"), atom("I(2)")],
        &Atom::num(3),
        Prescription::PlusI0,
        "analytic hyperbolic solution on a regular physical domain",
    )
}
fn point(n: i64) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(names().0, Atom::num(n))])
}
fn source(flow: &RustFlow, at: i64, last: i32) -> Result<CachedBoundary> {
    let p = Precision::decimal(80)?;
    let range = EpsilonRange::new(-1, last)?;
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(point(at)),
        kind: PointKind::Physical,
        range,
        coefficients: (-1..=last)
            .map(|power| vec![p.zero(), p.i(i64::from(power == 0 || power == 1))])
            .collect(),
        accuracy: BoundaryAccuracy::supplied(
            40,
            p.bits,
            (-1..=last)
                .map(|power| vec![p.tolerance(56), p.tolerance((55 - power) as u32)])
                .collect(),
            "Declared global pole bound minus one; exact source I=(0,1+epsilon), with conservative coefficient errors",
        )?,
    })
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 40,
        series_order: 64,
        ..Default::default()
    }
}
fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}

#[test]
fn sheared_cache_transports_pole_and_finite_terms_with_extra_source_order() -> Result<()> {
    let original = flow(atom("1/eps"))?;
    let regular = original.regularize_epsilon(&RunContext::default())?;
    assert_eq!(regular.shearing().weights(), &[-1, 0]);
    let requested = EpsilonRange::new(-1, 0)?;
    let needed = regular.required_range(requested)?;
    assert_eq!(needed, EpsilonRange::new(-1, 1)?);
    assert!(matches!(
        regular.to_sheared_boundary(&source(&original, 0, 0)?, needed.last),
        Err(Error::InvalidInput(_))
    ));
    let input = source(&original, 0, 1)?;
    let converted = regular.to_sheared_boundary(&input, needed.last)?;
    assert_eq!(
        converted.accuracy.verified_digits(),
        input.accuracy.verified_digits()
    );
    assert_eq!(
        converted.accuracy.working_bits(),
        input.accuracy.working_bits()
    );
    assert_eq!(
        converted.accuracy.comparison_errors()[1][0],
        input.accuracy.comparison_errors()[0][0]
    );
    assert_eq!(
        converted.accuracy.comparison_errors()[2][1],
        input.accuracy.comparison_errors()[2][1]
    );
    let p = Precision::decimal(80)?;
    assert_eq!(converted.accuracy.comparison_errors()[0][0], p.real(0));
    let mut bank = RustFlowCache::default();
    bank.insert(input)?;
    bank.insert(converted)?;
    let result = regular.flow().evaluate_to(
        &mut bank,
        &point(1),
        needed,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(result.transport.is_some());
    assert!(result.inserted_points > 0);
    let restored = regular.to_original_boundary(&result.boundary, requested)?;
    let exponential = p.exp(&p.i(1));
    let inverse = p.exp(&p.i(-1));
    let sinh = p.div(&p.sub(&exponential, &inverse), &p.i(2));
    let cosh = p.div(&p.add(&exponential, &inverse), &p.i(2));
    for (values, expected) in restored
        .coefficients
        .iter()
        .zip([vec![sinh.clone(), p.zero()], vec![sinh, cosh]])
    {
        for (value, expected) in values.iter().zip(expected) {
            assert!(p.close(value, &expected, 20), "{value} != {expected}");
        }
    }
    assert_eq!(restored.identity.key(), original.identity().key());
    assert_eq!(
        restored.accuracy.verified_digits(),
        result.boundary.accuracy.verified_digits()
    );
    let directory = std::env::temp_dir().join(format!("epsilon-cache-{}", std::process::id()));
    bank.save(&directory)?;
    let mut loaded = RustFlowCache::load(&directory)?;
    let hit = regular.flow().evaluate_to(
        &mut loaded,
        &point(1),
        needed,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.boundary.identity.key(), regular.flow().identity().key());
    std::fs::remove_dir_all(directory)?;
    Ok(())
}

#[test]
fn epsilon_conversion_authenticates_basis_and_retains_original_domains() -> Result<()> {
    let original = flow(atom("(s^2-1)/(eps*(s-1)) + 1/(s-eps)"))?;
    let regular = original.regularize_epsilon(&RunContext::default())?;
    for at in [0, 1] {
        assert!(
            regular
                .to_sheared_boundary(&source(&original, at, 1)?, 1)
                .is_err()
        );
    }
    let valid = source(&original, 2, 1)?;
    let converted = regular.to_sheared_boundary(&valid, 1)?;
    let mut bank = RustFlowCache::default();
    bank.insert(converted.clone())?;
    assert!(regular.to_sheared_boundary(&converted, 1).is_err());
    assert!(
        regular
            .to_original_boundary(&valid, EpsilonRange::new(-1, 0)?)
            .is_err()
    );
    let unscaled_labels = RustFlow::new(
        regular.flow().system().clone(),
        &[atom("I(1)"), atom("I(2)")],
        &Atom::num(3),
        Prescription::PlusI0,
        "analytic hyperbolic solution on a regular physical domain",
    )?;
    assert_ne!(
        unscaled_labels.identity().key(),
        regular.flow().identity().key()
    );
    let before = bank.entries().len();
    assert!(
        regular
            .flow()
            .evaluate_to(
                &mut bank,
                &point(0),
                converted.range,
                &options(),
                &RunContext::default(),
                &policy(),
            )
            .is_err()
    );
    assert_eq!(before, bank.entries().len());
    Ok(())
}

#[test]
fn epsilon_rescaling_preserves_prescription_and_homotopy_identity() -> Result<()> {
    let continuation = |side| PhysicalContinuation {
        prescriptions: vec![PolynomialPrescription {
            polynomial: atom("s-10"),
            prescription: side,
        }],
        unprescribed_side: side,
        domain: "specified physical continuation without extra windings".into(),
    };
    let a =
        flow(atom("1/eps"))?.with_prescribed_continuation(continuation(Prescription::PlusI0))?;
    let b =
        flow(atom("1/eps"))?.with_prescribed_continuation(continuation(Prescription::MinusI0))?;
    let regular = a.regularize_epsilon(&RunContext::default())?;
    let other = b.regularize_epsilon(&RunContext::default())?;
    assert_ne!(
        regular.flow().identity().key(),
        other.flow().identity().key()
    );
    let seed = regular.to_sheared_boundary(&source(&a, 0, 1)?, 1)?;
    assert!(other.to_sheared_boundary(&source(&a, 0, 1)?, 1).is_err());
    let mut bank = RustFlowCache::default();
    bank.insert(seed.clone())?;
    // Prescribed identities still require the explicit prescribed entrypoint.
    assert!(
        regular
            .flow()
            .evaluate_to(
                &mut bank,
                &point(1),
                seed.range,
                &options(),
                &RunContext::default(),
                &policy(),
            )
            .is_err()
    );
    let result = regular.flow().evaluate_prescribed_to(
        &mut bank,
        &point(1),
        seed.range,
        &options(),
        &RunContext::default(),
        &policy(),
        &|_: &CachedBoundary, _: &CachedPoint, _: &physical_transport::PhysicalRoute| Ok(true),
    )?;
    assert!(result.transport.is_some());
    regular
        .to_original_boundary(&result.boundary, EpsilonRange::new(-1, 0)?)?
        .validate()?;
    Ok(())
}
