use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};
use symbolica::prelude::*;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::physical_transport::BoundaryAttemptOutcome;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{Error, FlowOptions, Precision, Prescription, Result, RunContext, RustFlow};

fn coordinate() -> Symbol {
    symbol!("accuracy_fallback::s")
}
fn point(value: Atom) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(coordinate(), value)])
}
fn weak_point() -> Atom {
    Atom::num(Rational::from((19, 10)))
}
fn flow() -> Result<RustFlow> {
    RustFlow::new(
        KinematicSystem {
            epsilon: symbol!("accuracy_fallback::eps"),
            derivatives: BTreeMap::from([(
                coordinate(),
                vec![vec![Atom::num(Complex::new(
                    Rational::from(0),
                    Rational::from(1),
                ))]],
            )]),
        },
        &[Atom::var(symbol!("accuracy_fallback::Y"))],
        &Atom::one(),
        Prescription::PlusI0,
        "real s; exp(i*s)",
    )
}
fn boundary(flow: &RustFlow, coordinate: Atom, digits: u32) -> Result<CachedBoundary> {
    let p = Precision::decimal(80)?;
    let s = p.eval(&coordinate, &Default::default())?;
    let value = p.exp(&p.mul(&s, &p.parse("0", "1")?));
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(point(coordinate)),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 0)?,
        coefficients: vec![vec![value]],
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.real(0)]],
            "independent analytic exp(i*s); conservative supplied evidence",
        )?,
    })
}
fn bank(flow: &RustFlow, far_digits: u32) -> Result<RustFlowCache> {
    let mut bank = RustFlowCache::default();
    bank.insert_many(vec![
        boundary(flow, weak_point(), 20)?,
        boundary(flow, Atom::new(), far_digits)?,
    ])?;
    Ok(bank)
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 30,
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
fn assert_unchanged(bank: &RustFlowCache, old: &RustFlowCache) -> Result<()> {
    assert_eq!(bank.len(), old.len());
    for (a, b) in bank.entries().iter().zip(old.entries()) {
        assert_eq!(a.identity.key(), b.identity.key());
        assert_eq!(
            a.point.restart_coordinates()?,
            b.point.restart_coordinates()?
        );
        assert_eq!(a.coefficients, b.coefficients);
        assert_eq!(a.accuracy.verified_digits(), b.accuracy.verified_digits());
        assert_eq!(
            a.accuracy.comparison_errors(),
            b.accuracy.comparison_errors()
        );
    }
    Ok(())
}

#[test]
fn farther_accurate_source_succeeds_after_nearby_uncertainty_failure() -> Result<()> {
    let flow = flow()?;
    let mut bank = bank(&flow, 60)?;
    let range = EpsilonRange::new(0, 0)?;
    let result = flow.evaluate_to(
        &mut bank,
        &point(Atom::num(2)),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert_eq!(
        result.starting_point.restart_coordinates()?,
        point(Atom::new())
    );
    assert_eq!(result.boundary_attempts.len(), 2);
    let attempts = &result.boundary_attempts;
    assert_eq!(
        attempts[0].starting_point.restart_coordinates()?,
        point(weak_point())
    );
    assert!(matches!(
        attempts[0].outcome,
        BoundaryAttemptOutcome::AccuracyRejected { .. }
    ));
    assert_eq!(attempts[0].source_verified_digits, 20);
    assert_eq!(attempts[1].source_verified_digits, 60);
    assert!(matches!(
        attempts[1].outcome,
        BoundaryAttemptOutcome::Accepted
    ));
    assert!(attempts[0].cost < attempts[1].cost);
    let p = Precision::decimal(80)?;
    let expected = p.exp(&p.parse("0", "2")?);
    assert!(p.close(&result.boundary.coefficients[0][0], &expected, 20));
    assert!(bank.len() > 2);
    let hit = flow.evaluate_to(
        &mut bank,
        &point(Atom::num(2)),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.boundary_attempts.len(), 1);
    assert!(matches!(
        hit.boundary_attempts[0].outcome,
        BoundaryAttemptOutcome::Accepted
    ));
    Ok(())
}

#[test]
fn budget_and_exhausted_accuracy_leave_the_original_bank_unchanged() -> Result<()> {
    let flow = flow()?;
    let mut bank = bank(&flow, 60)?;
    let old = bank.clone();
    let limited = FlowOptions {
        max_boundary_attempts: 1,
        ..options()
    };
    assert!(
        matches!(flow.evaluate_to(&mut bank,&point(Atom::num(2)),EpsilonRange::new(0,0)?,&limited,&RunContext::default(),&policy()),Err(Error::Limit(message)) if message.contains("budget 1 exhausted"))
    );
    assert_unchanged(&bank, &old)?;
    let mut weak = self::bank(&flow, 20)?;
    let old = weak.clone();
    assert!(
        matches!(flow.evaluate_to(&mut weak,&point(Atom::num(2)),EpsilonRange::new(0,0)?,&options(),&RunContext::default(),&policy()),Err(Error::Accuracy(message)) if message.contains("none of 2 compatible"))
    );
    assert_unchanged(&weak, &old)?;
    let invalid = FlowOptions {
        max_boundary_attempts: 0,
        ..options()
    };
    assert!(matches!(invalid.validate(), Err(Error::InvalidInput(_))));
    Ok(())
}

#[test]
fn cancellation_between_source_attempts_does_not_commit_failed_trajectories() -> Result<()> {
    let flow = flow()?;
    let mut bank = bank(&flow, 60)?;
    let old = bank.clone();
    let context = RunContext::default();
    let token = context.cancellation.clone();
    let seen = Arc::new(Mutex::new(Vec::new()));
    let recorded = seen.clone();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: move |source: &CachedBoundary, _: &CachedPoint| {
            let coordinate = source.point.restart_coordinates()?[&coordinate()].clone();
            recorded.lock().unwrap().push(coordinate.clone());
            if coordinate.is_zero() {
                token.cancel();
            }
            Ok(true)
        },
    };
    assert!(matches!(
        flow.evaluate_to(
            &mut bank,
            &point(Atom::num(2)),
            EpsilonRange::new(0, 0)?,
            &options(),
            &context,
            &policy
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(*seen.lock().unwrap(), vec![weak_point(), Atom::new()]);
    assert_unchanged(&bank, &old)?;
    Ok(())
}
