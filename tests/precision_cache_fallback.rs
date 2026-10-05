use std::cell::Cell;
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::physical_transport::BoundaryAttemptOutcome;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{
    CancellationToken, Error, FlowOptions, Precision, Prescription, Result, RunContext, RustFlow,
};

const SUPPLIED_DIGITS: u32 = 180;

fn coordinate() -> Symbol {
    symbol!("precision_cache_fallback::s")
}
fn large() -> Rational {
    Integer::from(2).pow(400).into()
}
fn slope() -> Rational {
    Rational::from((1, 3)) - large()
}
fn near() -> Rational {
    Rational::one() - Rational::one() / large()
}
fn exact_value(s: &Rational) -> Rational {
    large() + slope() * s
}
fn point(s: Rational) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(coordinate(), Atom::num(s))])
}
fn flow() -> Result<RustFlow> {
    RustFlow::new(
        KinematicSystem {
            epsilon: symbol!("precision_cache_fallback::eps"),
            derivatives: BTreeMap::from([(
                coordinate(),
                vec![
                    vec![Atom::zero(), Atom::num(slope())],
                    vec![Atom::zero(), Atom::zero()],
                ],
            )]),
        },
        &[
            Atom::var(symbol!("precision_cache_fallback::Y1")),
            Atom::var(symbol!("precision_cache_fallback::Y2")),
        ],
        &Atom::one(),
        Prescription::PlusI0,
        "entire nilpotent connection; Y(s)=[2^400+(1/3-2^400)*s,1]",
    )
}
fn boundary(flow: &RustFlow, s: Rational) -> Result<CachedBoundary> {
    let p = Precision::decimal(220)?;
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(point(s.clone())),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 0)?,
        coefficients: vec![vec![p.rational(&exact_value(&s)), p.i(1)]],
        accuracy: BoundaryAccuracy::supplied(
            SUPPLIED_DIGITS,
            p.bits,
            vec![vec![p.real(0); 2]],
            "exact rational analytic solution rounded directly at 220 decimal digits; conservative 180-digit supplied evidence",
        )?,
    })
}
fn bank(flow: &RustFlow) -> Result<RustFlowCache> {
    let mut bank = RustFlowCache::default();
    bank.insert_many(vec![
        boundary(flow, Rational::zero())?,
        boundary(flow, near())?,
    ])?;
    Ok(bank)
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 40,
        series_order: 16,
        max_steps: 100,
        max_precision_attempts: 1,
        ..Default::default()
    }
}

// Deliberately rank the ill-conditioned route first. Compare coordinates
// exactly: 1-2^-400 would round to the target at the initial working precision.
#[derive(Default)]
struct ExhaustingSourceFirst {
    near_prices: Cell<usize>,
    cancel_before_second_attempt: Option<CancellationToken>,
}
impl TransportCost for ExhaustingSourceFirst {
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        let coordinates = source.point.restart_coordinates()?;
        if coordinates == point(near()) {
            let count = self.near_prices.get() + 1;
            self.near_prices.set(count);
            if count == 2
                && let Some(token) = &self.cancel_before_second_attempt
            {
                token.cancel();
            }
        }
        Ok(Some(p.real(
            if coordinates == target.restart_coordinates()?
                || coordinates == point(Rational::zero())
            {
                0
            } else {
                1
            },
        )))
    }
}
fn assert_boundary_unchanged(a: &CachedBoundary, b: &CachedBoundary) -> Result<()> {
    assert_eq!(a.identity.key(), b.identity.key());
    assert_eq!(
        a.point.restart_coordinates()?,
        b.point.restart_coordinates()?
    );
    assert_eq!(a.coefficients, b.coefficients);
    assert_eq!(a.accuracy.verified_digits(), b.accuracy.verified_digits());
    assert_eq!(a.accuracy.working_bits(), b.accuracy.working_bits());
    assert_eq!(
        a.accuracy.comparison_errors(),
        b.accuracy.comparison_errors()
    );
    Ok(())
}
fn assert_unchanged(bank: &RustFlowCache, old: &RustFlowCache) -> Result<()> {
    assert_eq!(bank.len(), old.len());
    for (a, b) in bank.entries().iter().zip(old.entries()) {
        assert_boundary_unchanged(a, b)?;
    }
    Ok(())
}

#[test]
fn exact_gaussian_solution_validates_both_boundaries_and_target() {
    let s = Atom::var(coordinate());
    let solution = Atom::num(large()) + Atom::num(slope()) * &s;
    assert!((solution.derivative(coordinate()) - Atom::num(slope())).is_zero());
    assert_eq!(exact_value(&Rational::zero()), large());
    assert_eq!(exact_value(&Rational::one()), Rational::from((1, 3)));
    assert_eq!(
        exact_value(&near()),
        Rational::from((4, 3)) - Rational::one() / (Rational::from(3) * large())
    );
    // The near route has O(1) normalized coefficients, without changing the
    // connection, target, or independently supplied boundary uncertainty.
    assert_eq!(
        slope() * (Rational::one() - near()),
        -Rational::one() + Rational::one() / (Rational::from(3) * large())
    );
}

#[test]
fn exhausted_precision_source_falls_back_to_a_conditioned_cached_route() -> Result<()> {
    let flow = flow()?;
    let mut bank = bank(&flow)?;
    let original = bank.clone();
    let range = EpsilonRange::new(0, 0)?;
    let policy = ExhaustingSourceFirst::default();
    let result = flow.evaluate_to(
        &mut bank,
        &point(Rational::one()),
        range,
        &options(),
        &RunContext::default(),
        &policy,
    )?;
    let attempts = &result.boundary_attempts;
    assert_eq!(attempts.len(), 2);
    assert_eq!(
        attempts[0].starting_point.restart_coordinates()?,
        point(Rational::zero())
    );
    assert!(matches!(
        &attempts[0].outcome,
        BoundaryAttemptOutcome::AccuracyRejected { message }
            if message.contains("insufficient working precision")
                && message.contains("precision retry budget exhausted")
    ));
    assert_eq!(
        attempts[1].starting_point.restart_coordinates()?,
        point(near())
    );
    assert!(matches!(
        attempts[1].outcome,
        BoundaryAttemptOutcome::Accepted
    ));
    assert!(attempts[0].cost < attempts[1].cost);
    assert_eq!(attempts[0].source_verified_digits, SUPPLIED_DIGITS);
    assert_eq!(attempts[1].source_verified_digits, SUPPLIED_DIGITS);
    assert_eq!(result.starting_point.restart_coordinates()?, point(near()));
    assert!(result.boundary.accuracy.verified_digits() >= options().digits);
    assert!(result.boundary.accuracy.verified_digits() <= SUPPLIED_DIGITS);
    let p = Precision::decimal(220)?;
    assert!(p.close(
        &result.boundary.coefficients[0][0],
        &p.rational(&Rational::from((1, 3))),
        20
    ));
    assert!(p.close(&result.boundary.coefficients[0][1], &p.i(1), 20));
    for source in original.entries() {
        let unchanged = bank
            .entries()
            .iter()
            .find(|entry| {
                entry.point.restart_coordinates().ok() == source.point.restart_coordinates().ok()
            })
            .unwrap();
        assert_boundary_unchanged(unchanged, source)?;
    }
    assert!(bank.len() > original.len());
    let hit = flow.evaluate_to(
        &mut bank,
        &point(Rational::one()),
        range,
        &options(),
        &RunContext::default(),
        &policy,
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.boundary_attempts.len(), 1);
    Ok(())
}

#[test]
fn precision_fallback_respects_attempt_cap_and_preserves_failed_cache() -> Result<()> {
    let flow = flow()?;
    let mut bank = bank(&flow)?;
    let old = bank.clone();
    let limited = FlowOptions {
        max_boundary_attempts: 1,
        ..options()
    };
    let error = flow
        .evaluate_to(
            &mut bank,
            &point(Rational::one()),
            EpsilonRange::new(0, 0)?,
            &limited,
            &RunContext::default(),
            &ExhaustingSourceFirst::default(),
        )
        .err()
        .expect("transport must fail");
    assert!(
        matches!(error, Error::Limit(message) if message.contains("budget 1 exhausted") && message.contains("insufficient working precision"))
    );
    assert_unchanged(&bank, &old)?;
    let mut only_exhausting = RustFlowCache::default();
    only_exhausting.insert(boundary(&flow, Rational::zero())?)?;
    let old = only_exhausting.clone();
    let error = flow
        .evaluate_to(
            &mut only_exhausting,
            &point(Rational::one()),
            EpsilonRange::new(0, 0)?,
            &options(),
            &RunContext::default(),
            &ExhaustingSourceFirst::default(),
        )
        .err()
        .expect("transport must fail");
    assert!(
        matches!(error, Error::Accuracy(message) if message.contains("none of 1 compatible") && message.contains("insufficient working precision"))
    );
    assert_unchanged(&only_exhausting, &old)?;
    Ok(())
}

#[test]
fn cancellation_after_precision_exhaustion_does_not_commit_a_trajectory() -> Result<()> {
    let flow = flow()?;
    let mut bank = bank(&flow)?;
    let old = bank.clone();
    let context = RunContext::default();
    let policy = ExhaustingSourceFirst {
        near_prices: Cell::new(0),
        cancel_before_second_attempt: Some(context.cancellation.clone()),
    };
    let error = flow
        .evaluate_to(
            &mut bank,
            &point(Rational::one()),
            EpsilonRange::new(0, 0)?,
            &options(),
            &context,
            &policy,
        )
        .err()
        .expect("transport must fail");
    assert!(matches!(error, Error::Cancelled));
    assert_eq!(policy.near_prices.get(), 2);
    assert_unchanged(&bank, &old)?;
    Ok(())
}
