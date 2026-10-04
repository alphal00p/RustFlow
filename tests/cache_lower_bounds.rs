use std::{cell::RefCell, collections::BTreeMap};
use symbolica::prelude::*;
use symbolica_amflow::{kinematics::KinematicSystem, transport_cache::*, *};

fn fixture() -> (Symbol, BoundaryIdentity, Precision) {
    let (s, epsilon) = symbol!("cache_bounds::s", "cache_bounds::epsilon");
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]])]),
    };
    let identity = BoundaryIdentity::new(
        &system,
        &[parse!("cache_bounds::one")],
        &Atom::one(),
        Prescription::PlusI0,
        "constant analytic solution",
    )
    .unwrap();
    (s, identity, Precision::decimal(60).unwrap())
}
fn entry(identity: &BoundaryIdentity, s: Symbol, coordinate: Atom, digits: u32) -> CachedBoundary {
    let p = Precision::decimal(90).unwrap();
    CachedBoundary {
        identity: identity.clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, coordinate)])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 0).unwrap(),
        coefficients: vec![vec![p.i(1)]],
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.real(0)]],
            "Exact analytic constant one, not an estimate from working precision",
        )
        .unwrap(),
    }
}
fn index(boundary: &CachedBoundary) -> i64 {
    let value = boundary
        .point
        .restart_coordinates()
        .unwrap()
        .into_values()
        .next()
        .unwrap();
    (1..=8).find(|&i| value == Atom::num(i)).unwrap()
}
struct Costs<B, V> {
    bounds: B,
    values: V,
    calls: RefCell<Vec<i64>>,
}
impl<B, V> TransportCost for Costs<B, V>
where
    B: Fn(i64, Precision) -> Result<Option<Float>>,
    V: Fn(i64, Precision) -> Result<Option<Float>>,
{
    fn lower_bound(
        &self,
        source: &CachedBoundary,
        _: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        (self.bounds)(index(source), p)
    }
    fn cost(
        &self,
        source: &CachedBoundary,
        _: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        self.calls.borrow_mut().push(index(source));
        (self.values)(index(source), p)
    }
}
fn costs<B, V>(bounds: B, values: V) -> Costs<B, V> {
    Costs {
        bounds,
        values,
        calls: RefCell::new(Vec::new()),
    }
}
fn selected(
    cache: &RustFlowCache,
    identity: &BoundaryIdentity,
    s: Symbol,
    p: Precision,
    policy: &dyn TransportCost,
) -> Result<i64> {
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::new())]));
    let query = BoundaryQuery::new(identity, &target, EpsilonRange::new(0, 0)?, 20)?;
    Ok(index(cache.best(&query, policy, p)?.unwrap().boundary))
}
fn bank(identity: &BoundaryIdentity, s: Symbol, coordinates: &[i64]) -> RustFlowCache {
    let mut bank = RustFlowCache::default();
    for &coordinate in coordinates {
        bank.insert(entry(identity, s, Atom::num(coordinate), 40))
            .unwrap();
    }
    bank
}

#[test]
fn lower_bounds_price_before_guards_and_preserve_all_rounded_ties() {
    let (s, identity, p) = fixture();
    let cache = bank(&identity, s, &[4, 3, 2, 1]);
    let calls = RefCell::new(Vec::new());
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |source: &CachedBoundary, _: &CachedPoint| {
            let i = index(source);
            calls.borrow_mut().push(i);
            Ok(i != 1)
        },
    };
    assert_eq!(selected(&cache, &identity, s, p, &policy).unwrap(), 2);
    assert_eq!(*calls.borrow(), vec![1, 2]);

    let gap = Atom::num(10).pow(-120);
    let imaginary = Atom::num(symbolica::domains::float::Complex::new(
        Rational::from(0),
        Rational::from(1),
    ));
    let near = Atom::one() + &imaginary;
    let far = Atom::one() + &imaginary * (Atom::one() + gap);
    let mut complex = RustFlowCache::default();
    complex.insert(entry(&identity, s, far, 80)).unwrap();
    complex
        .insert(entry(&identity, s, near.clone(), 30))
        .unwrap();
    let calls = RefCell::new(0);
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| {
            *calls.borrow_mut() += 1;
            Ok(true)
        },
    };
    let target = CachedPoint::Exact(BTreeMap::from([(s, Atom::new())]));
    let query =
        BoundaryQuery::new(&identity, &target, EpsilonRange::new(0, 0).unwrap(), 20).unwrap();
    assert_eq!(
        complex
            .best(&query, &policy, p)
            .unwrap()
            .unwrap()
            .boundary
            .point
            .restart_coordinates()
            .unwrap()[&s],
        near
    );
    assert_eq!(*calls.borrow(), 2);
}

#[test]
fn conservative_and_unknown_bounds_do_not_act_as_exact_costs() {
    let (s, identity, p) = fixture();
    let cache = bank(&identity, s, &[1, 2, 3, 4]);
    let policy = costs(
        |i, p: Precision| Ok((i != 3).then(|| p.real(i))),
        |i, p: Precision| {
            Ok(Some(p.real(match i {
                1 => 20,
                2 => 2,
                3 => 10,
                _ => 4,
            })))
        },
    );
    assert_eq!(selected(&cache, &identity, s, p, &policy).unwrap(), 2);
    assert_eq!(*policy.calls.borrow(), vec![3, 1, 2]);
    let unknown = costs(
        |_, _: Precision| Ok(None),
        |i, p: Precision| Ok(Some(p.real(5 - i))),
    );
    assert_eq!(selected(&cache, &identity, s, p, &unknown).unwrap(), 4);
    assert_eq!(*unknown.calls.borrow(), vec![1, 2, 3, 4]);
}

#[test]
fn equal_actual_costs_keep_accuracy_and_insertion_ties_despite_bound_order() {
    let (s, identity, p) = fixture();
    let mut cache = bank(&identity, s, &[3, 1, 2]);
    let policy = costs(
        |i, p: Precision| Ok(Some(p.real(i))),
        |_, p: Precision| Ok(Some(p.real(8))),
    );
    assert_eq!(selected(&cache, &identity, s, p, &policy).unwrap(), 3);
    assert_eq!(*policy.calls.borrow(), vec![1, 2, 3]);
    cache.insert(entry(&identity, s, Atom::num(2), 50)).unwrap();
    assert_eq!(selected(&cache, &identity, s, p, &policy).unwrap(), 2);
    let equal = costs(
        |_, p: Precision| Ok(Some(p.real(8))),
        |_, p: Precision| Ok(Some(p.real(8))),
    );
    assert_eq!(selected(&cache, &identity, s, p, &equal).unwrap(), 2);
    assert_eq!(*equal.calls.borrow(), vec![3, 1, 2]);
}

#[test]
fn invalid_lower_bounds_and_observed_overestimates_are_errors() {
    let (s, identity, p) = fixture();
    let cache = bank(&identity, s, &[1, 2]);
    for value in [
        p.real(-1),
        Float::parse("nan", Some(p.bits)).unwrap(),
        Float::parse("inf", Some(p.bits)).unwrap(),
    ] {
        let policy = costs(
            |i, p: Precision| Ok(Some(if i == 1 { p.real(0) } else { value.clone() })),
            |_, p: Precision| Ok(Some(p.real(0))),
        );
        assert!(matches!(
            selected(&cache, &identity, s, p, &policy),
            Err(Error::InvalidInput(_))
        ));
        assert!(
            policy.calls.borrow().is_empty(),
            "all advertised bounds must be validated before pruning"
        );
    }
    let overestimate = costs(
        |_, p: Precision| Ok(Some(p.real(2))),
        |_, p: Precision| Ok(Some(p.real(1))),
    );
    assert!(
        matches!(selected(&cache, &identity, s, p, &overestimate), Err(Error::InvalidInput(message)) if message.contains("lower bound"))
    );
    assert!(selected(&cache, &identity, s, Precision { bits: 1 }, &overestimate).is_err());
}

#[test]
fn lower_bound_and_inspected_route_errors_propagate_without_cache_mutation() {
    let (s, identity, p) = fixture();
    let cache = bank(&identity, s, &[1, 2]);
    let cancelled = costs(
        |_, _: Precision| Err(Error::Cancelled),
        |_, p: Precision| Ok(Some(p.real(1))),
    );
    assert!(matches!(
        selected(&cache, &identity, s, p, &cancelled),
        Err(Error::Cancelled)
    ));
    let cancelled = costs(
        |_, p: Precision| Ok(Some(p.real(0))),
        |_, _: Precision| Err(Error::Cancelled),
    );
    assert!(matches!(
        selected(&cache, &identity, s, p, &cancelled),
        Err(Error::Cancelled)
    ));
    let failed = costs(
        |_, _: Precision| Err(Error::Numerical("bound failure".into())),
        |_, p: Precision| Ok(Some(p.real(1))),
    );
    assert!(matches!(
        selected(&cache, &identity, s, p, &failed),
        Err(Error::Numerical(_))
    ));
    assert_eq!(cache.len(), 2);
}

#[test]
fn insufficient_boundaries_are_filtered_before_bound_and_route_callbacks() {
    let (s, identity, p) = fixture();
    let mut cache = RustFlowCache::default();
    cache.insert(entry(&identity, s, Atom::num(1), 10)).unwrap();
    cache.insert(entry(&identity, s, Atom::num(2), 30)).unwrap();
    let policy = costs(
        |i, p: Precision| {
            assert_eq!(i, 2);
            Ok(Some(p.real(1)))
        },
        |i, p: Precision| {
            assert_eq!(i, 2);
            Ok(Some(p.real(1)))
        },
    );
    assert_eq!(selected(&cache, &identity, s, p, &policy).unwrap(), 2);
    assert_eq!(*policy.calls.borrow(), vec![2]);
}

#[test]
fn cancellation_during_cache_pricing_cannot_return_an_exact_hit() {
    let (s, identity, p) = fixture();
    let flow = RustFlow::new(
        KinematicSystem {
            epsilon: symbol!("cache_bounds::epsilon"),
            derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]])]),
        },
        &[parse!("cache_bounds::one")],
        &Atom::one(),
        Prescription::PlusI0,
        "constant analytic solution",
    )
    .unwrap();
    assert_eq!(flow.identity().key(), identity.key());
    let mut cache = bank(&identity, s, &[1]);
    let context = RunContext::default();
    let cancellation = context.cancellation.clone();
    let policy = costs(
        |_, _: Precision| {
            cancellation.cancel();
            Ok(Some(p.real(0)))
        },
        |_, _: Precision| Ok(Some(p.real(0))),
    );
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(1))]),
        EpsilonRange::new(0, 0).unwrap(),
        &FlowOptions::default(),
        &context,
        &policy,
    );
    assert!(matches!(result, Err(Error::Cancelled)));
    assert_eq!(cache.len(), 1);
}
