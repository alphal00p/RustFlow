//! Reuse physical boundaries for Y(s, epsilon) = s^epsilon on the positive axis.
//! Run with `cargo run --release --example physical_cache`.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::physical_transport::RustFlow;
use symbolica_amflow::transport_cache::{
    BoundaryAccuracy, CachedBoundary, CachedPoint, EpsilonRange, PointKind, RustFlowCache,
    ScaledDistance,
};
use symbolica_amflow::{FlowOptions, Precision, Prescription, Result, RunContext};

fn main() -> Result<()> {
    let (s, epsilon) = symbol!("physical_cache_example::s", "physical_cache_example::eps");
    let system = KinematicSystem::canonical_dlog(
        epsilon,
        &[s],
        &[Atom::var(s)],
        &[vec![vec![Atom::num(1)]]],
    )?;
    let transport = RustFlow::new(
        system,
        &[Atom::var(symbol!("physical_cache_example::Y"))],
        &Atom::num(1),
        Prescription::PlusI0,
        "positive real s; principal logarithm",
    )?;
    let range = EpsilonRange::new(0, 4)?;
    let p = Precision::decimal(90)?;
    let mut coefficients = vec![vec![p.zero()]; 5];
    coefficients[0][0] = p.i(1);
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary {
        identity: transport.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, Atom::num(1))])),
        kind: PointKind::Physical,
        range,
        coefficients,
        accuracy: BoundaryAccuracy::supplied(
            80,
            p.bits,
            vec![vec![p.real(0)]; 5],
            "exact analytic seed Y(1, epsilon)=1",
        )?,
    })?;

    // A straight path between positive real endpoints stays away from s=0
    // and remains on this logarithm branch. General problems need their own
    // path-admissibility checks; proximity alone is insufficient.
    let policy = ScaledDistance {
        scales: BTreeMap::from([(s, Atom::num(1))]),
        admissible: |source: &CachedBoundary, target: &CachedPoint| {
            let p = Precision::decimal(60)?;
            for point in [&source.point, target] {
                let value = &point.evaluate(p)?[&s];
                if value.im != p.real(0) || value.re <= p.real(0) {
                    return Ok(false);
                }
            }
            Ok(true)
        },
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let first = transport.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(16))]),
        range,
        &options,
        &context,
        &policy,
    )?;
    let cold_steps = first.transport.as_ref().unwrap().diagnostics.steps;
    println!(
        "s=16: start={}, steps={}, retained physical boundaries={}",
        first.starting_point.evaluate(p)?[&s],
        cold_steps,
        cache.len()
    );

    let destination = BTreeMap::from([(s, Atom::num(17))]);
    let nearby =
        transport.evaluate_to(&mut cache, &destination, range, &options, &context, &policy)?;
    let reused_steps = nearby.transport.as_ref().unwrap().diagnostics.steps;
    println!(
        "s=17: selected cached start={}, steps={}, verified digits={}",
        nearby.starting_point.evaluate(p)?[&s],
        reused_steps,
        nearby.boundary.accuracy.verified_digits()
    );
    assert!(reused_steps < cold_steps);
    assert!(p.close(
        &nearby.boundary.coefficients[1][0],
        &p.log(&p.i(17)),
        options.digits
    ));

    let repeated =
        transport.evaluate_to(&mut cache, &destination, range, &options, &context, &policy)?;
    assert!(repeated.transport.is_none());
    assert_eq!(repeated.boundary.coefficients, nearby.boundary.coefficients);
    println!(
        "s=17 again: selected cached start={}, exact-coordinate hit, steps=0",
        repeated.starting_point.evaluate(p)?[&s]
    );
    Ok(())
}
