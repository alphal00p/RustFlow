use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{transport_cache::*, *};
#[test]
fn unequal_mass_two_loop_amf_seed_transports_and_matches_direct_evaluation() -> Result<()> {
    let (x, y, eps) = symbol!(
        "generic_multimass::x",
        "generic_multimass::y",
        "generic_multimass::eps"
    );
    let family = IntegralFamily {
        name: "generic_unequal_sunset".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::var(x), &[])?,
            Propagator::quadratic(&[0, 1], &[], Atom::var(y), &[])?,
            Propagator::quadratic(&[1, -1], &[], Atom::new(), &[])?,
        ],
        physical_propagators: 3,
        epsilon: eps,
        dimension: 4,
    };
    let context = RunContext::default();
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let targets = [Integral(vec![1, 1, 1]), Integral(vec![2, 1, 1])];
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &targets,
        &[x, y],
        &backend,
        &options,
        "positive unequal vacuum masses",
        &context,
    )?;
    let mut cache = RustFlowCache::default();
    let source = BTreeMap::from([(x, Atom::one()), (y, Atom::num(2))]);
    let destination = BTreeMap::from([
        (x, Atom::num(Rational::from((3, 2)))),
        (y, Atom::num(Rational::from((5, 2)))),
    ]);
    let range = EpsilonRange::new(-4, 0)?;
    let required = prepared.required_master_range(&destination, range)?;
    let seed = prepared.seed_cache(
        &mut cache,
        &source,
        required.last,
        16,
        &options,
        &backend,
        &context,
    )?;
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let checked = FlowOptions {
        digits: 28,
        ..options.clone()
    };
    let result = prepared.flow().evaluate_to(
        &mut cache,
        &destination,
        seed.range,
        &checked,
        &context,
        &policy,
    )?;
    let projected = prepared.project_targets(&result.boundary, range, 20)?;
    let point = KinematicPoint(
        destination
            .iter()
            .map(|(k, v)| (Atom::var(*k), v.clone()))
            .collect(),
    );
    let direct = solve_integrals(&family, &targets, &point, 0, &options, &backend, &context)?;
    let p = Precision::decimal(80)?;
    // The same public reduction and series interfaces handle this connected
    // family; no supplied master values or family-specific transport kernels.
    for (i, (a, b)) in projected.iter().zip(&direct).enumerate() {
        assert!(a.verified_digits >= 20);
        assert!(b.verified_digits.is_some_and(|d| d >= 20));
        for (power, value) in &a.coefficients {
            let independent = &b.coefficients[power];
            assert!(
                p.close(value, independent, 20),
                "target {i}, epsilon^{power}: {value} != {independent}"
            );
        }
        // These vacuum integrals have a double pole, not the conservative 1/eps^4
        // bound used when fitting a general two-loop family.
        for power in [-4, -3] {
            assert!(p.close(&a.coefficients[&power], &p.zero(), 20));
        }
    }
    assert!(p.close(&projected[0].coefficients[&-2], &p.i(2), 20));
    assert!(p.close(&projected[1].coefficients[&-2], &p.scale(&p.i(1), 1, 2), 20));
    assert!(cache.len() > 1);
    let count = cache.len();
    let hit = prepared.flow().evaluate_to(
        &mut cache,
        &destination,
        seed.range,
        &checked,
        &context,
        &policy,
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(cache.len(), count);
    Ok(())
}
