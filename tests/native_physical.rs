use feynkit_graph::symbols;
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use std::{collections::BTreeMap, sync::Arc};
use symbolica::prelude::*;
use symbolica_amflow::{hepkit::GraphIntegral, transport_cache::*, *};

#[test]
fn native_graph_amf_seed_drives_a_progressively_filled_physical_cache() {
    let model =
        Arc::new(Model::from_json(include_str!("../fixtures/hepkit/massless_phi3.json")).unwrap());
    let s = symbol!("native_physical::s");
    let epsilon = symbol!("native_physical::eps");
    let kinematics = Kinematics::in_dimension(&parse!("native_physical::D"))
        .unwrap()
        .with_mass_squared(&symbols::external_momentum().call(1), Atom::var(s))
        .unwrap();
    let graph = GraphIntegral::from_dot(
        model,
        include_str!("../fixtures/hepkit/massless_bubble.dot"),
        &kinematics,
    )
    .unwrap();
    let context = RunContext::default();
    let groups = graph
        .integral_groups(&KinematicPoint::default(), epsilon, 4, 1000, &context)
        .unwrap();
    assert_eq!(groups.len(), 1);
    let (family, weights) = &groups[0];
    assert_eq!(
        weights,
        &BTreeMap::from([(Integral(vec![1, 1]), Atom::one())])
    );
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let prepared = PreparedPhysicalFamily::new(
        family,
        &[Integral(vec![1, 1]), Integral(vec![2, 1])],
        &[s],
        &backend,
        &options,
        "negative real external invariant",
        &context,
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![1, 1])]);
    assert_eq!(
        prepared
            .required_master_range(
                &BTreeMap::from([(s, Atom::num(-1))]),
                EpsilonRange::new(-2, 0).unwrap()
            )
            .unwrap(),
        EpsilonRange::new(-2, 0).unwrap()
    );
    let mut cache = RustFlowCache::default();
    let seed = prepared
        .seed_cache(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(-1))]),
            0,
            16,
            &options,
            &backend,
            &context,
        )
        .unwrap();
    assert!(seed.accuracy.verified_digits() >= 36);
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let p = Precision::decimal(80).unwrap();
    let gamma = ComplexFloat::new(
        Float::from_raw(rug::Float::with_val(p.bits, rug::float::Constant::Euler)),
        p.real(0),
    );
    // Preserve a reserve for the uncertainty amplified by target reduction.
    let transport_options = FlowOptions {
        digits: 28,
        ..options.clone()
    };
    for destination in [-2, -3, -3] {
        let count = cache.len();
        let result = prepared
            .flow()
            .evaluate_to(
                &mut cache,
                &BTreeMap::from([(s, Atom::num(destination))]),
                seed.range,
                &transport_options,
                &context,
                &policy,
            )
            .unwrap();
        assert!(result.boundary.accuracy.verified_digits() >= 20);
        assert!(p.close(&result.boundary.coefficients[0][0], &p.zero(), 20));
        assert!(p.close(&result.boundary.coefficients[1][0], &p.i(1), 20));
        let finite = p.sub(&p.sub(&p.i(2), &gamma), &p.log(&p.i(-destination)));
        assert!(p.close(&result.boundary.coefficients[2][0], &finite, 20));
        let targets = prepared
            .project_targets(&result.boundary, seed.range, 20)
            .unwrap();
        assert!(p.close(&targets[0].coefficients[&0], &finite, 20));
        // I(2,1) = -(1-2*epsilon)/s I(1,1), including its finite shift.
        assert!(p.close(
            &targets[1].coefficients[&-1],
            &p.scale(&p.i(1), 1, -destination),
            20
        ));
        assert!(p.close(
            &targets[1].coefficients[&0],
            &p.scale(&p.sub(&finite, &p.i(2)), 1, -destination),
            20
        ));
        if result.transport.is_none() {
            assert_eq!(destination, -3);
            assert_eq!(cache.len(), count);
        } else {
            assert!(cache.len() > count);
        }
    }
}
