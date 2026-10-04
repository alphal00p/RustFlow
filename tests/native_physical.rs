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
        &weights.keys().cloned().collect::<Vec<_>>(),
        &[s],
        &backend,
        &options,
        "negative real external invariant",
        &context,
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![1, 1])]);
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
    for destination in [-2, -3, -3] {
        let count = cache.len();
        let result = prepared
            .flow()
            .evaluate_to(
                &mut cache,
                &BTreeMap::from([(s, Atom::num(destination))]),
                seed.range,
                &options,
                &context,
                &policy,
            )
            .unwrap();
        assert!(result.boundary.accuracy.verified_digits() >= 20);
        assert!(p.close(&result.boundary.coefficients[0][0], &p.zero(), 20));
        assert!(p.close(&result.boundary.coefficients[1][0], &p.i(1), 20));
        let finite = p.sub(&p.sub(&p.i(2), &gamma), &p.log(&p.i(-destination)));
        assert!(p.close(&result.boundary.coefficients[2][0], &finite, 20));
        if result.transport.is_none() {
            assert_eq!(destination, -3);
            assert_eq!(cache.len(), count);
        } else {
            assert!(cache.len() > count);
        }
    }
}
