use feynkit_graph::{
    DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramHalfEdge, DiagramVertex,
    ExternalLeg, ExternalState, FeynmanDiagram, symbols,
};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock};
use symbolica::prelude::*;
pub fn model() -> Arc<Model> {
    static MODEL: OnceLock<Arc<Model>> = OnceLock::new();
    MODEL
        .get_or_init(|| {
            Arc::new(
                Model::from_json(include_str!("../../fixtures/hepkit/massless_phi3.json")).unwrap(),
            )
        })
        .clone()
}
pub fn diagram() -> FeynmanDiagram {
    let model = model();
    let rule = model.vertex_rule_id("phi3").unwrap();
    let particle = model.particle_id("phi").unwrap();
    let mut builder = FeynmanDiagram::builder(model, "native_two_body_cut");
    let a = builder.add_vertex(DiagramVertex::interaction("a", rule));
    let b = builder.add_vertex(DiagramVertex::interaction("b", rule));
    let mut external = DiagramEdge::new(particle, false);
    external.external = Some(ExternalLeg {
        name: "p".into(),
        index: 0,
        state: ExternalState::Incoming,
        connection: 0,
    });
    let incoming = builder.add_edge(a, b, external).unwrap();
    let first = builder
        .add_edge(a, b, DiagramEdge::new(particle, false))
        .unwrap();
    let second = builder
        .add_edge(a, b, DiagramEdge::new(particle, false))
        .unwrap();
    let half = |edge, endpoint| DiagramHalfEdge { edge, endpoint };
    let side = |endpoint| DiagramCutSide {
        half_edges: [incoming, first, second]
            .map(|edge| half(edge, endpoint))
            .to_vec(),
        coupling_orders: BTreeMap::new(),
        loop_count: 0,
    };
    builder
        .cuts(vec![DiagramCut {
            cut: vec![
                half(first, DiagramEndpoint::Target),
                half(second, DiagramEndpoint::Target),
            ],
            left: side(DiagramEndpoint::Target),
            right: side(DiagramEndpoint::Source),
        }])
        .build()
        .unwrap()
}
pub fn kinematics() -> Kinematics {
    Kinematics::in_dimension(&Atom::var(symbol!("native_cut_test::D")))
        .unwrap()
        .with_mass_squared(&symbols::external_momentum().call(0), Atom::num(25))
        .unwrap()
}
