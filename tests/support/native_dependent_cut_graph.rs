use feynkit_graph::{
    DiagramEdge, DiagramEndpoint, DiagramHalfEdge, DiagramVertex, ExternalLeg, ExternalState,
    FeynmanDiagram,
};
use feynkit_model::Model;
use std::sync::Arc;

/// A two-point insertion on one massless final-state leg gives the dependent
/// smooth uncut pole r²-M² alongside cut r² and cut (P-r)².
pub fn diagram() -> FeynmanDiagram {
    let mut model: serde_json::Value =
        serde_json::from_str(include_str!("../../fixtures/hepkit/massless_phi3.json")).unwrap();
    model["name"] = "native_dependent_cut_model".into();
    model["parameters"].as_array_mut().unwrap().push(serde_json::json!({
        "name":"M","lhablock":null,"lhacode":null,"nature":"external","parameter_type":"real","value":[123.0,0.0],"expression":null
    }));
    let mut heavy = model["particles"][0].clone();
    heavy["name"] = "chi".into();
    heavy["antiname"] = "chi".into();
    heavy["pdg_code"] = 1001.into();
    heavy["mass"] = "M".into();
    model["particles"].as_array_mut().unwrap().push(heavy);
    let mut propagator = model["propagators"][0].clone();
    propagator["name"] = "chi_prop".into();
    propagator["particle"] = "chi".into();
    propagator["denominator"] = "(UFO::P(UFO::idx(1,1)))^2-UFO::M^2".into();
    model["propagators"]
        .as_array_mut()
        .unwrap()
        .push(propagator);
    model["lorentz_structures"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({"name":"scalar2","spins":[1,1],"structure":"1"}));
    model["vertex_rules"].as_array_mut().unwrap().extend([
        serde_json::json!({"name":"phichi","particles":["phi","chi"],"color_structures":["1"],"lorentz_structures":["scalar2"],"couplings":[[null]]}),
        serde_json::json!({"name":"phiphichi","particles":["phi","phi","chi"],"color_structures":["1"],"lorentz_structures":["scalar3"],"couplings":[[null]]})
    ]);
    let model = Arc::new(Model::from_json(&model.to_string()).unwrap());
    let phi = model.particle_id("phi").unwrap();
    let chi = model.particle_id("chi").unwrap();
    let mut builder = FeynmanDiagram::builder(model.clone(), "native_dependent_cut");
    let a = builder.add_vertex(DiagramVertex::interaction(
        "a",
        model.vertex_rule_id("phi3").unwrap(),
    ));
    let b = builder.add_vertex(DiagramVertex::interaction(
        "b",
        model.vertex_rule_id("phichi").unwrap(),
    ));
    let c = builder.add_vertex(DiagramVertex::interaction(
        "c",
        model.vertex_rule_id("phiphichi").unwrap(),
    ));
    let mut external = DiagramEdge::new(phi, false);
    external.external = Some(ExternalLeg {
        name: "p".into(),
        index: 0,
        state: ExternalState::Incoming,
        connection: 0,
    });
    let incoming = builder.add_edge(a, c, external).unwrap();
    let first = builder
        .add_edge(a, b, DiagramEdge::new(phi, false))
        .unwrap();
    let second = builder
        .add_edge(a, c, DiagramEdge::new(phi, false))
        .unwrap();
    let uncut = builder
        .add_edge(b, c, DiagramEdge::new(chi, false))
        .unwrap();
    let half = |edge, endpoint| DiagramHalfEdge { edge, endpoint };
    let mut left = [incoming, first, second]
        .map(|edge| half(edge, DiagramEndpoint::Target))
        .to_vec();
    left.extend([
        half(uncut, DiagramEndpoint::Source),
        half(uncut, DiagramEndpoint::Target),
    ]);
    let right = [incoming, first, second]
        .map(|edge| half(edge, DiagramEndpoint::Source))
        .to_vec();
    builder
        .build()
        .unwrap()
        .with_loop_momentum_edges(&[first])
        .unwrap()
        .with_cut_partitions(vec![(left, right)])
        .unwrap()
}
