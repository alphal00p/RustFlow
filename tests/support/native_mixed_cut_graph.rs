use feynkit_graph::{
    DiagramEdge, DiagramEndpoint, DiagramHalfEdge, DiagramVertex, ExternalLeg, ExternalState,
    FeynmanDiagram,
};
use feynkit_model::Model;
use std::sync::Arc;

/// Connected massive two-body cut with an internal massless virtual bubble.
pub fn connected_diagram() -> FeynmanDiagram {
    // A phi/chi scalar model: massive phi cuts and massless chi virtual lines.
    // Native model validation and native routing own this graph, as in user input.
    let mut model: serde_json::Value =
        serde_json::from_str(include_str!("../../fixtures/hepkit/massless_phi3.json")).unwrap();
    model["name"] = "native_mixed_cut_model".into();
    model["parameters"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({
            "name":"M", "lhablock":null,"lhacode":null,"nature":"external","parameter_type":"real",
            "value":[123.0,0.0],"expression":null
        }));
    let mut massless = model["particles"][0].clone();
    massless["name"] = "chi".into();
    massless["antiname"] = "chi".into();
    massless["pdg_code"] = 1001.into();
    model["particles"][0]["mass"] = "M".into();
    model["particles"].as_array_mut().unwrap().push(massless);
    let mut massless_prop = model["propagators"][0].clone();
    massless_prop["name"] = "chi_prop".into();
    massless_prop["particle"] = "chi".into();
    model["propagators"][0]["denominator"] = "(UFO::P(UFO::idx(1,1)))^2-UFO::M^2".into();
    model["propagators"]
        .as_array_mut()
        .unwrap()
        .push(massless_prop);
    model["lorentz_structures"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({"name":"scalar4","spins":[1,1,1,1],"structure":"1"}));
    model["vertex_rules"].as_array_mut().unwrap().extend([
        serde_json::json!({"name":"phichichi","particles":["phi","chi","chi"],"color_structures":["1"],"lorentz_structures":["scalar3"],"couplings":[[null]]}),
        serde_json::json!({"name":"phiphichichi","particles":["phi","phi","chi","chi"],"color_structures":["1"],"lorentz_structures":["scalar4"],"couplings":[[null]]}),
    ]);
    let model = Arc::new(Model::from_json(&model.to_string()).unwrap());
    let phi = model.particle_id("phi").unwrap();
    let chi = model.particle_id("chi").unwrap();
    let mut builder = FeynmanDiagram::builder(model.clone(), "connected_native_real_virtual");
    let a = builder.add_vertex(DiagramVertex::interaction(
        "a",
        model.vertex_rule_id("phi3").unwrap(),
    ));
    let b = builder.add_vertex(DiagramVertex::interaction(
        "b",
        model.vertex_rule_id("phichichi").unwrap(),
    ));
    let c = builder.add_vertex(DiagramVertex::interaction(
        "c",
        model.vertex_rule_id("phiphichichi").unwrap(),
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
    let virtual_first = builder
        .add_edge(b, c, DiagramEdge::new(chi, false))
        .unwrap();
    let virtual_second = builder
        .add_edge(b, c, DiagramEdge::new(chi, false))
        .unwrap();
    let half = |edge, endpoint| DiagramHalfEdge { edge, endpoint };
    let mut left = [incoming, first, second]
        .map(|edge| half(edge, DiagramEndpoint::Target))
        .to_vec();
    for edge in [virtual_first, virtual_second] {
        left.push(half(edge, DiagramEndpoint::Source));
        left.push(half(edge, DiagramEndpoint::Target));
    }
    let right = [incoming, first, second]
        .map(|edge| half(edge, DiagramEndpoint::Source))
        .to_vec();
    builder
        .build()
        .unwrap()
        .with_loop_momentum_edges(&[first, virtual_first])
        .unwrap()
        .with_cut_partitions(vec![(left, right)])
        .unwrap()
}
