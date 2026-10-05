use serde_json::{Value, json};
use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
};
use symbolica::prelude::*;
use symbolica_amflow::{ComplexFloat, Precision};

struct Temporary(PathBuf);
impl Temporary {
    fn new() -> Self {
        let id = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let result = Self(std::env::temp_dir().join(format!(
            "rustflow-cli-epsilon-shearing-{}-{id}",
            std::process::id()
        )));
        fs::create_dir_all(&result.0).unwrap();
        result
    }
    fn execute(&self, card: &Value) -> Output {
        let path = self.0.join("input.json");
        fs::write(&path, serde_json::to_vec(card).unwrap()).unwrap();
        Command::new(env!("CARGO_BIN_EXE_rustflow"))
            .arg("transport")
            .arg(path)
            .output()
            .unwrap()
    }
    fn run(&self, card: &Value) -> Value {
        let output = self.execute(card);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        serde_json::from_slice(&output.stdout).unwrap()
    }
    fn reject(&self, card: &Value, message: &str) {
        let output = self.execute(card);
        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}
impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn card() -> Value {
    let mut value: Value = serde_json::from_str(include_str!(
        "../examples/cli/epsilon-shearing-transport.json"
    ))
    .unwrap();
    value["cache_directory"] = "bank".into();
    value
}
fn number(v: &Value, p: Precision) -> ComplexFloat {
    p.parse(
        v["real"].as_str().unwrap(),
        v["imaginary"].as_str().unwrap(),
    )
    .unwrap()
}
fn coordinate() -> String {
    Atom::var(symbol!("rustflow_epsilon_cli::s")).to_canonical_string()
}

#[test]
fn epsilon_cli_restores_original_outputs_and_reuses_a_sheared_binary_bank() {
    let tmp = Temporary::new();
    let mut input = card();
    let output = tmp.run(&input);
    assert_eq!(output["epsilon_shearing"]["weights"], json!([-1, 0]));
    assert_eq!(
        output["identity"],
        output["epsilon_shearing"]["original_identity"]
    );
    assert_ne!(
        output["epsilon_shearing"]["cached_identity"],
        output["epsilon_shearing"]["original_identity"]
    );
    assert_eq!(
        output["epsilon_shearing"]["original_output_range"],
        json!({"leading":-1,"last":0})
    );
    assert_eq!(
        output["epsilon_shearing"]["requested_cached_range"],
        json!({"leading":-1,"last":1})
    );
    let p = Precision::decimal(80).unwrap();
    for (result, s) in output["results"].as_array().unwrap().iter().zip([1, 2, 2]) {
        assert_eq!(result["coefficients"].as_array().unwrap().len(), 2);
        assert_eq!(result["absolute_errors"].as_array().unwrap().len(), 2);
        let ex = p.exp(&p.i(s));
        let em = p.exp(&p.i(-s));
        let sinh = p.scale(&p.sub(&ex, &em), 1, 2);
        let cosh = p.scale(&p.add(&ex, &em), 1, 2);
        assert!(p.close(&number(&result["coefficients"][0][0], p), &sinh, 20));
        assert!(p.close(&number(&result["coefficients"][1][0], p), &sinh, 20));
        assert!(p.close(&number(&result["coefficients"][0][1], p), &p.zero(), 20));
        assert!(p.close(&number(&result["coefficients"][1][1], p), &cosh, 20));
    }
    assert_eq!(output["results"][1]["starting_point"][coordinate()], "1");
    assert_eq!(output["results"][2]["steps"], 0);
    input.as_object_mut().unwrap().remove("seeds");
    let restarted = tmp.run(&input);
    assert_eq!(
        restarted["loaded_boundaries"],
        output["retained_boundaries"]
    );
    assert_eq!(restarted["epsilon_shearing"], output["epsilon_shearing"]);
    assert!(
        restarted["results"]
            .as_array()
            .unwrap()
            .iter()
            .all(|r| r["steps"] == 0)
    );
    for (old, new) in output["results"]
        .as_array()
        .unwrap()
        .iter()
        .zip(restarted["results"].as_array().unwrap())
    {
        assert_eq!(old["coefficients"], new["coefficients"]);
        assert_eq!(old["absolute_errors"], new["absolute_errors"]);
    }
}

#[test]
fn epsilon_cli_rejects_missing_original_orders_even_when_a_cache_hit_exists() {
    let tmp = Temporary::new();
    let input = card();
    tmp.run(&input);
    let path = tmp.0.join("bank/physical-boundaries.bin");
    let saved = fs::read(&path).unwrap();
    let mut short = input.clone();
    short["seeds"][0]["last_epsilon_power"] = 0.into();
    short["seeds"][0]["coefficients"]
        .as_array_mut()
        .unwrap()
        .truncate(2);
    short["seeds"][0]["absolute_errors"]
        .as_array_mut()
        .unwrap()
        .truncate(2);
    tmp.reject(
        &short,
        "requires original seed coefficients through epsilon^1",
    );
    assert_eq!(fs::read(&path).unwrap(), saved);
    short["seeds"][0]
        .as_object_mut()
        .unwrap()
        .remove("last_epsilon_power");
    tmp.reject(
        &short,
        "requires original seed coefficients through epsilon^1",
    );
    assert_eq!(fs::read(&path).unwrap(), saved);
    let mut incompatible = input;
    incompatible["derivatives"] = json!({"s":[["1/eps","0"],["0","0"]]});
    tmp.reject(&incompatible, "negative constraint cycle");
    assert_eq!(fs::read(&path).unwrap(), saved);
}

#[test]
fn epsilon_cli_rejects_explicit_algebraic_and_canonical_shearing() {
    let tmp = Temporary::new();
    let mut algebraic = card();
    algebraic["roots"] = json!({"r":"s"});
    tmp.reject(&algebraic, "supported only for rational derivatives");
    assert!(!tmp.0.join("bank").exists());
    let mut canonical = card();
    canonical.as_object_mut().unwrap().remove("derivatives");
    canonical["canonical"] =
        json!({"variables":["s"],"letters":["s"],"matrices":[[["1","0"],["0","1"]]]});
    tmp.reject(&canonical, "algebraic and canonical inputs are unsupported");
    assert!(!tmp.0.join("bank").exists());
}

#[test]
fn epsilon_cli_preserves_prescribed_branch_values_and_metadata() {
    let tmp = Temporary::new();
    let mut input = card();
    input["derivatives"] = json!({"s":[["0","1/(4*eps*s)"],["eps/(4*s)","0"]]});
    input["seeds"][0]["coordinates"] = json!({"s":"1"});
    input["destinations"] = json!([{"s":"-1"},{"s":"-1"}]);
    input["continuation"] = json!({"kind":"prescribed_affine","domain":"principal log(s+i0) along the declared planner route","prescriptions":[{"polynomial":"s","side":"+i0"}],"unprescribed_side":"+i0","homotopy_admission":"all_planner_routes_in_declared_domain"});
    let result = tmp.run(&input);
    assert_eq!(
        result["identity"],
        result["epsilon_shearing"]["original_identity"]
    );
    assert_ne!(
        result["identity"],
        result["epsilon_shearing"]["cached_identity"]
    );
    assert_eq!(result["continuation"]["prescriptions"][0]["side"], "+i0");
    let p = Precision::decimal(80).unwrap();
    let half = ComplexFloat::new(p.real(2).sqrt() / 2, p.real(0));
    let positive = p.mul(&p.complex(0, 1), &half);
    for row in result["results"].as_array().unwrap() {
        assert!(p.close(&number(&row["coefficients"][0][0], p), &positive, 20));
        assert!(p.close(&number(&row["coefficients"][1][0], p), &positive, 20));
        assert!(p.close(&number(&row["coefficients"][1][1], p), &half, 20));
    }
    assert_eq!(result["results"][1]["steps"], 0);
    let file = tmp.0.join("bank/physical-boundaries.bin");
    let snapshot = fs::read(&file).unwrap();
    let mut opposite = input.clone();
    opposite.as_object_mut().unwrap().remove("seeds");
    opposite["continuation"]["prescriptions"][0]["side"] = "-i0".into();
    tmp.reject(&opposite, "no compatible cached physical boundary");
    assert_eq!(fs::read(&file).unwrap(), snapshot);
    opposite["seeds"] = input["seeds"].clone();
    let negative = tmp.run(&opposite);
    assert_ne!(negative["identity"], result["identity"]);
    assert!(p.close(
        &number(&negative["results"][0]["coefficients"][0][0], p),
        &p.neg(&positive),
        20
    ));
}
