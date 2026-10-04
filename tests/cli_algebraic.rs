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
            "rustflow-cli-algebraic-{}-{id}",
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
    let mut value: Value =
        serde_json::from_str(include_str!("../examples/cli/algebraic-transport.json")).unwrap();
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
    Atom::var(symbol!("rustflow_algebraic_cli::s")).to_canonical_string()
}
fn root() -> String {
    Atom::var(symbol!("rustflow_algebraic_cli::r")).to_canonical_string()
}

#[test]
fn algebraic_cli_progresses_restarts_and_keeps_exact_germs() {
    let tmp = Temporary::new();
    let mut input = card();
    let result = tmp.run(&input);
    assert_eq!(result["loaded_boundaries"], 0);
    assert!(result["retained_boundaries"].as_u64().unwrap() > 3);
    let points = result["results"].as_array().unwrap();
    assert!(points[0]["steps"].as_u64().unwrap() > 0);
    assert!(points[1]["steps"].as_u64().unwrap() > 0);
    assert_eq!(points[1]["starting_point"][coordinate()], "4");
    assert_eq!(points[2]["steps"], 0);
    let p = Precision::decimal(80).unwrap();
    for (point, coefficients) in points.iter().zip([[1, 2, 2], [1, 4, 8], [1, 4, 8]]) {
        assert!(point["verified_digits"].as_u64().unwrap() >= 20);
        assert_eq!(point["root_germ"][root()], "principal");
        assert_eq!(point["starting_root_germ"][root()], "principal");
        for (actual, expected) in point["coefficients"]
            .as_array()
            .unwrap()
            .iter()
            .zip(coefficients)
        {
            assert!(p.close(&number(&actual[0], p), &p.i(expected), 20));
        }
    }
    assert!(tmp.0.join("bank/physical-boundaries.bin").is_file());
    input.as_object_mut().unwrap().remove("seeds");
    let restarted = tmp.run(&input);
    assert_eq!(
        restarted["loaded_boundaries"],
        result["retained_boundaries"]
    );
    assert!(
        restarted["results"]
            .as_array()
            .unwrap()
            .iter()
            .all(|p| p["steps"] == 0)
    );
}

#[test]
fn algebraic_cli_never_substitutes_an_opposite_sheet_cache_hit() {
    let tmp = Temporary::new();
    let input = card();
    tmp.run(&input);
    let file = tmp.0.join("bank/physical-boundaries.bin");
    let original = fs::read(&file).unwrap();
    let mut opposite = input.clone();
    opposite.as_object_mut().unwrap().remove("seeds");
    opposite["destinations"] = json!([{"coordinates":{"s":"9"},"root_germ":{"r":"opposite"}}]);
    tmp.reject(&opposite, "no compatible cached physical boundary");
    assert_eq!(fs::read(&file).unwrap(), original);
    // An explicit independent opposite-sheet seed enables the requested sheet.
    // Both sheets have exact [1,0,0] at s=1 for their normalized exponential.
    opposite["seeds"] = input["seeds"].clone();
    opposite["seeds"][0]["root_germ"]["r"] = "opposite".into();
    opposite["seeds"][0]["provenance"] =
        "Analytic exp(-2*eps*(sqrt(s)-1)): exact coefficients at s=1".into();
    let result = tmp.run(&opposite);
    let point = &result["results"][0];
    assert!(point["steps"].as_u64().unwrap() > 0);
    assert_eq!(point["root_germ"][root()], "opposite");
    let p = Precision::decimal(80).unwrap();
    assert!(p.close(&number(&point["coefficients"][1][0], p), &p.i(-4), 20));
    assert!(p.close(&number(&point["coefficients"][2][0], p), &p.i(8), 20));
    let mut principal = input;
    principal.as_object_mut().unwrap().remove("seeds");
    principal["destinations"] = json!([{"coordinates":{"s":"9"},"root_germ":{"r":"principal"}}]);
    let result = tmp.run(&principal);
    assert_eq!(result["results"][0]["steps"], 0);
    assert!(p.close(
        &number(&result["results"][0]["coefficients"][1][0], p),
        &p.i(4),
        20
    ));
}

#[test]
fn algebraic_cli_rejects_missing_foreign_and_duplicate_germs() {
    let tmp = Temporary::new();
    let input = card();
    let mut missing = input.clone();
    missing["seeds"][0]
        .as_object_mut()
        .unwrap()
        .remove("root_germ");
    tmp.reject(&missing, "requires an explicit root_germ");
    let mut missing = input.clone();
    missing["destinations"] = json!([{"s":"4"}]);
    tmp.reject(&missing, "requires an explicit root_germ");
    let mut foreign = input.clone();
    foreign["seeds"][0]["root_germ"] = json!({"q":"principal"});
    tmp.reject(&foreign, "must name every registered root exactly once");
    let mut duplicate = input.clone();
    duplicate["seeds"][0]["root_germ"] =
        json!({"r":"principal","rustflow_algebraic_cli::r":"opposite"});
    tmp.reject(&duplicate, "names the same root twice");
    let mut no_roots = input;
    no_roots.as_object_mut().unwrap().remove("roots");
    no_roots["derivatives"] = json!({"s":[["eps/s"]]});
    tmp.reject(&no_roots, "requires a nonempty roots registry");
}

#[test]
fn algebraic_cli_preserves_completed_destinations_if_a_later_path_crosses_a_root() {
    let tmp = Temporary::new();
    let mut input = card();
    input["destinations"] = json!([
        {"coordinates":{"s":"4"},"root_germ":{"r":"principal"}},
        {"coordinates":{"s":"-1"},"root_germ":{"r":"principal"}},
    ]);
    tmp.reject(&input, "no compatible cached physical boundary");
    assert!(tmp.0.join("bank/physical-boundaries.bin").is_file());
    input.as_object_mut().unwrap().remove("seeds");
    input["destinations"] = json!([{"coordinates":{"s":"4"},"root_germ":{"r":"principal"}}]);
    let result = tmp.run(&input);
    assert!(result["loaded_boundaries"].as_u64().unwrap() > 1);
    assert_eq!(result["results"][0]["steps"], 0);
}

#[test]
fn algebraic_cli_extension_keeps_legacy_rational_input_and_output() {
    let tmp = Temporary::new();
    let mut input: Value =
        serde_json::from_str(include_str!("../examples/cli/physical-transport.json")).unwrap();
    input["cache_directory"] = "rational-bank".into();
    let result = tmp.run(&input);
    let p = Precision::decimal(70).unwrap();
    for (point, s) in result["results"].as_array().unwrap().iter().zip([4, 5, 5]) {
        assert!(point.get("root_germ").is_none());
        assert!(point.get("starting_root_germ").is_none());
        assert!(p.close(
            &number(&point["coefficients"][1][0], p),
            &p.log(&p.i(s)),
            20
        ));
    }
}
