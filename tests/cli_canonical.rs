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
            "rustflow-cli-canonical-{}-{id}",
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
        serde_json::from_str(include_str!("../examples/cli/canonical-transport.json")).unwrap();
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
    Atom::var(symbol!("rustflow_canonical_cli::s")).to_canonical_string()
}
fn root() -> String {
    Atom::var(symbol!("rustflow_canonical_cli::r")).to_canonical_string()
}

#[test]
fn canonical_cli_reports_native_identity_and_restarts_progressive_root_transport() {
    use symbolica_amflow::algebraic::{CanonicalAlgebraicSystem, SquareRoot};
    use symbolica_amflow::{Prescription, RustFlow};
    let tmp = Temporary::new();
    let mut input = card();
    let result = tmp.run(&input);
    assert_eq!(result["representation"], "canonical");
    let eps = symbol!("rustflow_canonical_cli::eps");
    let s = symbol!("rustflow_canonical_cli::s");
    let r = symbol!("rustflow_canonical_cli::r");
    let flow = RustFlow::new_canonical(
        CanonicalAlgebraicSystem::new(
            eps,
            &[s],
            &[Atom::var(s) + Atom::var(r) + 1],
            &[vec![vec![Atom::one()]]],
            vec![SquareRoot {
                symbol: r,
                radicand: Atom::var(s),
            }],
        )
        .unwrap(),
        &[Atom::var(symbol!("rustflow_canonical_cli::I"))],
        &Atom::one(),
        Prescription::PlusI0,
        input["branch_domain"].as_str().unwrap(),
    )
    .unwrap();
    assert_eq!(result["identity"], flow.identity().key());
    let points = result["results"].as_array().unwrap();
    assert_eq!(points[1]["starting_point"][coordinate()], "4");
    assert_eq!(points[2]["steps"], 0);
    let p = Precision::decimal(80).unwrap();
    for (point, ratio) in points.iter().zip([7, 13, 13]) {
        let logarithm = p.log(&p.div(&p.i(ratio), &p.i(3)));
        assert_eq!(point["root_germ"][root()], "principal");
        assert!(p.close(&number(&point["coefficients"][1][0], p), &logarithm, 20));
        assert!(p.close(
            &number(&point["coefficients"][2][0], p),
            &p.div(&p.mul(&logarithm, &logarithm), &p.i(2)),
            20
        ));
    }
    input.as_object_mut().unwrap().remove("seeds");
    let restarted = tmp.run(&input);
    assert_eq!(restarted["identity"], result["identity"]);
    assert!(
        restarted["results"]
            .as_array()
            .unwrap()
            .iter()
            .all(|p| p["steps"] == 0)
    );
}

#[test]
fn rootless_canonical_cli_keeps_distinct_identity_and_equivalent_dense_values() {
    let tmp = Temporary::new();
    let mut input = card();
    input.as_object_mut().unwrap().remove("roots");
    input["canonical"]["letters"] = json!(["s"]);
    input["seeds"][0]
        .as_object_mut()
        .unwrap()
        .remove("root_germ");
    input["destinations"] = json!([{"s":"4"},{"s":"9"},{"s":"9"}]);
    let mut dense = input.clone();
    dense.as_object_mut().unwrap().remove("canonical");
    dense["derivatives"] = json!({"s":[["eps/s"]]});
    let legacy = tmp.run(&dense);
    assert!(legacy.get("representation").is_none());
    assert!(legacy.get("identity").is_none());
    let mut unseeded = input.clone();
    unseeded.as_object_mut().unwrap().remove("seeds");
    tmp.reject(&unseeded, "no compatible cached physical boundary");
    let canonical = tmp.run(&input);
    let p = Precision::decimal(80).unwrap();
    for (a, b) in canonical["results"]
        .as_array()
        .unwrap()
        .iter()
        .zip(legacy["results"].as_array().unwrap())
    {
        assert!(a.get("root_germ").is_none());
        for (a, b) in a["coefficients"]
            .as_array()
            .unwrap()
            .iter()
            .zip(b["coefficients"].as_array().unwrap())
        {
            assert!(p.close(&number(&a[0], p), &number(&b[0], p), 20));
        }
    }
    let restarted = tmp.run(&unseeded);
    assert!(
        restarted["results"]
            .as_array()
            .unwrap()
            .iter()
            .all(|p| p["steps"] == 0)
    );
}

#[test]
fn canonical_cli_rejects_ambiguous_missing_or_invalid_representations() {
    let tmp = Temporary::new();
    let input = card();
    let mut both = input.clone();
    both["derivatives"] = json!({"s":[["eps/s"]]});
    tmp.reject(&both, "exactly one of derivatives or canonical");
    both["derivatives"] = json!({});
    tmp.reject(&both, "exactly one of derivatives or canonical");
    let mut neither = input.clone();
    neither.as_object_mut().unwrap().remove("canonical");
    tmp.reject(&neither, "exactly one of derivatives or canonical");
    let mut missing_germ = input.clone();
    missing_germ["destinations"] = json!([{"s":"4"}]);
    tmp.reject(&missing_germ, "requires an explicit root_germ");
    let mut invalid = input;
    invalid["canonical"]["matrices"] = json!([[["s"]]]);
    tmp.reject(&invalid, "undeclared symbols");
}
