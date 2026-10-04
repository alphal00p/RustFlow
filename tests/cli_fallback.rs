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

#[test]
fn cli_boundary_budget_and_accuracy_retry_preserve_existing_binary_bank() {
    let tmp = Temporary::new();
    let p = Precision::decimal(80).unwrap();
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let near = p.exp(
        &p.eval(
            &(imaginary.clone() * Atom::num(Rational::from((19, 10)))),
            &Default::default(),
        )
        .unwrap(),
    );
    let seed = |coordinate: &str, value: &ComplexFloat, digits| {
        json!({
            "coordinates":{"s":coordinate},"working_digits":80,"verified_digits":digits,
            "coefficients":[[{"real":value.re.as_raw().to_string(),"imaginary":value.im.as_raw().to_string()}]],
            "absolute_errors":[["1e-70"]],"provenance":"Independent analytic exp(i*s), rounded at80 digits with conservative source evidence",
        })
    };
    let mut input = json!({
        "schema_version":1,"namespace":"rustflow_fallback_cli","derivatives":{"s":[[imaginary.to_canonical_string()]]},
        "basis":["Y"],"normalization":"1","branch_domain":"real s; exp(i*s)","leading_epsilon_power":0,"last_epsilon_power":0,
        "cache_directory":"bank","seeds":[seed("19/10",&near,20),seed("0",&p.i(1),60)],
        "destinations":[{"s":"0"}],"options":{"digits":20,"guard_digits":30,"series_order":64,"max_boundary_attempts":1},
    });
    tmp.run(&input);
    let file = tmp.0.join("bank/physical-boundaries.bin");
    let original = fs::read(&file).unwrap();
    input.as_object_mut().unwrap().remove("seeds");
    input["destinations"] = json!([{"s":"2"}]);
    tmp.reject(&input, "budget 1 exhausted");
    assert_eq!(original, fs::read(&file).unwrap());
    input["options"]["max_boundary_attempts"] = 8.into();
    let result = tmp.run(&input);
    let point = &result["results"][0];
    let attempts = point["boundary_attempts"].as_array().unwrap();
    assert_eq!(attempts.len(), 2);
    assert_eq!(attempts[0]["outcome"]["status"], "accuracy_rejected");
    assert_eq!(attempts[1]["outcome"]["status"], "accepted");
    assert_eq!(attempts[0]["source_verified_digits"], 20);
    assert_eq!(attempts[1]["source_verified_digits"], 60);
    let actual = &point["coefficients"][0][0];
    let actual = p
        .parse(
            actual["real"].as_str().unwrap(),
            actual["imaginary"].as_str().unwrap(),
        )
        .unwrap();
    assert!(p.close(&actual, &p.exp(&p.parse("0", "2").unwrap()), 20));
    let repeated = tmp.run(&input);
    assert_eq!(repeated["results"][0]["steps"], 0);
    assert!(repeated["results"][0].get("boundary_attempts").is_none());
}
