use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
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
        let result =
            Self(std::env::temp_dir().join(format!("rustflow-cli-{}-{id}", std::process::id())));
        fs::create_dir_all(&result.0).unwrap();
        result
    }
    fn run(&self, command: &str, card: &Value) -> Value {
        let path = self.0.join("input.json");
        fs::write(&path, serde_json::to_vec(card).unwrap()).unwrap();
        let result = Command::new(env!("CARGO_BIN_EXE_rustflow"))
            .arg(command)
            .arg(path)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        serde_json::from_slice(&result.stdout).unwrap()
    }
}
impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn number(value: &Value, p: Precision) -> ComplexFloat {
    p.parse(
        value["real"].as_str().unwrap(),
        value["imaginary"].as_str().unwrap(),
    )
    .unwrap()
}

#[test]
fn native_dot_graph_cli_matches_the_massless_bubble() {
    let tmp = Temporary::new();
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut card: Value =
        serde_json::from_str(include_str!("../examples/cli/massless-bubble.json")).unwrap();
    card["model"] = root
        .join("fixtures/hepkit/massless_phi3.json")
        .display()
        .to_string()
        .into();
    card["diagram"] = root
        .join("fixtures/hepkit/massless_bubble.dot")
        .display()
        .to_string()
        .into();
    let result = tmp.run("graph", &card);
    assert_eq!(result["verified_digits"], 20);
    let p = Precision::decimal(80).unwrap();
    let gamma = ComplexFloat::new(
        Float::from_raw(rug::Float::with_val(p.bits, rug::float::Constant::Euler)),
        p.real(0),
    );
    let finite = p.sub(&p.sub(&p.i(2), &gamma), &p.log(&p.i(2)));
    assert!(p.close(&number(&result["coefficients"]["-2"], p), &p.zero(), 20));
    assert!(p.close(&number(&result["coefficients"]["-1"], p), &p.i(1), 20));
    assert!(p.close(&number(&result["coefficients"]["0"], p), &finite, 20));
}

#[test]
fn cli_reuses_a_growing_binary_bank_across_processes() {
    let tmp = Temporary::new();
    let mut card: Value =
        serde_json::from_str(include_str!("../examples/cli/physical-transport.json")).unwrap();
    card["cache_directory"] = "bank".into();
    let result = tmp.run("transport", &card);
    assert_eq!(result["loaded_boundaries"], 0);
    assert!(result["retained_boundaries"].as_u64().unwrap() > 3);
    let points = result["results"].as_array().unwrap();
    assert!(points[0]["steps"].as_u64().unwrap() > 0);
    assert!(points[1]["steps"].as_u64().unwrap() > 0);
    assert_eq!(points[2]["steps"], 0);
    let coordinate = Atom::var(symbol!("feynkit_graph::s")).to_canonical_string();
    assert_eq!(points[1]["starting_point"][&coordinate], "4");
    let p = Precision::decimal(80).unwrap();
    for (point, s) in points.iter().zip([4, 5, 5]) {
        assert!(point["verified_digits"].as_u64().unwrap() >= 20);
        let log = p.log(&p.i(s));
        assert!(p.close(&number(&point["coefficients"][0][0], p), &p.i(1), 20));
        assert!(p.close(&number(&point["coefficients"][1][0], p), &log, 20));
        assert!(p.close(
            &number(&point["coefficients"][2][0], p),
            &p.scale(&p.mul(&log, &log), 1, 2),
            20
        ));
    }
    assert!(tmp.0.join("bank/physical-boundaries.bin").is_file());
    card.as_object_mut().unwrap().remove("seeds");
    let restarted = tmp.run("transport", &card);
    assert_eq!(
        restarted["loaded_boundaries"],
        result["retained_boundaries"]
    );
    for point in restarted["results"].as_array().unwrap() {
        assert_eq!(point["steps"], 0);
    }
}
