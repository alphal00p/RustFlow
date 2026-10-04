use serde_json::{Value, json};
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
    fn execute(&self, command: &str, card: &Value) -> std::process::Output {
        let path = self.0.join("input.json");
        fs::write(&path, serde_json::to_vec(card).unwrap()).unwrap();
        Command::new(env!("CARGO_BIN_EXE_rustflow"))
            .arg(command)
            .arg(path)
            .output()
            .unwrap()
    }
    fn run(&self, command: &str, card: &Value) -> Value {
        let result = self.execute(command, card);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        serde_json::from_slice(&result.stdout).unwrap()
    }
    fn reject(&self, card: &Value, message: &str) {
        let result = self.execute("transport", card);
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr).contains(message),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
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
    assert!(result.get("continuation").is_none());
    assert!(result.get("identity").is_none());
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

fn prescribed_card() -> Value {
    let mut value: Value =
        serde_json::from_str(include_str!("../examples/cli/prescribed-transport.json")).unwrap();
    value["cache_directory"] = "bank".into();
    value
}
fn prescribed_coordinate() -> String {
    Atom::var(symbol!("rustflow_prescribed_cli::s")).to_canonical_string()
}
fn prescribed_root() -> String {
    Atom::var(symbol!("rustflow_prescribed_cli::r")).to_canonical_string()
}

#[test]
fn prescribed_root_cli_crosses_threshold_reuses_and_restarts_without_seeds() {
    let tmp = Temporary::new();
    let mut card = prescribed_card();
    let result = tmp.run("transport", &card);
    assert_eq!(result["continuation"]["kind"], "prescribed_affine");
    assert_eq!(
        result["continuation"]["prescriptions"][0]["polynomial"],
        prescribed_coordinate()
    );
    assert_eq!(result["continuation"]["prescriptions"][0]["side"], "-i0");
    assert_eq!(
        result["continuation"]["homotopy_admission"],
        "all_planner_routes_in_declared_domain"
    );
    let points = result["results"].as_array().unwrap();
    let p = Precision::decimal(80).unwrap();
    let expected = [
        [p.i(1), p.complex(-2, -2), p.complex(0, 4)],
        [p.i(1), p.complex(-2, -4), p.complex(-6, 8)],
        [p.i(1), p.complex(-2, -4), p.complex(-6, 8)],
    ];
    for (point, coefficients) in points.iter().zip(expected) {
        assert!(point["verified_digits"].as_u64().unwrap() >= 20);
        assert_eq!(point["root_germ"][prescribed_root()], "opposite");
        for (actual, wanted) in point["coefficients"]
            .as_array()
            .unwrap()
            .iter()
            .zip(coefficients)
        {
            assert!(p.close(&number(&actual[0], p), &wanted, 20));
        }
    }
    assert!(points[0]["steps"].as_u64().unwrap() > 0);
    assert!(points[1]["steps"].as_u64().unwrap() > 0);
    assert_eq!(points[1]["starting_point"][prescribed_coordinate()], "-1");
    assert_eq!(
        points[1]["starting_root_germ"][prescribed_root()],
        "opposite"
    );
    assert_eq!(points[2]["steps"], 0);
    assert!(result["retained_boundaries"].as_u64().unwrap() > 3);
    card.as_object_mut().unwrap().remove("seeds");
    let restart = tmp.run("transport", &card);
    assert_eq!(restart["identity"], result["identity"]);
    assert_eq!(restart["loaded_boundaries"], result["retained_boundaries"]);
    for point in restart["results"].as_array().unwrap() {
        assert_eq!(point["steps"], 0);
    }
}

fn rational_prescribed_card() -> Value {
    let mut card = prescribed_card();
    card.as_object_mut().unwrap().remove("roots");
    card["derivatives"] = json!({"s":[["eps/s"]]});
    card["seeds"][0]
        .as_object_mut()
        .unwrap()
        .remove("root_germ");
    card["destinations"] = json!([{"s":"-1"}]);
    card["branch_domain"] =
        "Log(s) on the explicitly prescribed route without extra winding".into();
    card
}

#[test]
fn prescribed_rational_and_rootless_canonical_cli_keep_native_representations() {
    let p = Precision::decimal(80).unwrap();
    let pi = Float::from_raw(rug::Float::with_val(p.bits, rug::float::Constant::Pi));
    let logarithm = ComplexFloat::new(p.real(0), -pi);
    let mut identities = Vec::new();
    for canonical in [false, true] {
        let tmp = Temporary::new();
        let mut card = rational_prescribed_card();
        if canonical {
            card.as_object_mut().unwrap().remove("derivatives");
            card["canonical"] = json!({"variables":["s"],"letters":["s"],"matrices":[[["1"]]]});
        }
        let result = tmp.run("transport", &card);
        let point = &result["results"][0];
        for (row, expected) in point["coefficients"].as_array().unwrap().iter().zip([
            p.i(1),
            logarithm.clone(),
            p.scale(&p.mul(&logarithm, &logarithm), 1, 2),
        ]) {
            assert!(p.close(&number(&row[0], p), &expected, 20));
        }
        assert!(point.get("root_germ").is_none());
        if canonical {
            assert_eq!(result["representation"], "canonical");
        }
        identities.push(result["identity"].as_str().unwrap().to_string());
        card.as_object_mut().unwrap().remove("seeds");
        let restart = tmp.run("transport", &card);
        assert_eq!(restart["results"][0]["steps"], 0);
    }
    assert_ne!(identities[0], identities[1]);
}

#[test]
fn prescribed_rootful_canonical_cli_tracks_the_declared_sheet() {
    let tmp = Temporary::new();
    let mut card = prescribed_card();
    card.as_object_mut().unwrap().remove("derivatives");
    card["canonical"] = json!({"variables":["s"],"letters":["r"],"matrices":[[["1"]]]});
    card["destinations"] = json!([{"coordinates":{"s":"-1"},"root_germ":{"r":"opposite"}}]);
    let result = tmp.run("transport", &card);
    assert_eq!(result["representation"], "canonical");
    let point = &result["results"][0];
    assert_eq!(point["root_germ"][prescribed_root()], "opposite");
    let p = Precision::decimal(80).unwrap();
    let logarithm = p.scale(
        &ComplexFloat::new(
            p.real(0),
            -Float::from_raw(rug::Float::with_val(p.bits, rug::float::Constant::Pi)),
        ),
        1,
        2,
    );
    assert!(p.close(&number(&point["coefficients"][1][0], p), &logarithm, 20));
    assert!(p.close(
        &number(&point["coefficients"][2][0], p),
        &p.scale(&p.mul(&logarithm, &logarithm), 1, 2),
        20
    ));
}

#[test]
fn prescribed_cli_rejects_missing_admission_invalid_sides_and_nonphysical_polynomials() {
    for (field, message) in [
        ("homotopy_admission", "homotopy_admission"),
        ("unprescribed_side", "unprescribed_side"),
    ] {
        let tmp = Temporary::new();
        let mut card = prescribed_card();
        card["continuation"].as_object_mut().unwrap().remove(field);
        tmp.reject(&card, message);
        assert!(!tmp.0.join("bank/physical-boundaries.bin").exists());
    }
    for (value, message) in [
        (json!("guessed"), "unknown variant"),
        (json!(true), "invalid type"),
    ] {
        let tmp = Temporary::new();
        let mut card = prescribed_card();
        card["continuation"]["homotopy_admission"] = value;
        tmp.reject(&card, message);
    }
    for side in ["upper", "principal", ""] {
        let tmp = Temporary::new();
        let mut card = prescribed_card();
        card["continuation"]["prescriptions"][0]["side"] = side.into();
        tmp.reject(&card, "unknown variant");
    }
    for polynomial in ["eps", "r", "foreign", "sqrt(s)", "1/s"] {
        let tmp = Temporary::new();
        let mut card = prescribed_card();
        card["continuation"]["prescriptions"][0]["polynomial"] = polynomial.into();
        let output = tmp.execute("transport", &card);
        assert!(!output.status.success(), "accepted {polynomial}");
        assert!(!tmp.0.join("bank/physical-boundaries.bin").exists());
    }
    let tmp = Temporary::new();
    let mut card = prescribed_card();
    card["continuation"]["unknown_option"] = true.into();
    tmp.reject(&card, "unknown field");
}

#[test]
fn prescribed_cli_wrong_germ_does_not_mutate_a_completed_bank() {
    let tmp = Temporary::new();
    let mut card = prescribed_card();
    card["destinations"].as_array_mut().unwrap().truncate(1);
    let good = tmp.run("transport", &card);
    let bank = tmp.0.join("bank/physical-boundaries.bin");
    let before = fs::read(&bank).unwrap();
    card.as_object_mut().unwrap().remove("seeds");
    card["destinations"][0]["root_germ"]["r"] = "principal".into();
    tmp.reject(&card, "germ");
    assert_eq!(fs::read(&bank).unwrap(), before);
    card["destinations"][0]["root_germ"]["r"] = "opposite".into();
    let repeat = tmp.run("transport", &card);
    assert_eq!(repeat["retained_boundaries"], good["retained_boundaries"]);
    assert_eq!(repeat["results"][0]["steps"], 0);
}

#[test]
fn prescribed_cli_side_policy_changes_identity_and_cannot_reuse_another_sheet() {
    let tmp = Temporary::new();
    let mut card = rational_prescribed_card();
    let lower = tmp.run("transport", &card);
    let saved_seeds = card.as_object_mut().unwrap().remove("seeds").unwrap();
    let bank = tmp.0.join("bank/physical-boundaries.bin");
    let before = fs::read(&bank).unwrap();
    card["continuation"]["prescriptions"][0]["side"] = "+i0".into();
    card["continuation"]["unprescribed_side"] = "+i0".into();
    tmp.reject(&card, "compatible");
    assert_eq!(fs::read(&bank).unwrap(), before);
    card["seeds"] = saved_seeds;
    let upper = tmp.run("transport", &card);
    assert_ne!(upper["identity"], lower["identity"]);
    assert!(upper["results"][0]["steps"].as_u64().unwrap() > 0);
    let p = Precision::decimal(80).unwrap();
    let upper_log = number(&upper["results"][0]["coefficients"][1][0], p);
    let lower_log = number(&lower["results"][0]["coefficients"][1][0], p);
    assert!(p.close(&upper_log, &p.scale(&lower_log, -1, 1), 20));
    assert!(
        upper["retained_boundaries"].as_u64().unwrap()
            > lower["retained_boundaries"].as_u64().unwrap()
    );
}
