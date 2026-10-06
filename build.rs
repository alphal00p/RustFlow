#[path = "build_support/source_fingerprint.rs"]
mod source_fingerprint;

use source_fingerprint::{
    hash_file, hash_package, hash_tree, own_package, resolved_packages, validate_features,
};
use std::{env, fs, path::PathBuf, process::Command};

fn main() {
    for variable in [
        "RUSTFLOW_WORKSPACE_MANIFEST",
        "RUSTFLOW_WORKSPACE_FEATURES",
        "RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES",
        "CARGO_HOME",
        "RUSTRED_RUNTIME_ARITIES",
    ] {
        println!("cargo:rerun-if-env-changed={variable}");
    }
    let own_root = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let own_manifest = own_root.join("Cargo.toml");
    let host_manifest = env::var_os("RUSTFLOW_WORKSPACE_MANIFEST")
        .map(PathBuf::from)
        .unwrap_or_else(|| own_manifest.clone())
        .canonicalize()
        .expect("RUSTFLOW_WORKSPACE_MANIFEST must name an existing Cargo.toml");
    let mut command = Command::new(env::var_os("CARGO").unwrap());
    command
        .current_dir(host_manifest.parent().unwrap())
        .args([
            "metadata",
            "--offline",
            "--locked",
            "--format-version",
            "1",
            "--manifest-path",
        ])
        .arg(&host_manifest)
        .args(["--filter-platform", &env::var("TARGET").unwrap()]);
    // The standalone invocation must also resolve enabled optional Python crates.
    // In a host, its dependency declaration supplies these features.
    let standalone = host_manifest == own_manifest.canonicalize().unwrap();
    if standalone {
        let features = [
            "native",
            "automatic",
            "wasm",
            "python_api",
            "python",
            "python_wasm",
            "python_stubgen",
        ]
        .into_iter()
        .filter(|feature| {
            env::var_os(format!("CARGO_FEATURE_{}", feature.to_uppercase())).is_some()
        })
        .collect::<Vec<_>>();
        if !features.is_empty() {
            command.args(["--features", &features.join(",")]);
        }
    }
    if let Ok(features) = env::var("RUSTFLOW_WORKSPACE_FEATURES")
        && !features.is_empty()
    {
        command.args(["--features", &features]);
    }
    let no_defaults = match env::var("RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES").as_deref() {
        Ok("1") => true,
        Ok("0" | "") | Err(_) => false,
        Ok(_) => panic!("RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES must be 0 or 1"),
    };
    if no_defaults || (standalone && env::var_os("CARGO_FEATURE_DEFAULT").is_none()) {
        command.arg("--no-default-features");
    }
    let output = command.output().expect("run Cargo dependency metadata");
    assert!(
        output.status.success(),
        "Cannot fingerprint the resolved dependency sources: {}. Fetch dependencies first; embedded builds must set RUSTFLOW_WORKSPACE_MANIFEST to the owning workspace manifest.",
        String::from_utf8_lossy(&output.stderr)
    );
    let metadata: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("parse Cargo dependency metadata");
    let actual_features = own_package(&metadata, &own_manifest)["features"]
        .as_object()
        .unwrap()
        .keys()
        .filter(|feature| {
            env::var_os(format!(
                "CARGO_FEATURE_{}",
                feature.to_uppercase().replace('-', "_")
            ))
            .is_some()
        })
        .map(String::as_str)
        .collect();
    validate_features(&metadata, &own_manifest, &actual_features);
    let workspace = PathBuf::from(metadata["workspace_root"].as_str().unwrap());
    let mut dependencies = blake3::Hasher::new();
    dependencies.update(b"rustflow-resolved-sources-v2");
    for path in [workspace.join("Cargo.toml"), workspace.join("Cargo.lock")] {
        hash_file(
            &path,
            std::path::Path::new(path.file_name().unwrap()),
            &mut dependencies,
        );
    }
    // Track Cargo configuration too, including package source overrides.
    for ancestor in workspace.ancestors() {
        for name in [".cargo/config.toml", ".cargo/config"] {
            let config = ancestor.join(name);
            if config.is_file() {
                hash_file(&config, std::path::Path::new(name), &mut dependencies);
            }
        }
    }
    let cargo_home = env::var_os("CARGO_HOME")
        .map(PathBuf::from)
        .or_else(|| env::var_os("HOME").map(|home| PathBuf::from(home).join(".cargo")));
    if let Some(home) = cargo_home {
        for name in ["config.toml", "config"] {
            let path = home.join(name);
            if path.is_file() {
                hash_file(
                    &path,
                    &PathBuf::from(format!("cargo-home/{name}")),
                    &mut dependencies,
                );
            }
        }
    }
    dependencies.update(serde_json::to_string(&actual_features).unwrap().as_bytes());
    let mut rustred = blake3::Hasher::new();
    rustred.update(b"rustflow-resolved-rustred-v2");
    // RustRed's existing build-time registry setting affects available backends.
    let arities = env::var("RUSTRED_RUNTIME_ARITIES").ok();
    let arities = serde_json::to_vec(&arities).unwrap();
    dependencies.update(b"RUSTRED_RUNTIME_ARITIES");
    dependencies.update(&arities);
    rustred.update(b"RUSTRED_RUNTIME_ARITIES");
    rustred.update(&arities);
    for package in resolved_packages(&metadata, &own_manifest) {
        hash_package(package, &mut dependencies);
        let node = metadata["resolve"]["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|node| node["id"] == package["id"])
            .unwrap();
        dependencies.update(serde_json::to_string(&node["features"]).unwrap().as_bytes());
        if matches!(package["name"].as_str(), Some("rustred" | "rustred-order")) {
            hash_package(package, &mut rustred);
            rustred.update(serde_json::to_string(&node["features"]).unwrap().as_bytes());
        }
    }
    println!(
        "cargo:rustc-env=DEPENDENCY_SOURCE_DIGEST={}",
        dependencies.finalize()
    );
    println!(
        "cargo:rustc-env=RUSTRED_SOURCE_DIGEST={}",
        rustred.finalize()
    );
    let mut own = blake3::Hasher::new();
    hash_tree(&own_root.join("src"), &own_root, &mut own);
    hash_tree(&own_root.join("build_support"), &own_root, &mut own);
    // These files are embedded by the native amplitude module with include_str!.
    hash_tree(&own_root.join("fixtures/gg-hg"), &own_root, &mut own);
    hash_file(
        &own_root.join("build.rs"),
        std::path::Path::new("build.rs"),
        &mut own,
    );
    hash_file(&own_manifest, std::path::Path::new("Cargo.toml"), &mut own);
    println!("cargo:rustc-env=PORT_SOURCE_DIGEST={}", own.finalize());
    // Retain the resolved graph for inspection without exposing environment values.
    fs::write(
        PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("dependency-metadata.json"),
        output.stdout,
    )
    .expect("write dependency provenance");
}
