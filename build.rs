#[path = "build_support/source_fingerprint.rs"]
mod source_fingerprint;

use source_fingerprint::{hash_file, hash_package, hash_tree, resolved_packages};
use std::{env, fs, path::PathBuf, process::Command};

fn main() {
    println!("cargo:rerun-if-env-changed=RUSTFLOW_WORKSPACE_MANIFEST");
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
    if host_manifest == own_manifest.canonicalize().unwrap() {
        if env::var_os("CARGO_FEATURE_PYTHON_STUBGEN").is_some() {
            command.args(["--features", "python_stubgen"]);
        } else if env::var_os("CARGO_FEATURE_PYTHON").is_some() {
            command.args(["--features", "python"]);
        }
    }
    let output = command.output().expect("run Cargo dependency metadata");
    assert!(
        output.status.success(),
        "Cannot fingerprint the resolved dependency sources: {}. Fetch dependencies first; embedded builds must set RUSTFLOW_WORKSPACE_MANIFEST to the owning workspace manifest.",
        String::from_utf8_lossy(&output.stderr)
    );
    let metadata: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("parse Cargo dependency metadata");
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
    let mut rustred = blake3::Hasher::new();
    rustred.update(b"rustflow-resolved-rustred-v2");
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
