#[path = "../build_support/source_fingerprint.rs"]
mod source_fingerprint;

use serde_json::json;
use source_fingerprint::{hash_package, resolved_packages};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let root = std::env::temp_dir().join(format!(
            "rustflow-source-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(root.join("src")).unwrap();
        fs::write(
            root.join("Cargo.toml"),
            "[package]\nname='fixture'\nversion='0.1.0'\n",
        )
        .unwrap();
        Self(root)
    }
    fn package(&self) -> serde_json::Value {
        json!({"name":"fixture", "version":"0.1.0", "source":null,"manifest_path":self.0.join("Cargo.toml")})
    }
    fn digest(&self) -> blake3::Hash {
        let mut hash = blake3::Hasher::new();
        hash_package(&self.package(), &mut hash);
        hash.finalize()
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

#[test]
fn fingerprint_tracks_dirty_sources_and_additions_but_not_build_products() {
    let source = Scratch::new();
    fs::write(source.0.join("src/lib.rs"), "pub fn identity() {}\n").unwrap();
    let initial = source.digest();
    fs::create_dir_all(source.0.join("target")).unwrap();
    fs::write(source.0.join("target/output"), "not source").unwrap();
    assert_eq!(initial, source.digest());
    fs::write(source.0.join("src/added.rs"), "// new module\n").unwrap();
    let added = source.digest();
    assert_ne!(initial, added);
    fs::write(source.0.join("src/lib.rs"), "pub fn changed() {}\n").unwrap();
    assert_ne!(added, source.digest());
}

#[test]
fn package_fingerprint_is_independent_of_checkout_location() {
    let first = Scratch::new();
    let second = Scratch::new();
    fs::write(first.0.join("src/lib.rs"), "// identical\n").unwrap();
    fs::write(second.0.join("src/lib.rs"), "// identical\n").unwrap();
    assert_eq!(first.digest(), second.digest());
}

#[test]
fn source_graph_follows_resolved_dependency_ids_and_rejects_wrong_host() {
    let own = Scratch::new();
    let dep = Scratch::new();
    let unrelated = Scratch::new();
    let mut own_package = own.package();
    own_package["id"] = json!("own");
    let mut dep_package = dep.package();
    dep_package["id"] = json!("dep");
    let mut unrelated_package = unrelated.package();
    unrelated_package["id"] = json!("unrelated");
    let metadata = json!({"packages":[own_package,dep_package,unrelated_package],"resolve":{"nodes":[
        {"id":"own","dependencies":["dep"]},{"id":"dep","dependencies":[]},{"id":"unrelated","dependencies":[]}
    ]}});
    let selected = resolved_packages(&metadata, &own.0.join("Cargo.toml"));
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0]["id"], "dep");
    let wrong = Scratch::new();
    assert!(
        std::panic::catch_unwind(|| resolved_packages(&metadata, &wrong.0.join("Cargo.toml")))
            .is_err()
    );
}

#[test]
fn mismatched_python_features_fail_instead_of_hashing_the_default_host_graph() {
    use std::collections::BTreeSet;
    let own = Scratch::new();
    let mut package = own.package();
    package["id"] = json!("own");
    package["features"] = json!({"default":[],"python":[],"python_stubgen":[]});
    let mut metadata = json!({"packages":[package], "resolve":{"nodes":[
        {"id":"own", "dependencies":[], "features":["python"]}
    ]}});
    let actual = BTreeSet::from(["python", "python_stubgen"]);
    assert!(
        std::panic::catch_unwind(|| source_fingerprint::validate_features(
            &metadata,
            &own.0.join("Cargo.toml"),
            &actual
        ))
        .is_err()
    );
    metadata["resolve"]["nodes"][0]["features"] = json!(["python", "python_stubgen"]);
    source_fingerprint::validate_features(&metadata, &own.0.join("Cargo.toml"), &actual);
}

#[test]
fn embedded_fixture_and_cargo_home_override_changes_affect_the_digest() {
    use source_fingerprint::{hash_file, hash_tree};
    let own = Scratch::new();
    let fixture = own.0.join("fixtures/gg-hg");
    fs::create_dir_all(&fixture).unwrap();
    fs::write(fixture.join("published-test.json"), "{\"value\":1}").unwrap();
    let home = own.0.join("cargo-home");
    fs::create_dir_all(&home).unwrap();
    fs::write(home.join("config.toml"), "[patch.crates-io]\n").unwrap();
    let digest = || {
        let mut h = blake3::Hasher::new();
        hash_tree(&fixture, &own.0, &mut h);
        hash_file(
            &home.join("config.toml"),
            std::path::Path::new("cargo-home/config.toml"),
            &mut h,
        );
        h.finalize()
    };
    let initial = digest();
    fs::write(fixture.join("published-test.json"), "{\"value\":2}").unwrap();
    let embedded_change = digest();
    assert_ne!(initial, embedded_change);
    fs::write(
        home.join("config.toml"),
        "[patch.crates-io]\n# different override\n",
    )
    .unwrap();
    assert_ne!(embedded_change, digest());
}
