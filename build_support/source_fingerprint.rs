//! Source-sensitive build identities for standalone and embedded workspaces.
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

/// Ignore build products and repository administration, never source files.
fn excluded(name: &std::ffi::OsStr) -> bool {
    [".git", "target", ".venv", "venv", "__pycache__", ".cache"]
        .iter()
        .any(|n| name == *n)
}

fn field(hasher: &mut blake3::Hasher, bytes: &[u8]) {
    hasher.update(&(bytes.len() as u64).to_le_bytes());
    hasher.update(bytes);
}

pub fn hash_file(path: &Path, label: &Path, hasher: &mut blake3::Hasher) {
    println!("cargo:rerun-if-changed={}", path.display());
    field(hasher, label.to_string_lossy().as_bytes());
    field(
        hasher,
        &fs::read(path).expect("read source fingerprint input"),
    );
}

pub fn hash_tree(path: &Path, root: &Path, hasher: &mut blake3::Hasher) {
    // Watching directories also invalidates the fingerprint on newly added files.
    println!("cargo:rerun-if-changed={}", path.display());
    let mut entries = fs::read_dir(path)
        .expect("read dependency source tree")
        .map(|e| e.expect("read source directory entry").path())
        .collect::<Vec<_>>();
    entries.sort();
    for entry in entries {
        if excluded(entry.file_name().expect("source filename")) {
            continue;
        }
        if entry.is_dir() {
            hash_tree(&entry, root, hasher);
        } else {
            hash_file(&entry, entry.strip_prefix(root).unwrap(), hasher);
        }
    }
}

/// Select the actual resolved dependency graph, rather than guessing sibling paths.
pub fn resolved_packages<'a>(
    metadata: &'a serde_json::Value,
    own_manifest: &Path,
) -> Vec<&'a serde_json::Value> {
    let packages = metadata["packages"].as_array().expect("Cargo package list");
    let own = own_package(metadata, own_manifest);
    let nodes = metadata["resolve"]["nodes"]
        .as_array()
        .expect("Cargo resolved graph");
    let graph = nodes
        .iter()
        .map(|n| (n["id"].as_str().unwrap(), n))
        .collect::<BTreeMap<_, _>>();
    let own_id = own["id"].as_str().unwrap();
    let mut visited = BTreeSet::new();
    let mut pending = vec![own_id];
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        let node = graph
            .get(id)
            .expect("RustFlow missing from active host dependency graph");
        for dep in node["dependencies"].as_array().unwrap() {
            pending.push(dep.as_str().unwrap());
        }
    }
    let mut selected = packages
        .iter()
        .filter(|p| {
            let id = p["id"].as_str().unwrap();
            visited.contains(id) && id != own_id
        })
        .collect::<Vec<_>>();
    selected.sort_by_key(|p| {
        (
            p["name"].as_str().unwrap(),
            p["version"].as_str().unwrap(),
            p["source"].as_str().unwrap_or("path"),
        )
    });
    selected
}

/// Hash package sources and ancestor workspace manifests, without absolute paths.
pub fn hash_package(package: &serde_json::Value, hasher: &mut blake3::Hasher) {
    field(hasher, package["name"].as_str().unwrap().as_bytes());
    field(hasher, package["version"].as_str().unwrap().as_bytes());
    field(
        hasher,
        package["source"].as_str().unwrap_or("path").as_bytes(),
    );
    let manifest = Path::new(package["manifest_path"].as_str().unwrap());
    let root = manifest.parent().unwrap();
    // Registry archives are immutable and their checksums are in the host lock.
    // Git and path sources can be edited locally, so hash their actual contents.
    if package["source"]
        .as_str()
        .is_some_and(|s| s.starts_with("registry+"))
    {
        return;
    }
    hash_tree(root, root, hasher);
    for (depth, parent) in root.ancestors().skip(1).enumerate() {
        let ancestor = parent.join("Cargo.toml");
        if ancestor.is_file() {
            hash_file(
                &ancestor,
                &PathBuf::from(format!("ancestor-{depth}/Cargo.toml")),
                hasher,
            );
        }
        // Worktrees have a .git file; ordinary checkouts have a .git directory.
        if parent.join(".git").exists() {
            break;
        }
    }
}

/// Locate this package even when the host uses a dependency alias.
pub fn own_package<'a>(
    metadata: &'a serde_json::Value,
    own_manifest: &Path,
) -> &'a serde_json::Value {
    let packages = metadata["packages"].as_array().expect("Cargo package list");
    let own_manifest = own_manifest
        .canonicalize()
        .expect("canonical crate manifest");
    packages
        .iter()
        .find(|p| {
            Path::new(p["manifest_path"].as_str().unwrap())
                .canonicalize()
                .is_ok_and(|path| path == own_manifest)
        })
        .expect("RUSTFLOW_WORKSPACE_MANIFEST does not resolve this RustFlow checkout")
}

pub fn resolved_features<'a>(
    metadata: &'a serde_json::Value,
    own_manifest: &Path,
) -> BTreeSet<&'a str> {
    let own = own_package(metadata, own_manifest);
    metadata["resolve"]["nodes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|node| node["id"] == own["id"])
        .expect("current package absent from resolved graph")["features"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f.as_str().unwrap())
        .collect()
}

pub fn validate_features(
    metadata: &serde_json::Value,
    own_manifest: &Path,
    actual: &BTreeSet<&str>,
) {
    assert_eq!(
        resolved_features(metadata, own_manifest),
        *actual,
        "Cargo metadata features differ from this build. Set RUSTFLOW_WORKSPACE_FEATURES and RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES to match the embedding Cargo invocation"
    );
}
