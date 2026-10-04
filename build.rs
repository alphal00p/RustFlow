use std::{fs, path::Path};

fn hash_tree(path: &Path, root: &Path, hasher: &mut blake3::Hasher) {
    let mut entries = fs::read_dir(path)
        .expect("read dependency source tree")
        .map(|e| e.unwrap().path())
        .collect::<Vec<_>>();
    entries.sort();
    for entry in entries {
        if entry.is_dir() {
            hash_tree(&entry, root, hasher);
        } else {
            println!("cargo:rerun-if-changed={}", entry.display());
            hasher.update(
                entry
                    .strip_prefix(root)
                    .unwrap()
                    .to_string_lossy()
                    .as_bytes(),
            );
            hasher.update(&fs::read(&entry).expect("read dependency source"));
        }
    }
}
fn main() {
    let mut dependencies = blake3::Hasher::new();
    for path in ["Cargo.toml", "Cargo.lock"] {
        println!("cargo:rerun-if-changed={path}");
        dependencies.update(&fs::read(path).expect("read dependency fingerprint input"));
    }
    // Path dependencies carry no source revision in Cargo.lock. Include the
    // consumed native HEPKit implementation and workspace configuration in
    // symbolic/numerical cache compatibility, including local fixes.
    let hepkit = Path::new("../hepkit");
    let manifest = hepkit.join("Cargo.toml");
    println!("cargo:rerun-if-changed={}", manifest.display());
    dependencies.update(&fs::read(&manifest).expect("read HEPKit workspace manifest"));
    for name in [
        "feynkit-graph",
        "feynkit-model",
        "feynkit-kinematics",
        "feynkit-tensor",
        "linnet",
        "spenso",
        "idenso",
        "spenso-macros",
        "symbolica-utils",
    ] {
        dependencies.update(name.as_bytes());
        let package = hepkit.join("crates").join(name);
        let manifest = package.join("Cargo.toml");
        println!("cargo:rerun-if-changed={}", manifest.display());
        dependencies.update(&fs::read(&manifest).expect("read HEPKit package manifest"));
        hash_tree(&package.join("src"), &package, &mut dependencies);
        if name == "feynkit-model" {
            hash_tree(&package.join("data"), &package, &mut dependencies);
        }
    }
    println!(
        "cargo:rustc-env=DEPENDENCY_SOURCE_DIGEST={}",
        dependencies.finalize()
    );
    let mut own = blake3::Hasher::new();
    hash_tree(Path::new("src"), Path::new(""), &mut own);
    println!("cargo:rustc-env=PORT_SOURCE_DIGEST={}", own.finalize());
    let mut hasher = blake3::Hasher::new();
    // The core delegates integral ordering to a sibling crate. Both affect
    // reduction replay, and workspace settings supply inherited package data.
    let rustred = Path::new("../rustred");
    let manifest = rustred.join("Cargo.toml");
    println!("cargo:rerun-if-changed={}", manifest.display());
    hasher.update(&fs::read(&manifest).expect("read RustRed workspace manifest"));
    for name in ["rustred-core", "rustred-order"] {
        hasher.update(name.as_bytes());
        let package = rustred.join("crates").join(name);
        let manifest = package.join("Cargo.toml");
        println!("cargo:rerun-if-changed={}", manifest.display());
        hasher.update(&fs::read(&manifest).expect("read RustRed package manifest"));
        let build_script = package.join("build.rs");
        if build_script.exists() {
            println!("cargo:rerun-if-changed={}", build_script.display());
            hasher.update(&fs::read(&build_script).expect("read RustRed build script"));
        }
        hash_tree(&package.join("src"), &package, &mut hasher);
    }
    println!(
        "cargo:rustc-env=RUSTRED_SOURCE_DIGEST={}",
        hasher.finalize()
    );
}
