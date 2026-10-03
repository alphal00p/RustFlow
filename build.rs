use std::{fs, path::Path};

fn hash_tree(path: &Path, root: &Path, hasher: &mut blake3::Hasher) {
    let mut entries = fs::read_dir(path)
        .expect("read RustRed source tree")
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
            hasher.update(&fs::read(&entry).expect("read RustRed source"));
        }
    }
}
fn main() {
    let mut dependencies = blake3::Hasher::new();
    for path in ["Cargo.toml", "Cargo.lock"] {
        println!("cargo:rerun-if-changed={path}");
        dependencies.update(&fs::read(path).expect("read dependency fingerprint input"));
    }
    println!(
        "cargo:rustc-env=DEPENDENCY_SOURCE_DIGEST={}",
        dependencies.finalize()
    );
    let mut own = blake3::Hasher::new();
    hash_tree(Path::new("src"), Path::new(""), &mut own);
    println!("cargo:rustc-env=PORT_SOURCE_DIGEST={}", own.finalize());
    let mut hasher = blake3::Hasher::new();
    hash_tree(
        Path::new("../rustred/crates/rustred-core/src"),
        Path::new("../rustred/crates/rustred-core"),
        &mut hasher,
    );
    println!("cargo:rerun-if-changed=../rustred/crates/rustred-core/Cargo.toml");
    hasher.update(&fs::read("../rustred/crates/rustred-core/Cargo.toml").unwrap());
    println!(
        "cargo:rustc-env=RUSTRED_SOURCE_DIGEST={}",
        hasher.finalize()
    );
}
