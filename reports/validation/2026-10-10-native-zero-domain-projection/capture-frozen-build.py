#!/usr/bin/env python3
"""Capture an already completed shared Cargo build; never invokes Cargo."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]

def read(path):
    return json.loads(path.read_text())

def digest(path):
    before = path.stat()
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    after = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
        raise ValueError("artifact changed during capture: " + str(path))
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--template", default="../2026-10-09-finite-density-native-assembly/active-zero-canonical-build-provenance.json")
    parser.add_argument("--extra-artifact", action="append", default=[], help="Additional exact executable path from the completed Cargo build log")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9-]+", args.prefix):
        raise ValueError("unsafe report prefix")
    snapshot_path = BASE / (args.prefix + "-source-hashes.json")
    output = BASE / (args.prefix + "-build-provenance.json")
    if output.exists():
        raise ValueError("preserve the existing build capture; choose a distinct prefix")
    resource_path = BASE / (args.prefix + "-build-resources.json")
    log_path = BASE / (args.prefix + "-build-resources.log")
    launch_path = BASE / (args.prefix + "-build-resources.provenance.json")
    resources = read(resource_path)
    if resources["exit_code"] != 0 or not re.search(r"Finished .*release", log_path.read_text()):
        raise ValueError("the shared release build has not completed successfully")
    from frozen_sources import verify_snapshot
    snapshot = verify_snapshot(snapshot_path)
    template = read(BASE / args.template)
    names = list(template["files"])
    if len(names) != 10:
        raise ValueError("expected the ten baseline shared build artifacts")
    names.extend(args.extra_artifact)
    if len(names) != len(set(names)):
        raise ValueError("duplicate build artifact")
    emitted = re.findall(r"^\s*Executable .+ \(([^()]+)\)$", log_path.read_text(), re.MULTILINE)
    captured_executables = [name for name in names if "/deps/" in name and not name.endswith(".rlib")]
    if len(emitted) != len(set(emitted)) or set(emitted) != set(captured_executables):
        raise ValueError("capture does not match the completed build executable inventory")
    files = {}
    for name in names:
        path = ROOT / name
        sha = digest(path)
        stat = path.stat()
        files[name] = {"sha256": sha, "size_bytes": stat.st_size,
            "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()}
    feature_paths = {
        "rustred_unit": "target/release/.fingerprint/rustred-258a1e04558e72f1/test-lib-rustred.json",
        "rustred_library": "target/release/.fingerprint/rustred-392075d0c6995fd1/lib-rustred.json",
        "rustflow": "target/release/.fingerprint/symbolica-amflow-6c032eef1649ac2e/test-lib-symbolica_amflow.json",
    }
    features = {key: json.loads(read(ROOT / path)["features"]) for key, path in feature_paths.items()}
    if features != template["features"]:
        raise ValueError("shared build features changed: " + repr(features))
    environment = ROOT / next(name for name in names if name.endswith("/output"))
    identities = {}
    for key in ("DEPENDENCY_SOURCE_DIGEST", "RUSTRED_SOURCE_DIGEST", "PORT_SOURCE_DIGEST"):
        values = re.findall(r"cargo:rustc-env=" + key + r"=([0-9a-f]{64})", environment.read_text())
        if len(values) != 1:
            raise ValueError("missing/ambiguous compile identity: " + key)
        identities[key] = values[0]
    verify_snapshot(snapshot_path)
    report = {
        "source_files_verified_unchanged": len(snapshot["files"]),
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "source_snapshot": snapshot_path.name,
        "source_snapshot_sha256": digest(snapshot_path),
        "features": features,
        "feature_fingerprints": {key: {"path": path, "sha256": digest(ROOT / path)} for key, path in feature_paths.items()},
        "checks": {"combined_release_build": "passed", "wall_seconds": resources["wall_seconds"],
            "measurement_scope": "shared host; build times are not isolated benchmarks; regression execution recorded separately"},
        "build_resources": {"path": str(resource_path.relative_to(ROOT)), "sha256": digest(resource_path)},
        "build_log": {"path": str(log_path.relative_to(ROOT)), "sha256": digest(log_path)},
        "build_launch": {"path": str(launch_path.relative_to(ROOT)), "sha256": digest(launch_path)},
        "files": files,
        "compile_identities": identities,
        "capture_script_sha256": digest(Path(__file__)),
    }
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "captured", "source_files": len(snapshot["files"]), "artifacts": len(files), "report": str(output.relative_to(ROOT))}))

if __name__ == "__main__":
    main()
