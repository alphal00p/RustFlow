#!/usr/bin/env python3
"""Run one of five root regression gates against a captured frozen build.

No Cargo, numeric acceptance run, native-owner unit gate, or reference access.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIXES = {"lib": "symbolica_amflow-", "flow-boundary": "finite_density_flow_boundary-",
    "massless-sources": "finite_density_massless_sources-", "source-fingerprint": "source_fingerprint-",
    "python": "python_finite_density-"}

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("gate", choices=tuple(PREFIXES))
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9-]+", args.prefix):
        raise ValueError("unsafe report prefix")
    snapshot_path = BASE / (args.prefix + "-source-hashes.json")
    snapshot = json.loads(snapshot_path.read_text())
    for name, expected in snapshot["files"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("source changed after frozen build: " + name)
    build_path = BASE / (args.prefix + "-build-provenance.json")
    build = json.loads(build_path.read_text())
    if build["source_snapshot_sha256"] != digest(snapshot_path):
        raise ValueError("build belongs to a different source snapshot")
    matches = [name for name in build["files"] if Path(name).name.startswith(PREFIXES[args.gate])]
    if len(matches) != 1:
        raise ValueError("build does not identify exactly one gate executable")
    executable = matches[0]
    if digest(ROOT / executable) != build["files"][executable]["sha256"]:
        raise ValueError("executable changed after build capture")
    tag = args.prefix + "-" + args.gate
    resources = BASE / (tag + "-resources.json")
    if any(BASE.glob(tag + "-resources.*")):
        raise ValueError("gate artifacts already exist; preserve them and use a distinct prefix")
    filters = ["finite_density::"] if args.gate == "lib" else []
    command = ["timeout", "180", executable, *filters, "--nocapture", "--test-threads=1"]
    result = subprocess.run([sys.executable, str(BASE / "run-resource-command.py"), str(resources), *command], cwd=ROOT)
    evidence = {"source_snapshot_sha256": digest(snapshot_path), "build_provenance_sha256": digest(build_path),
        "launcher_sha256": digest(Path(__file__)), "gate": args.gate, "command": command,
        "exit_code": result.returncode, "native_numerical_acceptance_run": False}
    (BASE / (tag + "-binding.json")).write_text(json.dumps(evidence, indent=2) + "\n")
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
