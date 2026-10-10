#!/usr/bin/env python3
"""Run one assigned regression gate against a captured frozen build.

No Cargo, numerical acceptance run, or reference access. Coordinate ownership before running native.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIXES = {"lib": "symbolica_amflow-", "flow-boundary": "finite_density_flow_boundary-",
    "massless-sources": "finite_density_massless_sources-", "source-fingerprint": "source_fingerprint-",
    "python": "python_finite_density-", "capacity": "finite_density_capacity-", "polynomial-sources": "finite_density_polynomial_sources-", "native": "rustred-"}

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
    from frozen_sources import verify_snapshot
    snapshot = verify_snapshot(snapshot_path)
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
    filters = ["finite_density::"] if args.gate == "lib" else ["solver::guarded"] if args.gate == "native" else []
    command = ["timeout", "180", executable, *filters, "--nocapture", "--test-threads=1"]
    environment = dict(os.environ)
    source_report = BASE / (args.prefix + "-actual62-source-equivalence.json")
    if args.gate == "polynomial-sources":
        if source_report.exists():
            raise ValueError("preserve existing exact source report")
        environment["RUSTFLOW_POLYNOMIAL_SOURCE_REPORT"] = str(source_report)
    result = subprocess.run([sys.executable, str(BASE / "run-resource-command.py"), str(resources), *command], cwd=ROOT, env=environment)
    verify_snapshot(snapshot_path)
    if digest(ROOT / executable) != build["files"][executable]["sha256"]:
        raise ValueError("executable changed during gate execution")
    evidence = {"source_snapshot_sha256": digest(snapshot_path), "build_provenance_sha256": digest(build_path),
        "launcher_sha256": digest(Path(__file__)), "gate": args.gate, "command": command,
        "exit_code": result.returncode, "native_numerical_acceptance_run": False}
    if args.gate == "polynomial-sources":
        evidence["report_environment"] = {"RUSTFLOW_POLYNOMIAL_SOURCE_REPORT": str(source_report)}
        if result.returncode == 0:
            if not source_report.is_file():
                raise ValueError("passing equivalence gate did not emit its exact report")
            evidence["source_equivalence_report"] = {"path": str(source_report.relative_to(ROOT)), "sha256": digest(source_report)}
    (BASE / (tag + "-binding.json")).write_text(json.dumps(evidence, indent=2) + "\n")
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
