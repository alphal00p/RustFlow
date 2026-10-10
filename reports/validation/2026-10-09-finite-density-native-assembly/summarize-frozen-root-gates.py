#!/usr/bin/env python3
"""Summarize the five completed root regression gates, without executing tests."""
import argparse
import hashlib
import json
from pathlib import Path
import re

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
GATES = ("lib", "flow-boundary", "massless-sources", "source-fingerprint", "python")

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def item(path):
    return {"path": str(path.relative_to(ROOT)), "sha256": digest(path)}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9-]+", args.prefix):
        raise ValueError("unsafe report prefix")
    build_path = BASE / (args.prefix + "-build-provenance.json")
    build = json.loads(build_path.read_text())
    snapshot_path = BASE / build["source_snapshot"]
    snapshot = json.loads(snapshot_path.read_text())
    for path, expected in snapshot["files"].items():
        if digest(ROOT / path) != expected:
            raise ValueError("source changed after gate execution: " + path)
    records = []
    for gate in GATES:
        stem = args.prefix + "-" + gate
        paths = {kind: BASE / (stem + suffix) for kind, suffix in {
            "resources": "-resources.json", "log": "-resources.log",
            "launch": "-resources.provenance.json", "binding": "-binding.json"}.items()}
        resources = json.loads(paths["resources"].read_text())
        binding = json.loads(paths["binding"].read_text())
        if resources["exit_code"] != 0 or binding["exit_code"] != 0:
            raise ValueError("gate failed: " + gate)
        if binding["source_snapshot_sha256"] != digest(snapshot_path) or binding["build_provenance_sha256"] != digest(build_path):
            raise ValueError("gate does not bind this build: " + gate)
        matches = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out; finished in ([0-9.]+)s", paths["log"].read_text())
        if len(matches) != 1 or int(matches[0][1]) != 0:
            raise ValueError("unexpected harness result: " + gate)
        passed, failed, ignored, measured, filtered = map(int, matches[0][:5])
        records.append({"gate": gate, "status": "passed", "passed": passed, "failed": failed,
            "ignored": ignored, "measured": measured, "filtered": filtered,
            "harness_seconds": float(matches[0][5]),
            "wall_seconds": resources["wall_seconds"], "peak_child_rss_kib": resources["peak_child_rss_kib"],
            "artifacts": {kind: item(path) for kind, path in paths.items()}})
    report = {"schema_version": 1, "status": "passed", "scope": "Five root regression gates only; excludes the separately owned native guarded suite and all numerical acceptance runs.",
        "source_snapshot": item(snapshot_path), "build_provenance": item(build_path),
        "unchanged_source_files": len(snapshot["files"]), "compile_identities": build["compile_identities"],
        "features": build["features"], "gates": records,
        "unique_test_count": sum(row["passed"] for row in records),
        "numerical_acceptance_runs": 0, "native_four_loop_acceptance": False,
        "three_loop_full_amplitude_validated": False,
        "timing_scope": "Prebuilt test execution; compilation and Nix startup excluded. Shared-host measurements are not isolated performance benchmarks.",
        "summary_script_sha256": digest(Path(__file__))}
    output = BASE / (args.prefix + "-root-gates.json")
    if output.exists():
        raise ValueError("preserve existing gate summary")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "passed", "gates": len(records), "tests": report["unique_test_count"], "report": str(output.relative_to(ROOT))}))

if __name__ == "__main__":
    main()
