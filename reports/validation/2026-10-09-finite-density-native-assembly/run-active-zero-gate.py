#!/usr/bin/env python3
"""Run frozen-source validation gates; never loads reference values."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("gate", choices=("native", "lib", "flow-boundary", "massless-sources",
    "source-fingerprint", "python", "two-loop-fixed", "two-loop-laurent",
    "three-loop-fixed", "three-loop-laurent", "three-loop-double"))
args = parser.parse_args()
snapshot = json.loads((BASE / "active-zero-source-hashes.json").read_text())
for name, expected in snapshot["files"].items():
    actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError("source changed after frozen build: " + name)
build = json.loads((BASE / "active-zero-build-provenance.json").read_text())
if build["source_snapshot_sha256"] != hashlib.sha256((BASE / "active-zero-source-hashes.json").read_bytes()).hexdigest():
    raise ValueError("build belongs to a different source snapshot")

prefixes = {"native": "rustred-", "lib": "symbolica_amflow-",
    "flow-boundary": "finite_density_flow_boundary-", "massless-sources": "finite_density_massless_sources-",
    "source-fingerprint": "source_fingerprint-", "python": "python_finite_density-"}
numeric = args.gate not in prefixes
prefix = "finite_density_runtime_flow-" if numeric else prefixes[args.gate]
matches = [name for name in build["files"] if Path(name).name.startswith(prefix)]
if len(matches) != 1:
    raise ValueError("build does not identify one executable for " + args.gate)
executable = matches[0]
if hashlib.sha256((ROOT / executable).read_bytes()).hexdigest() != build["files"][executable]["sha256"]:
    raise ValueError("executable changed after build capture")
tag = "active-zero-" + args.gate
if (BASE / (tag + "-resources.json")).exists() or (numeric and (BASE / tag).exists()):
    raise ValueError("gate artifacts already exist; preserve them and use a distinct run")
if numeric:
    three = args.gate.startswith("three-")
    laurent = args.gate.endswith("laurent")
    double = args.gate.endswith("double")
    settings = {
        "RUSTFLOW_WEIGHTED_INPUT": "examples/finite_density/massless_" + ("three_loop_chain" if three else "two_loop_sunset") + ".json",
        "RUSTFLOW_WEIGHTED_REQUESTED_RAYS": "true", "RUSTFLOW_WEIGHTED_FREE_VIRTUAL_ZEROS": "true",
        "RUSTFLOW_WEIGHTED_ZERO_FACES": "true", "RUSTFLOW_WEIGHTED_PRIORITIZE_REQUESTED": "false",
        "RUSTFLOW_WEIGHTED_REUSED_RULES": "65536", "RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL": "32",
        "RUSTFLOW_WEIGHTED_ACTIVE_TARGET_CLOSURE": "true", "RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS": "262144",
        "RUSTFLOW_WEIGHTED_GUARD_PASSES": "0", "RUSTFLOW_WEIGHTED_DEPTH": "3",
        "RUSTFLOW_WEIGHTED_DOMAINS": "8192", "RUSTFLOW_WEIGHTED_ROUNDS": "16",
        "RUSTFLOW_WEIGHTED_FRONTIER": "1024", "RUSTFLOW_WEIGHTED_REQUESTED": "4096",
        "RUSTFLOW_DENSITY_FLOW_EPSILON": "-5/4",
        "RUSTFLOW_DENSITY_FLOW_PROFILES": "18:60:8" if double else "18:60:8,28:60:8,28:80:8,28:80:12",
        "RUSTFLOW_DENSITY_FLOW_REPORT": str((BASE / tag).relative_to(ROOT)),
    }
    if double:
        settings["RUSTFLOW_WEIGHTED_CUTS"] = "0,3"
    test = "runtime_graph_occupied_flow" if double else "runtime_graph_full_laurent" if laurent else "runtime_graph_full_amplitude"
    command = ["timeout", "1800" if laurent and three else "600", "env",
        *(key + "=" + value for key, value in settings.items()), executable,
        "--ignored", "--exact", test, "--nocapture", "--test-threads=1"]
else:
    filters = {"native": ["solver::guarded"], "lib": ["finite_density::"]}
    command = ["timeout", "180", executable, *filters.get(args.gate, []), "--nocapture", "--test-threads=1"]
result = subprocess.run([sys.executable, str(BASE / "run-resource-command.py"),
    str(BASE / (tag + "-resources.json")), *command], cwd=ROOT)
sys.exit(result.returncode)
