#!/usr/bin/env python3
"""Fresh matched double-placement runs using immutable runtime inputs."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

BASE = Path(__file__).resolve().parent
mode = sys.argv[1]
assert mode in ["double-slots12", "double-slots24"]
manifest = json.loads((BASE / "capsule-manifest.json").read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
for item in manifest["files"]:
    p = BASE / item["capsule"]
    assert sha(p) == item["sha256"] and p.stat().st_size == item["bytes"]
out = BASE / mode
binding_path = BASE / f"{mode}-binding.json"
assert not out.exists() and not binding_path.exists()
env = {"PILOT_REQUESTS": "16384", "PILOT_ROUNDS": "24", "PILOT_FRONTIER": "4096",
       "PILOT_RULES": "65536", "PILOT_ZERO_ATTEMPTS": "8192", "PILOT_RAY_DOMAINS": "1",
       "PILOT_POINT_DOMAINS": "1", "PILOT_DEPTH": "3"}
command = ["timeout", "1200", str(BASE / "capsule/pilot"), str(BASE / "capsule/input.json"),
           str(out), mode, manifest["proof_binding_sha256"]]
report = {"scope": "Fresh isolated algebraic closure pilot; no historical program import, no numerical admission.",
          "command": command, "environment": env, "capsule_manifest_sha256": sha(BASE / "capsule-manifest.json"),
          "launcher_sha256": sha(Path(__file__)), "build_binding_sha256": sha(BASE / "capsule/build-binding.json"),
          "unchanged_algorithm_and_sources": True, "runtime_reads_shared_libraries": False}
binding_path.write_text(json.dumps(report, indent=2) + "\n")
result = subprocess.run([sys.executable, str(BASE / "capsule/run-resource-command.py"),
    str(BASE / f"{mode}-resources.json"), *command], cwd=BASE.parents[2], env={**os.environ, **env})
report["exit_code"] = result.returncode
report["capsule_unchanged"] = all(sha(BASE / i["capsule"]) == i["sha256"] for i in manifest["files"])
assert report["capsule_unchanged"]
binding_path.write_text(json.dumps(report, indent=2) + "\n")
raise SystemExit(result.returncode)
