#!/usr/bin/env python3
"""Build an isolated native ordering variant against the bound dependency graph."""
import difflib
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREVIOUS = json.loads((BASE / "build-binding.json").read_text())
REPORT = BASE / "cut-then-degree-build"
ISOLATED = Path("/tmp/rustflow-cut-then-occupation-degree-core-20261010")
OUTPUT = Path("/tmp/librustred-cut-then-occupation-degree-20261010.rlib")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert not REPORT.exists() and not ISOLATED.exists() and not OUTPUT.exists()
for name, expected in PREVIOUS["base_native_sources"].items():
    assert digest(ROOT / name) == expected, name
for dependency in PREVIOUS["externs"].values():
    assert digest(Path(dependency["path"])) == dependency["sha256"]
arity = PREVIOUS["generated_arities"]
assert digest(Path(arity["path"])) == arity["sha256"]
REPORT.mkdir()
shutil.copytree(ROOT / "vendor/rustred/crates/rustred-core", ISOLATED)
patch = ""
for relative, old, new in [
    (
        "src/solver/index.rs",
        "            .then(delta_order)\n            .then(absolute_right.cmp(&absolute_left))",
        "            .then(delta_order)\n"
        "            // DIAGNOSTIC ONLY: preserve sectors and required-cut lowering.\n"
        "            .then(occupation_right.cmp(&occupation_left))\n"
        "            .then(absolute_right.cmp(&absolute_left))",
    ),
    (
        "src/solver/guarded/persistence.rs",
        "rustred.guarded-source-program.v2",
        "rustred.guarded-source-program.experiment-cut-then-occupation-degree.v1",
    ),
]:
    path = ISOLATED / relative
    before = path.read_text()
    assert old in before
    after = before.replace(old, new)
    path.write_text(after)
    patch += "".join(difflib.unified_diff(
        before.splitlines(True), after.splitlines(True),
        fromfile="a/" + relative, tofile="b/" + relative,
    ))
(REPORT / "experimental.patch").write_text(patch)
command = [part.replace(
    "/tmp/rustflow-occupation-order-core-20261010", str(ISOLATED)
).replace(
    "/tmp/librustred-occupation-degree-first-20261010.rlib", str(OUTPUT)
).replace(
    "metadata=rustflow_occupation_degree_first_20261010",
    "metadata=rustflow_cut_then_occupation_degree_20261010",
) for part in PREVIOUS["command"]]
binding = {
    "scope": "isolated cut-priority-then-occupation-degree comparator; fresh native proofs only",
    "base_binding_sha256": digest(BASE / "build-binding.json"),
    "launcher_sha256": digest(Path(__file__)),
    "patch_sha256": digest(REPORT / "experimental.patch"),
    "isolated_sources": {
        str(p.relative_to(ISOLATED)): digest(p)
        for p in ISOLATED.rglob("*") if p.is_file()
    },
    "command": command,
    "output": str(OUTPUT),
    "resource_scope_correction": "Includes complete native compilation; excludes Nix startup, Python preparation and hashing. Shared resource wrapper generic compilation-excluded text does not apply.",
}
path = REPORT / "build-binding.json"
path.write_text(json.dumps(binding, indent=2) + "\n")
result = subprocess.run([
    "python3", str(ROOT / "reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py"),
    str(REPORT / "build-resources.json"), *command,
], cwd=ROOT)
binding["exit_code"] = result.returncode
binding["isolated_sources_unchanged"] = all(
    digest(ISOLATED / name) == expected
    for name, expected in binding["isolated_sources"].items()
)
assert binding["isolated_sources_unchanged"]
if OUTPUT.exists():
    binding["output_sha256"] = digest(OUTPUT)
path.write_text(json.dumps(binding, indent=2) + "\n")
raise SystemExit(result.returncode)
