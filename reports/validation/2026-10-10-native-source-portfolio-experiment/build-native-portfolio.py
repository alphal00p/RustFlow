#!/usr/bin/env python3
"""Build the isolated guarded source portfolio with its private diagnostic entry."""
import difflib
import hashlib
import json
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PRIOR = ROOT / "reports/validation/2026-10-10-occupation-order-experiment/build-binding.json"
PREVIOUS = json.loads(PRIOR.read_text())
ISOLATED = Path("/tmp/rustred-source-portfolio-20261010/crates/rustred-core")
REPORT = BASE / "native-build"
OUTPUT = Path("/tmp/librustred-source-portfolio-20261010.rlib")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert ISOLATED.is_dir() and not REPORT.exists() and not OUTPUT.exists()
for name, expected in PREVIOUS["base_native_sources"].items():
    assert digest(ROOT / name) == expected, name
for dependency in PREVIOUS["externs"].values():
    assert digest(Path(dependency["path"])) == dependency["sha256"]
arity = PREVIOUS["generated_arities"]
assert digest(Path(arity["path"])) == arity["sha256"]
REPORT.mkdir()
production = ROOT / "vendor/rustred/crates/rustred-core"
patch = ""
isolated_sources = {}
for path in sorted(ISOLATED.rglob("*")):
    if not path.is_file():
        continue
    relative = path.relative_to(ISOLATED)
    isolated_sources[str(relative)] = digest(path)
    original = production / relative
    if original.exists() and digest(original) == digest(path):
        continue
    before = original.read_text() if original.exists() else ""
    patch += "".join(difflib.unified_diff(
        before.splitlines(True), path.read_text().splitlines(True),
        fromfile="a/" + str(relative) if original.exists() else "/dev/null",
        tofile="b/" + str(relative),
    ))
(REPORT / "experimental.patch").write_text(patch)
command = [part.replace(
    "/tmp/rustflow-occupation-order-core-20261010", str(ISOLATED)
).replace(
    "/tmp/librustred-occupation-degree-first-20261010.rlib", str(OUTPUT)
).replace(
    "metadata=rustflow_occupation_degree_first_20261010",
    "metadata=rustflow_guarded_source_portfolio_20261010",
) for part in PREVIOUS["command"]]
command.insert(command.index("rustc"), "CARGO_CRATE_NAME=rustred")
command[command.index("--crate-type"):command.index("--crate-type")] = [
    "--cfg", "rustflow_portfolio_controls",
]
binding = {
    "scope": "isolated point-only guarded source portfolio with cfg-gated private controls; no production changes",
    "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "base_binding_sha256": digest(PRIOR),
    "launcher_sha256": digest(Path(__file__)),
    "patch_sha256": digest(REPORT / "experimental.patch"),
    "isolated_sources": isolated_sources,
    "externs": PREVIOUS["externs"],
    "generated_arities": arity,
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
    digest(ISOLATED / name) == expected for name, expected in isolated_sources.items()
)
binding["production_sources_unchanged"] = all(
    digest(ROOT / name) == expected
    for name, expected in PREVIOUS["base_native_sources"].items()
)
assert binding["isolated_sources_unchanged"] and binding["production_sources_unchanged"]
if OUTPUT.exists():
    binding["output_sha256"] = digest(OUTPUT)
path.write_text(json.dumps(binding, indent=2) + "\n")
raise SystemExit(result.returncode)
