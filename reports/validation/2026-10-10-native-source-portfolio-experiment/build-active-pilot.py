#!/usr/bin/env python3
"""Build the isolated native active-frontier diagnostic, without Cargo."""
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
OUT = BASE / "active-pilot"
SOURCE = OUT / "pilot.rs"
EXE = Path("/tmp/rustred-source-portfolio-active-pilot-20261010")
LIBRARY = Path("/tmp/librustred-source-portfolio-20261010.rlib")
BUILD = OUT / "build-binding.json"


def meta(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


assert not EXE.exists() and not BUILD.exists()
native = json.loads((BASE / "native-build/build-binding.json").read_text())
assert native["exit_code"] == 0 and native["output_sha256"] == meta(LIBRARY)["sha256"]
dependencies = {
    "rustred": LIBRARY,
    "symbolica": ROOT / "target/release/deps/libsymbolica-02ed301b2d4bd7f2.rlib",
    "serde_json": ROOT / "target/release/deps/libserde_json-84255895f7dd0f90.rlib",
    "bincode": ROOT / "target/release/deps/libbincode-30fae4fd0bc9d1b5.rlib",
}
command = ["rustc", "--edition=2024", "--crate-name", "portfolio_active_pilot",
           "-C", "opt-level=0", "-C", "debuginfo=0"]
for name, path in dependencies.items():
    command += ["--extern", name + "=" + str(path)]
command += ["-L", "dependency=" + str(ROOT / "target/release/deps")]
for directory in ["gmp-mpfr-sys-2e6d668657be102d", "gmp-mpfr-sys-c789de89ed67d3db"]:
    command += ["-L", "native=" + str(ROOT / "target/release/build" / directory / "out/lib")]
command += [str(SOURCE), "-o", str(EXE)]
record = {
    "source": meta(SOURCE), "dependencies": {k: meta(v) for k, v in dependencies.items()},
    "command": command, "cargo_invoked": False,
    "scope": "Standalone diagnostic wrapper at opt-level 0 with release native dependencies; no production or performance acceptance claim.",
}
BUILD.write_text(json.dumps(record, indent=2) + "\n")
start = time.monotonic()
with (OUT / "build.log").open("w") as log:
    result = subprocess.run(command, cwd=ROOT,
                            env={**os.environ, "CARGO_CRATE_NAME": "portfolio_active_pilot"},
                            stdout=log, stderr=subprocess.STDOUT)
record["exit_code"] = result.returncode
record["wall_seconds"] = time.monotonic() - start
assert meta(SOURCE) == record["source"]
assert all(meta(v) == record["dependencies"][k] for k, v in dependencies.items())
if result.returncode == 0:
    record["executable"] = meta(EXE)
BUILD.write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps({"exit_code": result.returncode, "wall_seconds": record["wall_seconds"]}))
raise SystemExit(result.returncode)
