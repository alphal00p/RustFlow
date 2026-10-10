#!/usr/bin/env python3
"""Run one bounded source-portfolio pilot after the small controls have passed."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("corpus", type=Path)
parser.add_argument("tag")
parser.add_argument("--engine", choices=["baseline", "proposal"], default="proposal")
args = parser.parse_args()
assert args.tag and all(c.isalnum() or c == "-" for c in args.tag)
base = Path(__file__).resolve().parent
root = base.parents[2]
out = base / "active-pilot" / args.tag
assert not out.exists()
build_path = (base / "active-pilot/build-binding.json" if args.engine == "proposal" else
              root / "reports/validation/2026-10-10-global-boost-source-experiment/native-active-pilot/pilot-build.json")
build = json.loads(build_path.read_text())
assert build["exit_code"] == 0
executable = Path(build["executable"]["path"])
corpus = args.corpus.resolve()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert digest(executable) == build["executable"]["sha256"]
native_binding_path = (build_path if args.engine == "proposal" else
                       build_path.with_name("pilot-build-binding.json"))
native_binding = json.loads(native_binding_path.read_text())
assert digest(Path(native_binding["source"]["path"])) == native_binding["source"]["sha256"]
for dependency in native_binding["dependencies"].values():
    assert digest(Path(dependency["path"])) == dependency["sha256"]
settings = {
    "PILOT_REQUESTS": "4096", "PILOT_ROUNDS": "16", "PILOT_FRONTIER": "1024",
    "PILOT_RULES": "65536", "PILOT_ZERO_ATTEMPTS": "8192", "PILOT_DEPTH": "3",
    "PILOT_RAY_DOMAINS": "1", "PILOT_POINT_DOMAINS": "1", "PILOT_REPLAY_OLD": "false",
}
out.mkdir()
command = ["timeout", "300", str(executable), str(corpus), str(out)]
record = {
    "scope": "Bounded native source-only active closure diagnostic; complete original target sums and derivatives. No predictions or numerical acceptance.",
    "command": command, "settings": settings,
    "engine": args.engine,
    "executable_sha256": digest(executable), "source_corpus_sha256": digest(corpus),
    "build_binding_sha256": digest(build_path), "launcher_sha256": digest(Path(__file__)),
    "native_binding_sha256": digest(native_binding_path),
}
binding = out / "run-binding.json"
binding.write_text(json.dumps(record, indent=2) + "\n")
result = subprocess.run([
    sys.executable, str(root / "reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py"),
    str(out / "resources.json"), *command,
], cwd=root, env={**os.environ, **settings})
assert digest(corpus) == record["source_corpus_sha256"]
assert digest(executable) == record["executable_sha256"]
record["exit_code"] = result.returncode
record["post_run_inputs_unchanged"] = True
binding.write_text(json.dumps(record, indent=2) + "\n")
raise SystemExit(result.returncode)
