#!/usr/bin/env python3
"""Bind the completed standalone validation artifacts without changing them."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    comparison = json.loads((ROOT / "comparison.json").read_text())
    assert comparison["status"] == "passed"
    assert len(comparison["comparisons"]) == 3
    for name in ("native-resources.json", "reference-resources.json"):
        assert json.loads((ROOT / name).read_text())["exit_code"] == 0
    target = ROOT / "artifact-manifest.json"
    if target.exists():
        raise SystemExit("Refusing to overwrite the completed manifest")
    files = {
        str(path.relative_to(ROOT)): {
            "sha256": digest(path),
            "bytes": path.stat().st_size,
        }
        for path in sorted(ROOT.rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts and path != target
    }
    target.write_text(json.dumps({
        "status": "passed",
        "scope": "Standalone ordinary connected hard-boundary coefficient only",
        "claimed_reference_digits": 10,
        "full_three_loop_finite_density_acceptance": False,
        "production_sources_changed": False,
        "cargo_invoked": False,
        "files": files,
    }, indent=2) + "\n")
    print(json.dumps({"files": len(files), "manifest_sha256": digest(target)}))


if __name__ == "__main__":
    main()
