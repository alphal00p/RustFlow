#!/usr/bin/env python3
"""Verify lossless payload hashes and replay the saved old/new comparison."""
import gzip
import hashlib
import json
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "archive-manifest.json").read_text())
    payloads = {}
    for entry in manifest["archives"]:
        packed = (root / entry["archive"]).read_bytes()
        assert len(packed) == entry["archive_bytes"]
        assert digest(packed) == entry["archive_sha256"]
        assert packed[4:8] == bytes(4), "gzip timestamp must remain zero"
        raw = gzip.decompress(packed)
        assert len(raw) == entry["original_bytes"]
        assert digest(raw) == entry["original_sha256"]
        original = root / entry["original"]
        if original.exists():
            assert original.read_bytes() == raw
        payloads[entry["original"]] = json.loads(raw)
    for name, entry in manifest["retained_bindings"].items():
        raw = (root / name).read_bytes()
        assert len(raw) == entry["bytes"] and digest(raw) == entry["sha256"], name
    comparison = json.loads((root / "zero-projection/comparison.json").read_text())
    for control in comparison["controls"]:
        mode = control["mode"]
        old = payloads[f"{mode}.json"]["probes"]
        new = payloads[f"zero-projection/{mode}.json"]["probes"]
        assert len(old) == len(new) == control["points"] == 19
        assert [p["target"] for p in old] == [p["target"] for p in new]
        assert sum(a["rhs"] == b["rhs"] for a, b in zip(old, new)) == control["rhs_identical_to_old"]
        assert sum(a["conditions"] == b["conditions"] for a, b in zip(old, new)) == control["conditions_identical_to_old"]
        assert sum(len(p["rhs"]) for p in old) == control["old_rhs_terms"]
        assert sum(len(p["rhs"]) for p in new) == control["new_rhs_terms"]
    print(json.dumps({"status": "pass", "archives_verified": len(payloads), "matched_old_new_points": 76}))


if __name__ == "__main__":
    main()
