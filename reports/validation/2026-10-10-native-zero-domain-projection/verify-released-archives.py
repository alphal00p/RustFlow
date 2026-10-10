#!/usr/bin/env python3
"""Verify released gzip payloads and the unchanged pre-refinement manifest."""
import gzip
import hashlib
import json
from pathlib import Path


def fingerprint(stream):
    digest = hashlib.sha256()
    size = 0
    while block := stream.read(1024 * 1024):
        digest.update(block)
        size += len(block)
    return size, digest.hexdigest()


def main():
    directory = Path(__file__).resolve().parent
    root = directory.parents[2]
    manifest = json.loads((directory / "released-native-archives.json").read_text())
    frozen_path = root / manifest["frozen_manifest"]
    assert hashlib.sha256(frozen_path.read_bytes()).hexdigest() == manifest["frozen_manifest_sha256"]
    frozen = json.loads(frozen_path.read_text())
    mappings = {}
    for archive_manifest in [directory / "zero-projection-three-loop-archives.json", directory / "released-native-archives.json"]:
        for entry in json.loads(archive_manifest.read_text())["records"]:
            mappings[entry["original_path"]] = entry
    for entry in manifest["records"]:
        packed = root / entry["archive_path"]
        with packed.open("rb") as stream:
            assert fingerprint(stream) == (entry["archive_bytes"], entry["archive_sha256"])
        with packed.open("rb") as stream:
            assert stream.read(8)[4:] == bytes(4), "nonzero gzip timestamp"
        with gzip.open(packed, "rb") as stream:
            assert fingerprint(stream) == (entry["original_bytes"], entry["original_sha256"])
    verified = 0
    restored = 0
    for name, expected in {**frozen["files"], **frozen["external_artifacts"]}.items():
        path = root / name
        if path.is_file():
            stream = path.open("rb")
        else:
            entry = mappings[name]
            stream = gzip.open(root / entry["archive_path"], "rb")
            restored += 1
        with stream:
            assert fingerprint(stream) == (expected["size_bytes"], expected["sha256"]), name
        verified += 1
    print(json.dumps({"status": "pass", "new_archives": len(manifest["records"]), "frozen_entries_verified": verified, "frozen_entries_read_through_archives": restored, "frozen_manifest_unchanged": True}))


if __name__ == "__main__":
    main()
