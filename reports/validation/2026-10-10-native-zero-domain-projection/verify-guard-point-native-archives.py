#!/usr/bin/env python3
"""Verify lossless guard-point archives and unchanged human-readable evidence."""
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
    manifest = json.loads((directory / "guard-point-native-archives.json").read_text())
    for entry in manifest["records"]:
        archive = root / entry["archive_path"]
        with archive.open("rb") as stream:
            assert fingerprint(stream) == (entry["archive_bytes"], entry["archive_sha256"])
        with archive.open("rb") as stream:
            header = stream.read(10)
            assert header[3] == 0 and header[4:8] == bytes(4), "nondeterministic gzip header"
        with gzip.open(archive, "rb") as stream:
            assert fingerprint(stream) == (entry["original_bytes"], entry["original_sha256"])
        original = root / entry["original_path"]
        if original.exists():
            with original.open("rb") as stream:
                assert fingerprint(stream) == (entry["original_bytes"], entry["original_sha256"])
    for name, entry in manifest["preserved_files"].items():
        with (root / name).open("rb") as stream:
            assert fingerprint(stream) == (entry["bytes"], entry["sha256"]), name
    print(json.dumps({
        "status": "pass",
        "archives": len(manifest["records"]),
        "preserved_files": len(manifest["preserved_files"]),
        "totals": manifest["totals"],
        "scope": "byte-for-byte archive verification only; physical acceptance unchanged",
    }))


if __name__ == "__main__":
    main()
