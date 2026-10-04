#!/usr/bin/env python3
"""Fetch the pinned scientific inputs for the original DiffExp benchmarks.

Uses Python's standard library. Archive and member hashes come from the tracked
provenance report. Existing files must match; unrelated files are left alone.
No Mathematica code is executed, and no upstream checkout is modified.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reports/diffexp/ancillary-provenance.json"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify(path, expected):
    if path.stat().st_size != expected["bytes"]:
        raise ValueError(f"Unexpected size: {path}")
    if digest(path) != expected["sha256"]:
        raise ValueError(f"SHA256 mismatch: {path}")


def atomic_bytes(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def fetch_archive(output, record):
    archive = output / (record["version"] + ".tar.gz")
    if archive.exists():
        verify(archive, record)
        return archive
    request = urllib.request.Request(
        record["url"],
        headers={"User-Agent": "RustFlow scientific benchmark reproducibility"},
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        # These archives are under 6 MB. Bound the read by the pinned size.
        data = response.read(record["bytes"] + 1)
    if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise ValueError(f"Downloaded archive differs from the pin: {record['url']}")
    atomic_bytes(archive, data)
    return archive


def extract_members(archive, output, record):
    with tarfile.open(archive, "r:gz") as stream:
        for expected in record["members"]:
            relative = PurePosixPath(expected["member"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe manifest path: {relative}")
            member = stream.getmember(str(relative))
            if not member.isfile() or member.size != expected["size"]:
                raise ValueError(f"Unexpected archive member: {relative}")
            destination = output / record["version"] / Path(*relative.parts)
            if not destination.resolve().is_relative_to(output.resolve()):
                raise ValueError(f"Extraction would leave the output directory: {relative}")
            if destination.exists():
                verify(destination, {"bytes": expected["size"], "sha256": expected["sha256"]})
                continue
            with stream.extractfile(member) as source:
                data = source.read()
            if hashlib.sha256(data).hexdigest() != expected["sha256"]:
                raise ValueError(f"Member SHA256 mismatch: {relative}")
            atomic_bytes(destination, data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "target/diffexp-oracle/ancillary"
    )
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    for record in manifest["archives"]:
        archive = fetch_archive(args.output, record)
        extract_members(archive, args.output, record)
        print(f"Verified {record['version']}: {len(record['members'])} ancillary files")
    print("Known source-data caveat: v2 mzz point 7 has an incomplete epsilon-zero row.")


if __name__ == "__main__":
    main()
