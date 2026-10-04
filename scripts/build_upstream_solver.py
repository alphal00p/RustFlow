#!/usr/bin/env python3
"""Build the unmodified, pinned AMFlow 2.0 standalone C++ solver using Nix deps."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import urllib.request

COMMIT = "26005517a288086c4cb4d1b26d829691bc088485"
ARCHIVE_SHA256 = "895fdc59111752a610a1ce9c13cd29c81db005221c4a9e86cc9ad33679fd3249"
ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--nixpkgs", help="default: revision from repository flake.lock")
    args = parser.parse_args()
    if args.jobs < 1 or args.output.exists():
        parser.error("jobs must be positive and output must be a new directory")
    output = args.output.resolve()
    output.mkdir(parents=True)
    url = f"https://gitlab.com/multiloop-pku/amflow/-/archive/{COMMIT}/amflow-{COMMIT}.tar.gz"
    archive = output / "upstream.tar.gz"
    urllib.request.urlretrieve(url, archive)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise RuntimeError("upstream archive digest differs; refusing unverified source")
    with tarfile.open(archive) as source:
        source.extractall(output, filter="data")
    source = output / f"amflow-{COMMIT}" / "diffeq_solver"
    nixpkgs = args.nixpkgs or ("github:NixOS/nixpkgs/" + json.loads(
        (ROOT / "flake.lock").read_text())["nodes"]["nixpkgs"]["locked"]["rev"])
    packages = ["gmp.dev", "gmp.out", "mpfr.dev", "mpfr.out", "libmpc.out",
                "boost.dev", "yaml-cpp.out", "mpsolve.out"]
    paths = []
    for package in packages:
        path = subprocess.check_output(["nix", "build", "--no-link", "--print-out-paths",
                                        f"{nixpkgs}#{package}"], text=True).strip()
        if len(path.splitlines()) != 1:
            raise RuntimeError(f"unexpected Nix output for {package}")
        paths.append(path)
    config = {
        "AMFLOW_GMP_INCLUDE": paths[0] + "/include", "AMFLOW_GMP_LIB": paths[1] + "/lib",
        "AMFLOW_MPFR_INCLUDE": paths[2] + "/include", "AMFLOW_MPFR_LIB": paths[3] + "/lib",
        "AMFLOW_MPC_INCLUDE": paths[4] + "/include", "AMFLOW_MPC_LIB": paths[4] + "/lib",
        "AMFLOW_BOOST_INCLUDE": paths[5] + "/include",
        "AMFLOW_YAML_CPP_INCLUDE": paths[6] + "/include",
        "AMFLOW_YAML_CPP_LIB": paths[6] + "/lib", "MPSOLVE": paths[7] + "/bin/mpsolve",
    }
    (source / "performance.config").write_text("".join(f"{k}={v}\n" for k, v in config.items()))
    command = ["make", "-j", str(args.jobs), "desolver", "CONFIG_FILE=performance.config",
               "OPENMP=no", "QUAD=no", "STATIC=no", "DEBUG=no"]
    with (output / "build.log").open("w") as log:
        subprocess.run(command, cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
    provenance = {"upstream_commit": COMMIT, "archive_url": url,
                  "archive_sha256": ARCHIVE_SHA256, "nixpkgs": nixpkgs,
                  "dependencies": dict(zip(packages, paths)), "command": command,
                  "compiler": subprocess.check_output(["g++", "--version"], text=True),
                  "first_compile_command": next(line for line in
                      (output / "build.log").read_text().splitlines() if " -c src/" in line),
                  "native_arch_flag_filtered_by_nix": "Skipping impure flag -march=native" in
                      (output / "build.log").read_text(),
                  "binary_sha256": hashlib.sha256((source / "desolver").read_bytes()).hexdigest()}
    (output / "build.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(source / "desolver")


if __name__ == "__main__":
    main()
