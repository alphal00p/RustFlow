#!/usr/bin/env python3
"""Compare singular Frobenius matching against pinned original AMFlow C++.

Uses exact shared inputs and an analytic resonant solution. No Mathematica or
third-party Python packages are required. Integral values stay decimal strings.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import platform
import statistics

from benchmark_de import (
    COMMIT, ROOT, complex_components, digest, distance, inputs, metadata,
    run, rust_values, wl_list, wl_value,
)


CASE = {
    "name": "analytic_integer_resonance",
    "matrix": [["1/(1-eta)", "0"], ["1", "1/(eta*(1-eta))"]],
    "epsilon": "101/9999900",
    "start": ["0.5", "0"], "end": ["0.1", "0"],
    "boundary": [["1", "0"], ["0", "0"]],
}


def analytic_values():
    with localcontext() as context:
        context.prec = 150
        eta = Decimal(CASE["end"][0])
        first = 1 / (2 * (1 - eta))
        second = eta * (2 * eta).ln() / (2 * (1 - eta))
        return [[str(first), "0"], [str(second), "0"]]


def specifications(folder, digits, precision, order):
    specification = inputs(folder, CASE, digits, precision, order)
    (folder / "job.yaml").write_text(
        "Job: SingularMatching\nDiffEq: diffeq\nBoundaryValue: boundary\n"
        "Solution: solution\nExpansionPoint: 0\n"
        f"BoundaryPoint: {wl_value(CASE['start'])}\n"
        f"EstimationPoint: {wl_value(CASE['end'])}\n"
        f"WorkingPre: {precision}\nXOrder: {order}\nRationalizePre: 100\n"
        "ExtraXOrder: 20\nRunLength: 1000\nNThread: 1\nVariable: eta\n"
        "Direction: NegIm\nPrescription: 1\n"
        f"Numeric:\n  eps: {CASE['epsilon']}\n")
    return specification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rust", required=True, type=Path)
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--upstream-build-info", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rust-source-commit", required=True)
    parser.add_argument("--rust-library", type=Path,
                        help="immutable release rlib used for standalone linking")
    parser.add_argument("--rust-build-description",
                        help="compiler flags and link method used for the supplied driver")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--digits", type=int, default=20)
    args = parser.parse_args()
    if args.repeats < 1 or not 1 <= args.digits <= 60:
        parser.error("repeats must be positive; digits must be between 1 and 60")
    if args.cpu is not None and args.cpu not in os.sched_getaffinity(0):
        parser.error("requested CPU is outside the allowed affinity")
    if args.output.exists():
        parser.error("output must be a new directory")
    rust, upstream = str(args.rust.resolve()), str(args.upstream.resolve())
    upstream_build = json.loads(args.upstream_build_info.read_text())
    if (upstream_build["upstream_commit"] != COMMIT or
            upstream_build["binary_sha256"] != digest(Path(upstream))):
        parser.error("upstream executable does not match the pinned build record")
    result = {
        "schema": 1,
        "scope": "singular Frobenius construction and matching with supplied boundary",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "upstream_commit": COMMIT,
        "rust_source_commit": args.rust_source_commit,
        "rust_build_description": args.rust_build_description,
        "rust_executable_sha256": digest(Path(rust)),
        "upstream_executable_sha256": digest(Path(upstream)),
        "upstream_build": upstream_build,
        "benchmark_source_sha256": digest(ROOT / "examples/benchmark_singular.rs"),
        "harness_source_sha256": digest(Path(__file__)),
        "shared_harness_sha256": digest(ROOT / "scripts/benchmark_de.py"),
        "requested_digits": args.digits,
        "platform": platform.platform(),
        "cpu": args.cpu,
        "affinity": sorted(os.sched_getaffinity(0)),
        "cpu_model": next((v.split(":", 1)[1].strip() for v in
                           Path("/proc/cpuinfo").read_text().splitlines()
                           if v.startswith("model name")), "unknown"),
        "host_load_average_at_start": os.getloadavg(),
        "rustc": metadata(["rustc", "--version"]),
        "timing": "fresh processes; one excluded warmup; alternating order; no persistent cache",
        "threads": {"NThread": 1, "OMP_NUM_THREADS": 1, "RAYON_NUM_THREADS": 1},
        "analytic_solution": ["1/(2*(1-eta))", "eta*log(2*eta)/(2*(1-eta))"],
        "analytic_endpoint": [["0.5", "0"], ["0", "0"]],
        "comparison_limit": "Rust physical_limit is checked analytically; upstream is compared at eta=1/10. Neither run constructs an automatic integral boundary.",
        "runs": [],
    }
    if args.rust_library:
        result["rust_library"] = {
            "path": str(args.rust_library.resolve()),
            "sha256": digest(args.rust_library),
        }
    args.output.mkdir(parents=True)
    vectors = {}
    base_order = max(80, 4 * args.digits)
    for phase, precision, order in [
            ("base", args.digits + 40, base_order),
            ("refined", args.digits + 60, base_order + 32)]:
        folder = (args.output / phase).resolve()
        result[phase + "_input"] = specifications(folder, args.digits, precision, order)
        commands = {"rust": [rust, str(folder / "input.json")],
                    "upstream": [upstream, str(folder / "job.yaml")]}
        for index in range(-1, args.repeats if phase == "base" else 1):
            engines = ["rust", "upstream"] if index % 2 else ["upstream", "rust"]
            for engine in engines:
                label = f"{engine}-{'warmup' if index < 0 else index}"
                timing = run(commands[engine], folder, label, args.cpu)
                if engine == "rust":
                    parsed = rust_values(Path(timing["log"]))
                    values = parsed.pop("values")
                    timing["endpoint"] = parsed.pop("endpoint")
                    timing["rust_diagnostics"] = parsed
                    if parsed["max_log_power"] < 1:
                        raise RuntimeError("the benchmark did not exercise logarithmic solutions")
                else:
                    values = [complex_components(v) for v in
                              wl_list((folder / "solution").read_text())]
                    if "singular_matching" not in timing["upstream_timers_ns"]:
                        raise RuntimeError("missing original singular-matching phase timer")
                vectors[phase, engine] = values
                if index >= 0:
                    result["runs"].append(dict(timing, values=values, engine=engine,
                                               phase=phase, repetition=index))
    analytic = analytic_values()
    checks = {}
    for phase in ["base", "refined"]:
        checks[phase + "_cross_implementation"] = distance(
            vectors[phase, "rust"], vectors[phase, "upstream"])
    for engine in ["rust", "upstream"]:
        checks[engine + "_precision_order_refinement"] = distance(
            vectors["base", engine], vectors["refined", engine])
        checks[engine + "_all_repetitions"] = max(
            distance(r["values"], vectors["refined", engine])
            for r in result["runs"] if r["engine"] == engine)
        checks[engine + "_all_analytic"] = max(
            distance(r["values"], analytic)
            for r in result["runs"] if r["engine"] == engine)
    checks["rust_endpoint_all_analytic"] = max(
        distance(r["endpoint"], result["analytic_endpoint"])
        for r in result["runs"] if r["engine"] == "rust")
    result["errors"] = {key: str(value) for key, value in checks.items()}
    result["validated"] = all(v <= Decimal(10) ** -args.digits for v in checks.values())
    result["median_wall_ns"] = {engine: statistics.median(
        r["wall_ns"] for r in result["runs"]
        if r["engine"] == engine and r["phase"] == "base")
        for engine in ["rust", "upstream"]}
    result["median_singular_matching_ns"] = {
        engine: statistics.median(
            r["rust_diagnostics"]["singular_matching_ns"] if engine == "rust"
            else r["upstream_timers_ns"]["singular_matching"]
            for r in result["runs"] if r["engine"] == engine and r["phase"] == "base")
        for engine in ["rust", "upstream"]}
    (args.output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print("validated:", result["validated"], "median wall ns:", result["median_wall_ns"],
          "median singular matching ns:", result["median_singular_matching_ns"])
    if not result["validated"]:
        raise RuntimeError("accuracy check failed; timing comparison is not accepted")


if __name__ == "__main__":
    main()
