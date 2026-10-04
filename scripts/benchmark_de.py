#!/usr/bin/env python3
"""Compare Rust and pinned AMFlow C++ regular DE transport on identical input.

No Mathematica, Python numerical dependencies, evaluation of WL code, or binary
floating-point conversion of integral values is used. Timing uses OS counters.
Build both executables separately; see docs/performance.md.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "26005517a288086c4cb4d1b26d829691bc088485"
FIXTURES = ROOT / "fixtures/performance"
HASHES = {
    "upstream12-diffeq.wl": "828dd1ae41291007fd3b37b664e695d67e0daebde3e9330f0e178f3d39cc40e3",
    "upstream12-boundary.wl": "fcbe018f60ca9cf1bf909dec36d3608993d22fe5c8bb97a86df432870c1f92c6",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wl_list(text):
    """Only parse braces and comma separators; leave every scalar unevaluated."""
    text = re.sub(r"\s|\\\n", "", text)
    if not text.startswith("{") or not text.endswith("}"):
        raise ValueError("expected a Mathematica list")
    body, level, begin, items = text[1:-1], 0, 0, []
    for pos, char in enumerate(body):
        if char in "{(":
            level += 1
        elif char in "})":
            level -= 1
        elif char == "," and level == 0:
            items.append(body[begin:pos])
            begin = pos + 1
        if level < 0:
            raise ValueError("unbalanced input")
    if level != 0:
        raise ValueError("unbalanced input")
    items.append(body[begin:])
    return [wl_list(item) if item.startswith("{") else item for item in items]


def complex_components(text):
    text = re.sub(r"`[0-9.]*", "", text).replace("*^", "e")
    text = re.sub(r"\s|\\", "", text).replace("*I", "I")
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    match = re.fullmatch(f"({number})(?:({number})I)?", text)
    if not match:
        raise ValueError(f"not a real or complex number: {text!r}")
    return [match[1], match[2] or "0"]


def workloads():
    for name, expected in HASHES.items():
        if digest(FIXTURES / name) != expected:
            raise ValueError(f"modified pinned fixture: {name}")
    yield {
        "name": "analytic_log_chain",
        "matrix": [["1/(3*eta)", "1/eta", "0"],
                   ["0", "1/(3*eta)", "1/eta"],
                   ["0", "0", "1/(3*eta)"]],
        "epsilon": "101/9999900", "start": ["1", "0"], "end": ["4", "0"],
        "boundary": [["1", "0"]] * 3,
    }
    yield {
        "name": "upstream_12_master",
        "matrix": wl_list((FIXTURES / "upstream12-diffeq.wl").read_text()),
        "epsilon": "101/9999900", "start": ["0.5", "0"], "end": ["0.1", "0.2"],
        "boundary": [complex_components(v) for v in wl_list(
            (FIXTURES / "upstream12-boundary.wl").read_text())],
    }


def custom_case(path):
    """Read data only; never evaluate case strings as Python or Wolfram code."""
    case = json.loads(path.read_text())
    if not isinstance(case, dict) or not isinstance(case.get("name"), str):
        raise ValueError("case must be an object with a string name")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", case["name"]):
        raise ValueError("case name must contain only letters, digits, underscores or hyphens")
    matrix = case.get("matrix")
    if not isinstance(matrix, list) or not matrix:
        raise ValueError("case matrix must be a nonempty square array")
    size = len(matrix)
    if any(not isinstance(row, list) or len(row) != size or
           any(not isinstance(value, str) or not value.strip() for value in row)
           for row in matrix):
        raise ValueError("case matrix must be a square array of nonempty strings")
    boundary = case.get("boundary")
    if not isinstance(boundary, list) or len(boundary) != size:
        raise ValueError("case boundary dimension must match its matrix")
    for value in [case.get("start"), case.get("end"), *boundary]:
        if not isinstance(value, list) or len(value) != 2 or any(
                not isinstance(component, str) or not component.strip() for component in value):
            raise ValueError("complex values must be pairs of nonempty strings")
        for component in value:
            if not re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", component):
                raise ValueError("boundary and endpoint components must be finite decimal literals")
            # Decimal preserves the literal exactly; never convert values through float.
            if not Decimal(component).is_finite():
                raise ValueError("nonfinite boundary or endpoint component")
    epsilon = case.get("epsilon")
    if not isinstance(epsilon, str) or not re.fullmatch(r"[+-]?\d+(?:/[+-]?\d+)?", epsilon):
        raise ValueError("case epsilon must be an exact integer or integer fraction")
    if "/" in epsilon and int(epsilon.split("/", 1)[1]) == 0:
        raise ValueError("epsilon denominator cannot be zero")
    return case


def wl_value(value):
    real, imag = [s.replace("e", "*^").replace("E", "*^") for s in value]
    return real + ("+" if not imag.startswith(("-", "+")) else "") + imag + "*I"


def inputs(folder, case, digits, precision, order):
    folder.mkdir(parents=True, exist_ok=True)
    # Boost 1.89 multiprecision/detail/digits.hpp::digits10_2_2, used by upstream.
    bits = precision * 1000 // 301 + (2 if precision * 1000 % 301 else 1)
    data = dict(case, digits=digits, working_decimal_digits=precision,
                working_bits=bits, order=order)
    (folder / "input.json").write_text(json.dumps(data, indent=2) + "\n")
    (folder / "diffeq").write_text("{" + ",".join(
        "{" + ",".join(row) + "}" for row in case["matrix"]) + "}\n")
    (folder / "boundary").write_text("{" + ",".join(
        wl_value(v) for v in case["boundary"]) + "}\n")
    (folder / "job.yaml").write_text(
        "Job: RegularRun\nDiffEq: diffeq\nBoundaryValue: boundary\nSolution: solution\n"
        f"StartPoint: {wl_value(case['start'])}\nEndPoint: {wl_value(case['end'])}\n"
        f"WorkingPre: {precision}\nXOrder: {order}\nRationalizePre: 100\n"
        "ExtraXOrder: 20\nRunLength: 1000\nNThread: 1\nVariable: eta\n"
        f"Numeric:\n  eps: {case['epsilon']}\n")
    return data


def run(command, folder, label, cpu):
    env = dict(os.environ, OMP_NUM_THREADS="1", RAYON_NUM_THREADS="1")
    log = folder / f"{label}.log"
    affinity = None if cpu is None else lambda: os.sched_setaffinity(0, {cpu})
    with log.open("w") as output:
        start = time.monotonic_ns()
        process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT,
                                   env=env, preexec_fn=affinity)
        _, status, usage = os.wait4(process.pid, 0)
        elapsed = time.monotonic_ns() - start
        process.returncode = os.waitstatus_to_exitcode(status)
    if process.returncode:
        raise RuntimeError(f"failed command {command}; inspect {log}")
    record = {"wall_ns": elapsed, "user_seconds": usage.ru_utime,
            "system_seconds": usage.ru_stime, "max_rss_kib": usage.ru_maxrss,
            "command": command, "log": str(log)}
    if label.startswith("upstream"):
        record["upstream_timers_ns"] = {
            phase: int(Decimal(milliseconds) * 1_000_000)
            for milliseconds, phase in re.findall(
                r"Takes ([0-9.eE+-]+) ms on (\w+)\(\)", log.read_text())}
    return record


def rust_values(path):
    for line in reversed(path.read_text().splitlines()):
        if line.startswith('{"'):
            return json.loads(line)
    raise ValueError(f"no result in {path}")


def distance(left, right):
    """Normwise relative error, with absolute tolerance below unit magnitude."""
    if len(left) != len(right):
        raise ValueError("result vector dimensions differ")
    errors = []
    with localcontext() as context:
        context.prec = 160
        for a, b in zip(left, right):
            ar, ai, br, bi = map(Decimal, [*a, *b])
            if not all(v.is_finite() for v in [ar, ai, br, bi]):
                raise ValueError("nonfinite solution")
            scale = max(Decimal(1), (ar * ar + ai * ai).sqrt(), (br * br + bi * bi).sqrt())
            errors.append((((ar - br) ** 2 + (ai - bi) ** 2).sqrt() / scale))
        return max(errors)


def metadata(command):
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rust", required=True, type=Path)
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--upstream-build-info", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--rust-source-commit", required=True,
                        help="commit from which the Rust library was built")
    parser.add_argument("--digits", type=int, default=20)
    parser.add_argument("--case", type=Path, action="append",
                        help="use these supplied-boundary JSON cases instead of the two defaults")
    args = parser.parse_args()
    if args.repeats < 1 or args.digits < 1 or args.digits > 60:
        parser.error("repeats must be positive; digits must be between 1 and 60")
    cases = [custom_case(path) for path in args.case] if args.case else list(workloads())
    if len({case["name"] for case in cases}) != len(cases):
        parser.error("case names must be unique")
    if args.output.exists():
        parser.error("output must be a new directory, to preserve previous measurements")
    args.output.mkdir(parents=True)
    rust, upstream = str(args.rust.resolve()), str(args.upstream.resolve())
    upstream_build = json.loads(args.upstream_build_info.read_text())
    if (upstream_build["upstream_commit"] != COMMIT or
            upstream_build["binary_sha256"] != digest(Path(upstream))):
        parser.error("upstream executable does not match the pinned build provenance")
    result = {
        "schema": 1, "scope": "regular DE preparation and transport with supplied boundaries",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "upstream_commit": COMMIT, "rust_source_commit": args.rust_source_commit,
        "rust_executable_sha256": digest(Path(rust)),
        "upstream_executable_sha256": digest(Path(upstream)),
        "upstream_build": upstream_build,
        "benchmark_source_sha256": digest(ROOT / "examples/benchmark_de.rs"),
        "harness_source_sha256": digest(Path(__file__)),
        "fixture_sha256": ({str(path): digest(path) for path in args.case}
                           if args.case else HASHES), "requested_digits": args.digits,
        "platform": platform.platform(), "cpu": args.cpu,
        "affinity": sorted(os.sched_getaffinity(0)),
        "cpu_model": next((v.split(":", 1)[1].strip() for v in
                           Path("/proc/cpuinfo").read_text().splitlines()
                           if v.startswith("model name")), "unknown"),
        "host_load_average_at_start": os.getloadavg(),
        "rustc": metadata(["rustc", "--version"]),
        "cxx": metadata(["g++", "--version"]),
        "timing": "fresh process wall time and wait4 CPU/RSS; one excluded warmup; alternating order; no persistent cache",
        "threads": {"NThread": 1, "OMP_NUM_THREADS": 1, "RAYON_NUM_THREADS": 1},
        "cases": [],
    }
    for case in cases:
        report = {"name": case["name"], "dimension": len(case["matrix"]), "runs": []}
        vectors = {}
        for phase, precision, order in [("base", args.digits + 40, 80),
                                        ("refined", args.digits + 60, 112)]:
            folder = (args.output / case["name"] / phase).resolve()
            specification = inputs(folder, case, args.digits, precision, order)
            report[phase + "_input"] = specification
            commands = {"rust": [rust, str(folder / "input.json")],
                        "upstream": [upstream, str(folder / "job.yaml")]}
            for index in range(-1, args.repeats if phase == "base" else 1):
                engines = ["rust", "upstream"] if index % 2 else ["upstream", "rust"]
                for engine in engines:
                    label = f"{engine}-{'warmup' if index < 0 else index}"
                    timing = run(commands[engine], folder, label, args.cpu)
                    if engine == "rust":
                        parsed = rust_values(Path(timing["log"]))
                        value = parsed.pop("values")
                        timing["rust_diagnostics"] = parsed
                    else:
                        value = [complex_components(v) for v in
                                 wl_list((folder / "solution").read_text())]
                    vectors[phase, engine] = value
                    if index >= 0:
                        report["runs"].append(dict(timing, engine=engine, phase=phase,
                                                  repetition=index, values=value))
        checks = {
            "base_cross_implementation": distance(vectors["base", "rust"], vectors["base", "upstream"]),
            "refined_cross_implementation": distance(vectors["refined", "rust"], vectors["refined", "upstream"]),
            "rust_precision_order_refinement": distance(vectors["base", "rust"], vectors["refined", "rust"]),
            "upstream_precision_order_refinement": distance(vectors["base", "upstream"], vectors["refined", "upstream"]),
        }
        for engine in ["rust", "upstream"]:
            checks[engine + "_all_repetitions"] = max(
                distance(run["values"], vectors["refined", engine])
                for run in report["runs"] if run["engine"] == engine)
        if case["name"] == "analytic_log_chain":
            with localcontext() as context:
                context.prec = 150
                log = Decimal(4).ln()
                power = (log / 3).exp()
                analytic = [[str(v), "0"] for v in
                            [power * (1 + log + log * log / 2), power * (1 + log), power]]
            for engine in ["rust", "upstream"]:
                checks[engine + "_analytic"] = distance(vectors["base", engine], analytic)
        report["errors"] = {key: str(value) for key, value in checks.items()}
        report["validated"] = all(v <= Decimal(10) ** -args.digits for v in checks.values())
        report["median_wall_ns"] = {engine: statistics.median(
            r["wall_ns"] for r in report["runs"] if r["engine"] == engine and r["phase"] == "base")
            for engine in ["rust", "upstream"]}
        result["cases"].append(report)
        (args.output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
        print(case["name"], "validated:", report["validated"],
              "median wall ns:", report["median_wall_ns"], flush=True)
        if not report["validated"]:
            raise RuntimeError("accuracy check failed; timing is not a valid accuracy-matched comparison")


if __name__ == "__main__":
    main()
