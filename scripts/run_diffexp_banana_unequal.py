#!/usr/bin/env python3
"""Run unchanged original DiffExp with one CPU and polled time/RSS bounds.

This Linux runner requires an already licensed kernel. It installs nothing and
does not modify the upstream checkout. Memory is sampled every two seconds;
the observed peak can exceed the configured cap between samples.
"""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def process_tree_rss(pid):
    total = 0
    pending = [pid]
    seen = set()
    while pending:
        child = pending.pop()
        if child in seen:
            continue
        seen.add(child)
        try:
            status = Path(f"/proc/{child}/status").read_text()
            total += next(
                (int(line.split()[1]) * 1024 for line in status.splitlines()
                 if line.startswith("VmRSS:")), 0,
            )
            for task in Path(f"/proc/{child}/task").iterdir():
                try:
                    pending.extend(map(int, (task / "children").read_text().split()))
                except (FileNotFoundError, ProcessLookupError, PermissionError):
                    pass
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return total


def stop_group(process):
    def send(sig):
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            pass

    send(signal.SIGTERM)
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    # Clean remaining children even if the leader exited immediately.
    send(signal.SIGKILL)
    process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--route", choices=("full", "shortcut"), default="shortcut")
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--rss-bytes", type=int, default=16 * 1024**3)
    args = parser.parse_args()
    if args.timeout <= 0 or args.rss_bytes <= 0:
        parser.error("Resource limits must be positive")
    if args.cpu not in os.sched_getaffinity(0):
        parser.error("CPU is outside this process's available affinity")
    kernel = args.kernel.absolute()  # Preserve launcher symlinks.
    if not kernel.is_file() or not os.access(kernel, os.X_OK):
        parser.error("--kernel must name an executable file")
    upstream = args.upstream.absolute()
    package = upstream / "DiffExp.m"
    expected = "67fa23a0be32f747e292debcf3142b0ecef2cf526335b5b31f5ec7dd9cdeace7"
    if not package.is_file() or sha256(package) != expected:
        parser.error("Pinned DiffExp package mismatch")
    output = args.output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    driver = Path(__file__).with_name("diffexp_banana_unequal_oracle.wl").absolute()
    (output / "driver.wl").write_bytes(driver.read_bytes())
    (output / "runner.py").write_bytes(Path(__file__).read_bytes())
    command = [str(kernel), "-noinit", "-noprompt", "-script", str(driver)]
    env = dict(os.environ)
    env.update(
        DIFFEXP_ROOT=str(upstream), DIFFEXP_OUTPUT=str(output / "work"),
        DIFFEXP_BANANA_ROUTE=args.route, OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
    )
    os.sched_setaffinity(0, {args.cpu})
    start = time.monotonic()
    peak_rss = 0
    reason = None
    metadata = {
        "schema_version": 1, "command": command, "route": args.route,
        "upstream_commit": "784c8229bf92369a03f011a48e161522c8c54bbd",
        "package_sha256": expected, "driver_sha256": sha256(driver),
        "runner_sha256": sha256(Path(__file__)), "cpu_affinity": [args.cpu],
        "start_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "limits": {"seconds": args.timeout, "rss_bytes": args.rss_bytes},
        "poll_period_seconds": 2,
        "environment": {key: env[key] for key in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        )},
    }
    with (output / "kernel.log").open("w") as log, \
            (output / "monitor.jsonl").open("w") as events:
        process = subprocess.Popen(
            command, cwd=output, stdout=log, stderr=subprocess.STDOUT,
            env=env, start_new_session=True,
        )
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - start
                rss = process_tree_rss(process.pid)
                peak_rss = max(peak_rss, rss)
                events.write(json.dumps({
                    "elapsed_seconds": elapsed, "rss_bytes": rss,
                    "pid": process.pid,
                    "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                }) + "\n")
                events.flush()
                if elapsed >= args.timeout:
                    reason = "wall-clock limit"
                if rss > args.rss_bytes:
                    reason = "process-tree RSS limit"
                if reason:
                    stop_group(process)
                    break
                time.sleep(2)
        except BaseException:
            stop_group(process)
            raise
        code = process.wait()
    metadata["kernel_exit_code"] = code
    if reason is not None:
        code = 124
    if code == 0:
        try:
            result_file = output / "work/result.json"
            result = json.loads(result_file.read_text())
            if result.get("status") != "passed":
                raise ValueError("Missing successful numerical assertion")
            metadata["result_sha256"] = sha256(result_file)
        except (OSError, ValueError) as error:
            code = 2
            reason = str(error)
    metadata.update(
        exit_code=code, wall_seconds=time.monotonic() - start,
        peak_polled_process_tree_rss_bytes=peak_rss, termination_reason=reason,
        end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    )
    (output / "runner-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    return code if code >= 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
