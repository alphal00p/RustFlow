#!/usr/bin/env python3
"""Run a bounded analytic oracle with the original AMFlow, Kira and Wolfram.

No runtime/reducer installation or license configuration is performed. The
kernel executable must already work, including its spawned child kernels.
The private runtime wrapper explicitly limits both Wolfram thread pools to one.
"""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


COMMIT = "26005517a288086c4cb4d1b26d829691bc088485"
SOURCE_HASHES = {
    "AMFlow.m": "76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184",
    "diffeq_solver/DESolver.m": "0dfa7f91bd8b8d160412dbc3a8c076413e488745bcefe58bc434999d03f76a46",
    "ibp_interface/Kira/interface.m": "1c4476b32d8757ee18d6689d5b406bc00a2cfb00ad178926233f1768def7faaa",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wl_string(value):
    """Encode a local Unix filename as a Wolfram string, without shell parsing."""
    if any(ord(c) < 32 for c in value):
        raise ValueError("Control characters in executable paths are unsupported")
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def signal_group(pid, sig):
    try:
        os.killpg(pid, sig)
    except ProcessLookupError:
        # The kernel may finish between the timeout and our signal.
        pass


def executable(value):
    # Keep wrappers/symlinks: some licensed launchers use their own location.
    path = Path(value).absolute()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise argparse.ArgumentTypeError(f"Not an executable file: {path}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--kernel", required=True, type=executable)
    parser.add_argument("--kira", required=True, type=executable)
    parser.add_argument("--fermat", required=True, type=executable)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", choices=("bubble", "sunset"), default="bubble")
    parser.add_argument("--mode", choices=("sample", "laurent"), default="sample")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    affinity = os.sched_getaffinity(0)
    if args.cpu is not None and args.cpu not in affinity:
        parser.error("--cpu must be available in this process's CPU affinity")
    for name, expected in SOURCE_HASHES.items():
        path = args.upstream / name
        if not path.is_file() or sha256(path) != expected:
            parser.error(f"Pinned upstream source mismatch: {name}")
    output = args.output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    upstream = output / "upstream"
    shutil.copytree(args.upstream, upstream)
    install = upstream / "ibp_interface/Kira/install.m"
    install.write_text(
        "(* Local install configuration; original computational files unchanged. *)\n"
        f"$KiraExecutable = {wl_string(str(args.kira))};\n"
        f"$FermatExecutable = {wl_string(str(args.fermat))};\n"
        '$RatracerExecutable = "";\n'
    )
    work = output / "work"
    work.mkdir()
    driver = output / "oracle.wl"
    shutil.copyfile(Path(__file__).with_suffix(".wl"), driver)
    shutil.copyfile(Path(__file__), output / "runner.py")
    thread_launcher = output / "wolfram_one_thread.py"
    thread_loader = output / "wolfram_one_thread.wl"
    for destination in [thread_launcher, thread_loader]:
        shutil.copyfile(Path(__file__).with_name(destination.name), destination)
    thread_launcher.chmod(0o755)
    env = dict(os.environ)
    env.update(
        AMFLOW_ROOT=str(upstream),
        AMFLOW_CONTROLLED_KERNEL=str(thread_launcher),
        RUSTFLOW_WOLFRAM_KERNEL=str(args.kernel),
        AMFLOW_ORACLE_WORKDIR=str(work),
        AMFLOW_ORACLE_CASE=args.case,
        AMFLOW_ORACLE_MODE=args.mode,
        FERMATPATH=str(args.fermat),
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    version = subprocess.run(
        [str(args.kira), "-v"], cwd=output, env=env, text=True,
        capture_output=True, timeout=15, check=True,
    ).stdout
    command = [str(thread_launcher), "-noinit", "-noprompt", "-script", str(driver)]
    metadata = {
        "upstream_commit": COMMIT,
        "upstream_source": str(args.upstream.absolute()),
        "source_sha256": SOURCE_HASHES,
        "modified_upstream_file": "ibp_interface/Kira/install.m",
        "install_sha256": sha256(install),
        "driver_sha256": sha256(driver),
        "runner_sha256": sha256(Path(__file__)),
        "kernel_command": command,
        "licensed_kernel_launcher": str(args.kernel),
        "runtime_thread_controls": {
            "ParallelThreadNumber": 1,
            "MKLThreadNumber": 1,
            "scope": "Parent and AMFlow-generated Wolfram child scripts",
            "launcher_sha256": sha256(thread_launcher),
            "loader_sha256": sha256(thread_loader),
        },
        "kira_executable": str(args.kira),
        "kira_version": version.strip(),
        "fermat_executable": str(args.fermat),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "workers": 1,
        "timeout_seconds": args.timeout,
        "start_utc": utc_now(),
    }
    print(output, flush=True)
    start = time.monotonic()
    with (output / "kernel.log").open("w") as log:
        process = subprocess.Popen(
            command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            code = process.wait(timeout=args.timeout)
            metadata["timed_out"] = False
        except subprocess.TimeoutExpired:
            metadata["timed_out"] = True
            signal_group(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                pass
            # Reap descendants even if the leader exited promptly after TERM.
            signal_group(process.pid, signal.SIGKILL)
            process.wait()
            code = 124
    metadata.update(
        exit_code=code, end_utc=utc_now(), wall_seconds=time.monotonic()-start,
    )
    result = work / "oracle-metadata.json"
    if code == 0:
        try:
            parsed = json.loads(result.read_text())
            if parsed.get("analytic_check_passed") is not True:
                raise ValueError("The oracle did not assert analytic agreement")
            metadata["result_sha256"] = sha256(result)
        except (OSError, ValueError) as error:
            metadata["validation_error"] = str(error)
            metadata["exit_code"] = code = 1
    (output / "runner-metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(f"Oracle exit code: {code}; log: {output / 'kernel.log'}")
    return code


if __name__ == "__main__":
    sys.exit(main())
