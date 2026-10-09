"""Run one validation command and record process runtime/RSS, excluding builds.

Usage: python3 run-resource-command.py OUTPUT_JSON COMMAND [ARGS...]
Linux resource.ru_maxrss is reported in KiB. Child output is also saved as .log.
"""
import json
import pathlib
import resource
import subprocess
import sys
import time

output = pathlib.Path(sys.argv[1])
command = sys.argv[2:]
output.parent.mkdir(parents=True, exist_ok=True)
started = time.monotonic()
with output.with_suffix(".log").open("w") as log:
    completed = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True)
    for line in completed.stdout:
        log.write(line)
        sys.stdout.write(line)
        sys.stdout.flush()
    completed.wait()
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
output.write_text(json.dumps({
    "command": command,
    "exit_code": completed.returncode,
    "wall_seconds": time.monotonic() - started,
    "user_seconds": usage.ru_utime,
    "system_seconds": usage.ru_stime,
    "peak_child_rss_kib": usage.ru_maxrss,
    "scope": "one prebuilt validation process; compilation and Nix startup excluded",
}, indent=2) + "\n")
sys.exit(completed.returncode)
