import datetime, hashlib, json, os, pathlib, subprocess, time
root = pathlib.Path(__file__).resolve().parent
source = json.loads((root / "frozen-source-files.json").read_text())
def check_source():
    assert all(hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest() == v for k, v in source["files"].items()), "frozen source changed"
check_source()
artifacts = {}
for line in (root / "final-build.jsonl").read_text().splitlines():
    item = json.loads(line)
    if item.get("reason") == "compiler-artifact" and "amflow-dependent-cut-inventory" in item["manifest_path"] and item.get("executable") and item["profile"]["test"]:
        artifacts[item["target"]["name"]] = item["executable"]
checks = [("exact-unit", "symbolica_amflow", ["cut_inventory_tests"]), ("cut-private", "symbolica_amflow", ["cuts::tests"])]
checks += [(name, name, []) for name in ["dependent_cut_inventory", "native_cut_interfaces", "python_cut_interfaces", "partial_cut_boundaries", "massless_phase_space", "mixed_cut_flow", "cuts", "phase_space", "hepkit_cuts", "hepkit", "combinations"]]
results = {"started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "cpu_affinity": [44, 45], "test_threads": 2, "checks": []}
for label, name, filters in checks:
    check_source()
    binary = pathlib.Path(artifacts[name])
    args = ["taskset", "-c", "44,45", str(binary), *filters, "--nocapture", "--test-threads=2"]
    log = root / ("verified-" + label + ".log")
    started = time.monotonic()
    with log.open("wb") as output:
        run = subprocess.run(args, stdout=output, stderr=subprocess.STDOUT, timeout=1800)
    record = {"name": label, "command": args, "exit_code": run.returncode, "wall_seconds": time.monotonic() - started, "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(), "log": log.name, "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest()}
    results["checks"].append(record)
    (root / "verified-tests.json").write_text(json.dumps(results, indent=2) + "\n")
    print(label, run.returncode, round(record["wall_seconds"], 3), flush=True)
check_source()
results["completed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
results["passed"] = all(x["exit_code"] == 0 for x in results["checks"])
(root / "verified-tests.json").write_text(json.dumps(results, indent=2) + "\n")
raise SystemExit(0 if results["passed"] else 1)
