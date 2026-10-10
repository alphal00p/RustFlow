"""Verify the frozen exploration evidence without rerunning its numerical jobs."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REPORTS = HERE.parent
MANIFESTS = {
    "2026-10-10-integrated-angular-review": "artifact-manifest.json",
    "2026-10-10-integrated-certificate-preflight": "artifact-manifest.json",
    "2026-10-10-integrated-continuation-independent-audit": "artifact-manifest.json",
    "2026-10-10-integrated-flow-conditional-transport": "artifact-manifest.json",
    "2026-10-10-integrated-flow-conditional-validation": "artifact-manifest.json",
    "2026-10-10-integrated-flow-physical-fit": "artifact-manifest.json",
    "2026-10-10-integrated-flow-reconstruction": "implementation-manifest.json",
    "2026-10-10-integrated-flow-target-protocol": "artifact-manifest.json",
    "2026-10-10-integrated-moment-continuation": "artifact-manifest.json",
    "2026-10-10-integrated-moment-expansion": "artifact-manifest.json",
    "2026-10-10-integrated-moment-series-comparison": "artifact-manifest.json",
    "2026-10-10-integrated-moment-validation": "artifact-manifest.json",
}


def binding(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).splitlines()


def main():
    failures = []
    results = []
    for report, filename in MANIFESTS.items():
        directory = REPORTS / report
        path = directory / filename
        if not path.is_file():
            failures.append(f"Missing manifest: {path}")
            continue
        manifest = json.loads(path.read_text())
        files = manifest["files"]
        entries = files.items() if isinstance(files, dict) else (
            (item["path"], item) for item in files
        )
        count = 0
        for relative, expected in entries:
            actual = binding(directory / relative)
            for key in ("bytes", "sha256"):
                if actual[key] != expected[key]:
                    failures.append(f"{report}/{relative}: {key} mismatch")
            count += 1
        results.append({"path": str(path.relative_to(ROOT)),
                        **binding(path), "verified_files": count})

    base = json.loads((HERE / "plan.json").read_text())["base_commit"]
    changed = git("diff", "--name-only", base)
    untracked = git("ls-files", "--others", "--exclude-standard")
    unexpected = sorted(path for path in set(changed + untracked)
                        if not path.startswith("reports/validation/")
                        and path != "docs/finite-density-pause-checkpoint.md")
    failures.extend(f"Unexpected production change: {path}" for path in unexpected)
    output = {
        "schema": 1,
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Artifact integrity and unchanged production tree; mathematical and numerical evidence is in the bound reports",
        "base_commit": base,
        "head_before_final_checkpoint": git("rev-parse", "HEAD")[0],
        "verifier": binding(Path(__file__).resolve()),
        "manifests": results,
        "verified_files": sum(result["verified_files"] for result in results),
        "production_unchanged": not unexpected,
        "failures": failures,
        "passed": not failures,
    }
    (HERE / "checkpoint-verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: output[key] for key in
                      ("verified_files", "production_unchanged", "failures", "passed")}))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
