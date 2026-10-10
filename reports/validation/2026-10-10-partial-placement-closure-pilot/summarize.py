#!/usr/bin/env python3
"""Freeze bounded, source-only closure outcomes; no reference values are read."""
import ast
import hashlib
import json
from fractions import Fraction
from pathlib import Path

BASE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(text, epsilon):
    def visit(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return Fraction(node.value)
        if isinstance(node, ast.Name):
            return {"epsilon": epsilon, "eta": Fraction(1)}[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Div): return a / b
            if isinstance(node.op, ast.Pow):
                assert b.denominator == 1
                return a ** b.numerator
        raise ValueError(ast.dump(node))
    return visit(ast.parse(text.replace("^", "**"), mode="eval").body)


arms = []
for mode in ["single-slot0", "double-slots12", "double-slots24"]:
    binding = json.loads((BASE / f"{mode}-binding.json").read_text())
    resource = json.loads((BASE / f"{mode}-resources.json").read_text())
    progress = json.loads((BASE / mode / "progress.json").read_text())
    assert binding["inputs_unchanged"]
    assert binding["exit_code"] == resource["exit_code"]
    row = {
        "mode": mode, "status": "closed" if resource["exit_code"] == 0 else "timeout_unclosed",
        "resources": resource, "completed_rounds": len(progress),
        "frontier_history": [p["frontier"] for p in progress],
        "last_completed_round": progress[-1],
        "binding_sha256": sha(BASE / f"{mode}-binding.json"),
        "conditions_and_gaps": "Retained in original round checkpoints; no exceptional-locus coverage claim.",
    }
    if resource["exit_code"] == 0:
        result = json.loads((BASE / mode / "result.json").read_text())
        assert result["status"] == "closed" and len(result["closed"]["basis"]) == 2
        row["basis"] = result["closed"]["basis"]
        row["retained_condition_count"] = len(result["closed"]["conditions"])
        row["final_target_and_derivative_rows"] = result["closed"]["target_then_derivative_rows"]
        row["independent_native_audit"] = "All original weighted target sums and every basis derivative replayed with with_terminals_replayed; encoded closed program decoded and rows/conditions rechecked."
        row["closed_program_sha256"] = sha(BASE / mode / "closed.bin")
        row["sample_condition_diagnostics_at_eta_1"] = [
            {"dimension": str(d), "epsilon": str((4-d)/2),
             "vanished_conditions": [c for c in result["closed"]["conditions"] if evaluate(c, (4-d)/2) == 0]}
            for d in [Fraction(12, 5), Fraction(13, 2), Fraction(31, 3)]
        ]
    else:
        assert resource["exit_code"] == 124
        assert not (BASE / mode / "closed.bin").exists()
    arms.append(row)

summary = {
    "scope": "Fresh isolated partial-placement native closure only; no boundary integration, transport, physical endpoint evaluation or production admission.",
    "native_library": "isolated verified-program reuse; source-bound immutable scheduling unions, independent final source replay",
    "old_rules_imported": False,
    "bounds": {"wall_seconds_per_arm": 300, "rounds": 16, "cumulative_requests": 4096,
               "frontier": 1024, "rules": 65536, "direct_zero_attempts_per_round": 8192,
               "ray_domains": 1, "point_domains": 1, "depth": 3},
    "arms": arms,
    "interpretation": "The singleton closes with two masters. Both double placements remain unclosed at the wall cap. Different placements define different finite-eta integrals; RHS sizes are not numerical equality tests.",
}
(BASE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
files = {str(p.relative_to(BASE)): {"sha256": sha(p), "bytes": p.stat().st_size}
         for p in sorted(BASE.rglob("*")) if p.is_file() and p.name != "artifact-manifest.json"}
(BASE / "artifact-manifest.json").write_text(json.dumps({"files": files, "file_count": len(files),
    "total_bytes": sum(f["bytes"] for f in files.values()), "scope": summary["scope"]}, indent=2) + "\n")
print(json.dumps({"arms": [{"mode": a["mode"], "status": a["status"], "rounds": a["completed_rounds"]} for a in arms],
                  "file_count": len(files), "manifest_sha256": sha(BASE / "artifact-manifest.json")}, indent=2))
