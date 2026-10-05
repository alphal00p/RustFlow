#!/usr/bin/env python3
"""Prepare exact scalar HEFT inputs from an externally supplied pinned bridge.

No upstream implementation is bundled. The closed arithmetic expression reader
only evaluates the four external scalar assignments; tensor operations remain
owned by HEPKit/Idenso and the independent MG5/ALOHA oracle.
"""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

PINNED_BRIDGE_SHA256 = "e2df9f77089477730ca8b00ee591aa926bb6f242f1c7629880ab510357fad8f2"


def evaluate(node: ast.AST, variables: dict[str, Fraction]) -> Fraction:
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return Fraction(node.value)
    if isinstance(node, ast.Name):
        return variables[node.id]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -evaluate(node.operand, variables)
    if isinstance(node, ast.BinOp):
        left, right = evaluate(node.left, variables), evaluate(node.right, variables)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "pow"
        and len(node.args) == 2
        and not node.keywords
    ):
        power = evaluate(node.args[1], variables)
        if power.denominator == 1 and 0 <= power <= 8:
            return evaluate(node.args[0], variables) ** int(power)
    raise ValueError(f"unsupported exact scalar expression: {ast.dump(node)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--coherent-input", type=Path, required=True)
    parser.add_argument("--oracle-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("refusing to overwrite an existing HEFT input directory")
    raw = args.source.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    if source_hash != PINNED_BRIDGE_SHA256:
        raise ValueError("the external scalar bridge differs from the pinned validated source")
    coherent = json.loads(args.coherent_input.read_text())
    coordinates = coherent["physical_coordinates"]["physical_s_t_MH_squared"]
    if len(coordinates) != 3 or any(not isinstance(value, str) for value in coordinates):
        raise ValueError("phase coordinates must be exact strings")
    s, t, mass2 = map(Fraction, coordinates)
    assignments = re.findall(r"gghgHEFTTensor\[(\d+)\]\s*=([^;]+);", raw.decode())
    if [int(index) for index, _ in assignments] != [0, 1, 2, 3]:
        raise ValueError("external scalar assignments do not cover indices0 through3")
    values = {"s": s, "t": t, "mmH": mass2}
    coefficients = [
        {
            "index": int(index) + 1,
            "value": str(evaluate(ast.parse(expression.strip(), mode="eval").body, values)),
        }
        for index, expression in assignments
    ]
    source_report = json.loads(args.oracle_report.read_text())
    oracle = source_report["original"]["binary128"]
    if any(not isinstance(oracle[name], str) for name in ["heft_squared", "heft_ew_interference"]):
        raise ValueError("oracle outputs must retain decimal strings")
    exact = {
        "physical_s": str(s),
        "physical_t": str(t),
        "MH_squared": str(mass2),
        "coefficients": coefficients,
        "source": {"path": str(args.source.resolve()), "sha256": source_hash},
        "method": "Closed exact arithmetic AST over the four external scalar assignments; no binary float input or tensor operations",
    }
    provenance = {
        "source_sha256": source_hash,
        "coherent_input_sha256": hashlib.sha256(args.coherent_input.read_bytes()).hexdigest(),
        "oracle_report_sha256": hashlib.sha256(args.oracle_report.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.mkdir(parents=True)
    for filename, value in [("exact-heft.json", exact), ("heft-oracle.json", {"binary128": oracle}), ("provenance.json", provenance)]:
        (args.output / filename).write_text(json.dumps(value, indent=2) + "\n")


if __name__ == "__main__":
    main()
