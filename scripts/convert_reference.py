#!/usr/bin/env python3
"""Convert the pinned, simple upstream Laurent fixture, retaining precision.

This parser accepts only the numeric grammar present in this fixture. It does
not execute Mathematica or infer accuracy from printed mantissa length.
"""
import hashlib
import json
import re
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = root / "fixtures/amflow-2.0/two_loop_blade_auto.wl"
raw = source.read_bytes()
text = raw.decode()
records = list(re.finditer(r"j\[tt,([^\]]+)\]\s*->", text))
result = []
number = re.compile(r"([+-]?)(\d+\.\d+)`(\d+(?:\.\d*)?)(\*I)?")
for index, match in enumerate(records):
    powers = [int(s.strip()) for s in match[1].split(",")]
    expression = text[match.end():records[index + 1].start() if index + 1 < len(records) else len(text)]
    expression = re.sub(r"\s", "", expression).rstrip(",}")
    depth, starts = 0, [0]
    for pos, char in enumerate(expression):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char in "+-" and depth == 0 and pos > 0:
            starts.append(pos)
    starts.append(len(expression))
    coefficients = {}
    for begin, end in zip(starts, starts[1:]):
        term = expression[begin:end]
        sign = -1 if term.startswith("-") else 1
        term = term.lstrip("+-")
        parts = term.split("/eps")
        exponent = 0 if len(parts) == 1 else -int(parts[1][1:] or "1") if parts[1] else -1
        coefficient = parts[0].strip("()")
        components = {"re": "0", "im": "0", "re_precision": None, "im_precision": None}
        position = 0
        for value in number.finditer(coefficient):
            if value.start() != position:
                raise ValueError(f"Unrecognized coefficient: {coefficient}")
            position = value.end()
            field = "im" if value[4] else "re"
            if components[field + "_precision"] is not None:
                raise ValueError("Repeated component requires precision propagation")
            negative = sign * (-1 if value[1] == "-" else 1) < 0
            components[field] = ("-" if negative else "") + value[2]
            components[field + "_precision"] = value[3]
        if position != len(coefficient):
            raise ValueError(f"Unparsed coefficient: {coefficient}")
        coefficients[str(exponent)] = components
    if set(coefficients) != {str(k) for k in range(-4, 1)}:
        raise ValueError("Unexpected Laurent powers")
    result.append({"indices": powers, "coefficients": coefficients})
output = {
    "schema": 1,
    "upstream_commit": "26005517a288086c4cb4d1b26d829691bc088485",
    "upstream_path": "examples/automatic_vs_manual/backup/blade_auto",
    "source_sha256": hashlib.sha256(raw).hexdigest(),
    "precision_semantics": "Mathematica component decimal precision; null denotes exact zero",
    "targets": result,
}
(root / "fixtures/amflow-2.0/two_loop.json").write_text(json.dumps(output, indent=2) + "\n")
