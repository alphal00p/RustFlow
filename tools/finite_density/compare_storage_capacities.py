#!/usr/bin/env python3
"""Compare saved physical occupied-flow predictions across native capacities.

No reference answers, integral evaluation, or algebraic reduction is performed.
"""
import argparse
import json
from decimal import Decimal
from pathlib import Path

from compare_massive_reference import REPORT, ROOT, complex_decimal, digest, norm, read


def require(condition, message):
    if not condition:
        raise ValueError(message)


def label(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=REPORT / "native-capacity-physical")
    parser.add_argument("--output", type=Path, default=REPORT / "native-capacity-physical-comparison.json")
    args = parser.parse_args()
    root = args.predictions.resolve()
    provenance = [root / "input.json"]
    read(provenance[0])
    comparisons, connections = [], []
    for sector, physical, cuts in [("cut-0", 7, [0]), ("cut-01", 9, [0, 1])]:
        records, exact = [], []
        for capacity in (physical, 12):
            path = root / sector / f"prediction-capacity-{capacity}.json"
            connection = root / sector / f"connection-capacity-{capacity}.json"
            provenance.extend([path, connection])
            record, system = read(path), read(connection)
            require(record["physical_arity"] == system["physical_arity"] == physical, "physical arity mismatch")
            require(record["native_storage_capacity"] == system["native_storage_capacity"] == capacity, "native capacity mismatch")
            require(record["cut_slots"] == cuts and record["full_amplitude"] is False, "sector identity mismatch")
            require(record["epsilon"] == "4/5", "regulator mismatch")
            require(record["normalization"] == "unscaled Euclidean amplitude", "normalization mismatch")
            require(tuple(record[k] for k in ("digits", "guard_digits", "series_order", "occupied_start_scale")) == (18, 40, 60, 8), "numerical profile mismatch")
            require(record["independent_reference_comparisons"] == 0, "expected predictions saved before comparison")
            require(len(record["values"]) == 2, "target count mismatch")
            require(all(len(index) == physical for index in system["basis"]), "storage tail escaped into physical basis")
            require(all(len(term["indices"]) == physical for target in system["targets"] for term in target), "storage tail escaped into target weights")
            checkpoints = list((root / sector / f"native-capacity-{capacity}").glob("*-closed.json"))
            require(len(checkpoints) == 1, "expected one closed native checkpoint")
            checkpoint = read(checkpoints[0])
            require(checkpoint["physical_arity"] == physical and checkpoint["storage_capacity"] == capacity, "checkpoint layout mismatch")
            require(all(len(index) == capacity and all(i == 0 for i in index[physical:]) for index in checkpoint["frontier"]), "nonzero native storage tail")
            provenance.extend([checkpoints[0], checkpoints[0].with_suffix(".bin")])
            records.append(record)
            exact.append(system)
        connections.append({"sector": sector, "basis_sizes": [len(c["basis"]) for c in exact],
                            "identical_physical_connection": all(exact[0][key] == exact[1][key] for key in ("basis", "matrix", "targets", "nonzero_conditions"))})
        for target, (left, right) in enumerate(zip(records[0]["values"], records[1]["values"])):
            a, b = complex_decimal(left), complex_decimal(right)
            error = norm((a[0] - b[0], a[1] - b[1]))
            scale = max(norm(a), norm(b))
            tolerance = Decimal("1e-25") if scale < Decimal("1e-20") else scale * Decimal("1e-15")
            comparisons.append({"sector": sector, "target": target, "absolute_error": str(error),
                                "relative_error": str(error / scale) if scale else None,
                                "tolerance": str(tolerance), "passed": error <= tolerance})
    result = {"schema": 1, "passed": all(c["passed"] for c in comparisons),
              "scope": "Two occupied sectors, both scalar and raised medium-numerator targets; exact versus larger native storage; no independent reference or four-loop numerical comparison.",
              "connections": connections, "comparisons": comparisons,
              "sha256": {label(path): digest(path) for path in provenance}}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed": result["passed"], "comparisons": len(comparisons), "connections": connections}))
    require(result["passed"], "capacity predictions differ beyond tolerance")


if __name__ == "__main__":
    main()
