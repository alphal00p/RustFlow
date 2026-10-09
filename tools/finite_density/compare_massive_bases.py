#!/usr/bin/env python3
"""Validate saved old/new massive-sunset predictions at one matching profile.

This comparison does not derive a basis transformation, generate predictions,
or read an oracle/reference answer. It compares already saved native outputs.
"""
import argparse
import json
from decimal import Decimal
from pathlib import Path

from compare_massive_reference import REPORT, ROOT, ZERO, complex_decimal, digest, norm, read

SECTORS = ("vacuum", "cut_0", "cut_1", "cut_01", "total")
TARGETS = ("scalar", "raised_numerator")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, default=REPORT / "old-basis-high-guard")
    parser.add_argument("--new", type=Path, default=REPORT / "full-sunset-frontier-sectors")
    parser.add_argument("--output", type=Path, default=REPORT / "full-sunset-basis-comparison.json")
    args = parser.parse_args()
    directories = [args.old.resolve(), args.new.resolve()]
    paths, records, basis_sizes = [], [], []
    for directory in directories:
        prediction_path = directory / "prediction-18-60-8.json"
        input_path = directory / "input.json"
        paths.extend([input_path, prediction_path])
        record = read(prediction_path)
        records.append(record)
        require(tuple(record[k] for k in ("digits", "guard_digits", "series_order", "occupied_start_scale")) == (18, 40, 60, 8), "prediction numerical settings differ")
        require(record["epsilon"] == "4/5", "prediction dimension mismatch")
        require(record["normalization"] == "unscaled Euclidean amplitude", "prediction normalization mismatch")
        require(record["independent_reference_comparisons"] == 0, "expected predictions saved before comparison")
        require(len(record["contributions"]) == 4, "all vacuum and occupied contributions are required")
        sizes = []
        for cut in ("0", "1", "01"):
            files = sorted((directory / f"cut-{cut}").glob("*-closed.json"))
            require(len(files) == 1, "expected one closed source checkpoint per occupied sector")
            checkpoint = read(files[0])
            require(checkpoint["status"] == "closed", "unclosed basis record")
            sizes.append(len(checkpoint["frontier"]))
            paths.extend([files[0], files[0].with_suffix(".bin")])
        basis_sizes.append(sizes)
    require(read(paths[0]) == read(directories[1] / "input.json"), "physical input definitions differ")
    require(records[0]["search_frontier_sectors"] is False and records[1]["search_frontier_sectors"] is True, "expected only native frontier search to distinguish the two basis constructions")
    for key in ("guard_refinement", "double_cut_source_policy", "double_cut_positive_compact_energy_powers"):
        require(records[0][key] == records[1][key], f"source or guard setting differs: {key}")
    comparisons = []
    for sector, old_entries, new_entries in zip(SECTORS, records[0]["contributions"] + [records[0]["assembled"]], records[1]["contributions"] + [records[1]["assembled"]]):
        require(len(old_entries) == len(new_entries) == 2, "target count mismatch")
        for target, old_text, new_text in zip(TARGETS, old_entries, new_entries):
            old, new = complex_decimal(old_text), complex_decimal(new_text)
            difference = tuple(b - a for a, b in zip(old, new))
            absolute = norm(difference)
            scale = max(norm(old), norm(new))
            relative = absolute / scale if scale else None
            small = scale < Decimal("1e-20")
            tolerance = Decimal("1e-25") if small else Decimal("1e-12")
            comparisons.append({
                "target": target, "sector": sector,
                "old_prediction": old_text, "new_prediction": new_text,
                "new_minus_old_real": str(difference[0]), "new_minus_old_imaginary": str(difference[1]),
                "absolute_difference": str(absolute), "relative_difference": str(relative) if relative is not None else None,
                "criterion": "absolute" if small else "relative", "tolerance": str(tolerance),
                "passed": absolute <= tolerance if small else relative <= tolerance,
                "imaginary_zero_passed": max(abs(old[1]), abs(new[1])) <= Decimal("1e-25"),
            })
    require(len(comparisons) == 10, "vacuum, all cuts and total must be compared for both targets")
    output = args.output.resolve()
    require(output not in paths, "comparison output may not replace an input artifact")
    passed = all(c["passed"] and c["imaginary_zero_passed"] for c in comparisons)
    report = {
        "schema": 1, "status": "passed" if passed else "failed",
        "scope": "Numerical agreement of already saved native predictions at one shared regulator/precision/order/start setting; no exact master transformation or new independent-reference comparison.",
        "epsilon": "4/5", "normalization": "unscaled Euclidean amplitude",
        "digits": 18, "guard_digits": 40, "series_order": 60, "occupied_start_scale": 8,
        "old_basis_sizes": basis_sizes[0], "new_basis_sizes": basis_sizes[1],
        "old_search_frontier_sectors": False, "new_search_frontier_sectors": True,
        "comparison_count": len(comparisons), "independent_reference_values_read": 0,
        "supplied_oracle_numerical_records_compared": 0,
        "arithmetic": "Python Decimal, 100 decimal digits, original saved decimal strings",
        "inputs_sha256": {label(path): digest(path) for path in paths},
        "script_sha256": digest(Path(__file__)),
        "arithmetic_helper_sha256": digest(Path(__file__).with_name("compare_massive_reference.py")),
        "largest_relative_difference": str(max(Decimal(c["relative_difference"]) for c in comparisons if c["relative_difference"] is not None)),
        "comparisons": comparisons,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("status", "old_basis_sizes", "new_basis_sizes", "comparison_count", "largest_relative_difference")}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
