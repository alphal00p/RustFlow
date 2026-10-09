#!/usr/bin/env python3
"""Validation only: compare frozen AMF predictions with independent quadrature.

This utility never imports the production solver or changes either input artifact.
Run from any directory; the report directory is relative to the repository root.
"""

import hashlib
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly"
CONFIGURATIONS = [(18, 60, 8), (28, 60, 8), (28, 80, 8), (28, 80, 12)]
getcontext().prec = 100
ZERO = Decimal(0)
RELATIVE_TOLERANCE = Decimal("1e-12")
SMALL_THRESHOLD = Decimal("1e-20")
ABSOLUTE_TOLERANCE = Decimal("1e-25")
COMPLEX = re.compile(r"^\(([-+0-9.eE]+)\+([-+0-9.eE]+)i\)$")


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_decimal(value):
    match = COMPLEX.fullmatch(value)
    assert match is not None, f"unsupported saved complex literal: {value}"
    return tuple(Decimal(part) for part in match.groups())


def norm(value):
    return sum((part * part for part in value), ZERO).sqrt()


def difference(value, reference):
    delta = tuple(a - b for a, b in zip(value, reference))
    absolute = norm(delta)
    scale = norm(reference)
    relative = absolute / scale if scale else None
    small = scale < SMALL_THRESHOLD
    passed = absolute <= ABSOLUTE_TOLERANCE if small else relative <= RELATIVE_TOLERANCE
    return {
        "real_difference": str(delta[0]),
        "imaginary_difference": str(delta[1]),
        "absolute_difference": str(absolute),
        "relative_difference": str(relative) if relative is not None else None,
        "criterion": "absolute" if small else "relative",
        "tolerance": str(ABSOLUTE_TOLERANCE if small else RELATIVE_TOLERANCE),
        "passed": passed,
    }


def main():
    reference_path = REPORT / "independent-reference/quadrature-64-60.json"
    provenance_path = REPORT / "independent-reference/independent-reference.json"
    input_path = REPORT / "single-sunset-attempt-1/input.json"
    reference = read(reference_path)
    provenance = read(provenance_path)
    assert read(input_path) == provenance["definition"]
    assert (reference["order"], reference["digits"]) == (64, 60)
    assert provenance["status"] == "refinement_passed"
    assert provenance["epsilon"] == "4/5"
    assert provenance["normalization"] == "unscaled Euclidean amplitude"
    components = reference["contributions"]
    keys = {
        ("vacuum", "scalar"): ["scalar/vacuum"],
        ("vacuum", "raised_numerator"): ["raised_numerator/vacuum"],
        ("cut-0", "scalar"): ["scalar/cut_0"],
        ("cut-0", "raised_numerator"): [
            "raised_numerator/cut_0_bulk", "raised_numerator/cut_0_surface"
        ],
    }
    references = {
        key: (sum((Decimal(components[k]) for k in parts), ZERO), ZERO)
        for key, parts in keys.items()
    }
    comparisons = []
    refinements = []
    inputs = {
        str(p.relative_to(ROOT)): digest(p)
        for p in [reference_path, provenance_path, input_path]
    }
    for sector in ["vacuum", "cut-0"]:
        previous = None
        for config_index, (digits, order, start) in enumerate(CONFIGURATIONS):
            path = REPORT / f"single-sunset-attempt-1/prediction-{sector}-{digits}-{order}-{start}.json"
            record = read(path)
            assert record["sector"] == sector
            assert record["epsilon"] == provenance["epsilon"]
            assert record["normalization"] == provenance["normalization"]
            assert not record["full_amplitude"]
            assert (record["digits"], record["series_order"], record["occupied_start_scale"]) == (digits, order, start)
            assert len(record["values"]) == 2
            relative_path = str(path.relative_to(ROOT))
            inputs[relative_path] = digest(path)
            values = list(map(complex_decimal, record["values"]))
            for index, target in enumerate(["scalar", "raised_numerator"]):
                key = (sector, target)
                comparisons.append({
                    "sector": sector, "target": target,
                    "configuration": [digits, order, start],
                    "prediction_path": relative_path,
                    "prediction": record["values"][index],
                    "reference_component_keys": keys[key],
                    "reference_real": str(references[key][0]),
                    **difference(values[index], references[key]),
                    "imaginary_zero_absolute_difference": str(abs(values[index][1])),
                    "imaginary_zero_tolerance": str(ABSOLUTE_TOLERANCE),
                    "imaginary_zero_passed": abs(values[index][1]) <= ABSOLUTE_TOLERANCE,
                })
                if previous is not None:
                    refinements.append({
                        "sector": sector, "target": target,
                        "varied_setting": ["digits", "series_order", "occupied_start_scale"][config_index - 1],
                        "from_configuration": CONFIGURATIONS[config_index - 1],
                        "to_configuration": [digits, order, start],
                        **difference(previous[index], values[index]),
                    })
            previous = values
    passed = all(c["passed"] and c["imaginary_zero_passed"] for c in comparisons) and all(c["passed"] for c in refinements)
    report = {
        "schema": 1, "status": "passed" if passed else "failed",
        "scope": "massive sunset vacuum and cut-0 sectors only, two original targets at D=12/5",
        "full_amplitude": False, "laurent_coefficients_compared": 0,
        "supplied_oracle_numerical_records_compared": 0,
        "unique_reference_values_compared": 4,
        "reference_comparison_count": len(comparisons),
        "independent_refinement_comparison_count": len(refinements),
        "arithmetic": "Python Decimal, 100 decimal digits; exact saved decimal inputs, rounded square roots",
        "criterion": "complex absolute difference <= 1e-25 when reference magnitude < 1e-20; otherwise relative difference <= 1e-12",
        "reference_error_estimate": provenance["error_estimate"],
        "reference_caution": "agreement below the quadrature refinement changes is observed agreement, not a proven accuracy bound",
        "normalization": provenance["normalization"], "epsilon": provenance["epsilon"],
        "inputs_sha256": inputs,
        "script_sha256": digest(Path(__file__)),
        "comparisons": comparisons, "refinements": refinements,
        "largest_reference_relative_difference": str(max(Decimal(c["relative_difference"]) for c in comparisons)),
        "largest_refinement_relative_difference": str(max(Decimal(c["relative_difference"]) for c in refinements)),
    }
    output = REPORT / "single-sunset-reference-comparison.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["status", "reference_comparison_count", "independent_refinement_comparison_count", "largest_reference_relative_difference", "largest_refinement_relative_difference"]}, indent=2))
    for comparison in comparisons:
        if comparison["configuration"] == [28, 80, 12]:
            print(comparison["sector"], comparison["target"], "absolute", comparison["absolute_difference"], "relative", comparison["relative_difference"])
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
