#!/usr/bin/env python3
"""Compare saved complete native predictions, after generation, to independent data.

This validation-only program cannot generate or modify solver predictions. It
requires all four sample precision/order/start configurations, and a fifth
independent epsilon-grid configuration for Laurent data. All inputs are hashed.
"""

import argparse
import json
from decimal import Decimal
from pathlib import Path

from compare_massive_reference import (
    CONFIGURATIONS, REPORT, ROOT, ZERO, complex_decimal, digest, norm, read,
)

TARGETS = ("scalar", "raised_numerator")
SECTORS = ("vacuum", "cut_0", "cut_1", "cut_01", "total")
NORMALIZATION = "unscaled Euclidean amplitude"
LAURENT_CONFIGURATIONS = tuple((*config, 1000) for config in CONFIGURATIONS) + (
    (28, 80, 12, 2000),
)
VARIED_SETTINGS = ("digits", "series_order", "occupied_start_scale", "epsilon_grid_denominator")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def compare(value, reference, mode):
    absolute = norm(tuple(a - b for a, b in zip(value, reference)))
    scale = norm(reference)
    relative = absolute / scale if scale else None
    threshold = Decimal("1e-18" if mode == "laurent" else "1e-20")
    absolute_tolerance = Decimal("1e-20" if mode == "laurent" else "1e-25")
    small = scale < threshold
    passed = absolute <= absolute_tolerance if small else relative <= Decimal("1e-12")
    return {
        "absolute_difference": str(absolute),
        "relative_difference": str(relative) if relative is not None else None,
        "criterion": "absolute" if small else "relative",
        "tolerance": str(absolute_tolerance if small else Decimal("1e-12")),
        "passed": passed,
    }


def component_keys(target, sector):
    base = f"{target}/{sector}"
    if target == "raised_numerator" and sector in ("cut_0", "cut_01"):
        return [base + "_bulk", base + "_surface"]
    return [base]


def prediction_filename(configuration):
    digits, order, start = configuration[:3]
    grid = configuration[3] if len(configuration) == 4 else None
    suffix = f"-grid-{grid}" if grid is not None and grid != 1000 else ""
    return f"prediction-{digits}-{order}-{start}{suffix}.json"


def invariant_profile(record):
    guard_digits = record["guard_digits"]
    require(type(guard_digits) is int and 0 <= guard_digits < 2**32,
            "guard_digits must be an unsigned 32-bit integer")
    search_frontier = record["search_frontier_sectors"]
    require(type(search_frontier) is bool, "search_frontier_sectors must be a boolean")
    return {"guard_digits": guard_digits, "search_frontier_sectors": search_frontier}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("sample", "laurent"))
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prediction_dir = args.predictions.resolve()
    output = args.output.resolve()
    laurent = args.mode == "laurent"
    reference_dir = REPORT / ("independent-laurent-reference" if laurent else "independent-reference")
    provenance_path = reference_dir / ("independent-laurent-reference.json" if laurent else "independent-reference.json")
    reference_path = reference_dir / ("laurent-64-80-20000.json" if laurent else "quadrature-64-60.json")
    input_path = prediction_dir / "input.json"
    configurations = LAURENT_CONFIGURATIONS if laurent else CONFIGURATIONS
    paths = [input_path, provenance_path, reference_path] + [
        prediction_dir / prediction_filename(configuration)
        for configuration in configurations
    ]
    require(output not in [p.resolve() for p in paths], "output must not replace any input artifact")
    provenance = read(provenance_path)
    reference = read(reference_path)
    require(read(input_path) == provenance["definition"], "input definition differs from the reference")
    require(provenance["status"] == "refinement_passed", "reference refinement has not passed")
    require(provenance["normalization"] == NORMALIZATION, "unexpected reference normalization")
    if laurent:
        require((reference["nodes"], reference["digits"], reference["grid_denominator"]) == (64, 80, 20000), "unexpected Laurent reference configuration")
    else:
        require((reference["order"], reference["digits"]) == (64, 60), "unexpected sample reference configuration")
        require(provenance["epsilon"] == "4/5", "unexpected reference dimension")

    comparisons, refinements, assembly_checks = [], [], []
    previous = None
    common_profile = None
    for config_index, (configuration, path) in enumerate(zip(configurations, paths[3:])):
        digits, order, start = configuration[:3]
        record = read(path)
        require(record["normalization"] == NORMALIZATION, "unexpected prediction normalization")
        require((record["digits"], record["series_order"], record["occupied_start_scale"]) == (digits, order, start), "prediction configuration mismatch")
        require(record["independent_reference_comparisons"] == 0, "expected predictions saved before comparison")
        profile = invariant_profile(record)
        if common_profile is None:
            common_profile = profile
        require(profile == common_profile,
                "guard_digits and search_frontier_sectors must remain constant across comparison configurations")
        values = {}
        if laurent:
            require(type(record["epsilon_grid_denominator"]) is int
                    and record["epsilon_grid_denominator"] == configuration[3],
                    "prediction epsilon-grid denominator mismatch")
            require(record["full_amplitude"] is True, "Laurent data must represent the complete amplitude")
            require(len(record["expansions"]) == len(TARGETS), "target count mismatch")
            for target, expansion in zip(TARGETS, record["expansions"]):
                require(expansion["verified_digits"] >= digits, "production Laurent refinement failed")
                require(set(expansion["coefficients"]) == {"-2", "-1", "0"}, "requested Laurent orders missing or unexpected")
                for power, value in expansion["coefficients"].items():
                    values[(target, "total", power)] = complex_decimal(value)
        else:
            require(record["epsilon"] == "4/5", "prediction dimension mismatch")
            require(len(record["contributions"]) == 4, "all four vacuum/cut contributions required")
            for sector, entries in zip(SECTORS, record["contributions"] + [record["assembled"]]):
                require(len(entries) == len(TARGETS), "target count mismatch")
                for target, value in zip(TARGETS, entries):
                    values[(target, sector, None)] = complex_decimal(value)
            for target in TARGETS:
                total = tuple(sum((values[(target, sector, None)][part] for sector in SECTORS[:-1]), ZERO) for part in range(2))
                assembly_checks.append({"target": target, "configuration": list(configuration), **compare(total, values[(target, "total", None)], args.mode)})
        for (target, sector, power), value in values.items():
            keys = component_keys(target, sector)
            raw = reference["coefficients"] if laurent else reference["contributions"]
            expected = sum((Decimal(raw[key][power] if laurent else raw[key]) for key in keys), ZERO)
            absolute_tolerance = Decimal("1e-20" if laurent else "1e-25")
            comparisons.append({
                "target": target, "sector": sector, "power": power,
                "configuration": list(configuration), "prediction_path": label(path),
                "prediction_real": str(value[0]), "prediction_imaginary": str(value[1]),
                "reference_real": str(expected), "reference_component_keys": keys,
                **compare(value, (expected, ZERO), args.mode),
                "imaginary_zero_passed": abs(value[1]) <= absolute_tolerance,
                "imaginary_zero_absolute_tolerance": str(absolute_tolerance),
            })
            if previous is not None:
                refinements.append({
                    "target": target, "sector": sector, "power": power,
                    "varied_setting": VARIED_SETTINGS[config_index - 1],
                    "from_configuration": configurations[config_index - 1],
                    "to_configuration": list(configuration),
                    **compare(previous[(target, sector, power)], value, args.mode),
                })
        previous = values

    passed = all(c["passed"] and c["imaginary_zero_passed"] for c in comparisons) and all(c["passed"] for c in refinements + assembly_checks)
    report = {
        "schema": 2, "status": "passed" if passed else "failed", "mode": args.mode,
        "full_amplitude": True, "normalization": NORMALIZATION,
        "scope": "massive nonfactorized sunset, scalar and raised original medium numerator",
        "epsilon": None if laurent else "4/5", "laurent_orders": [-2, 0] if laurent else None,
        "configuration_fields": list(VARIED_SETTINGS[:4 if laurent else 3]),
        "configurations": configurations,
        "invariant_profile": common_profile,
        "reference_comparison_count": len(comparisons),
        "independent_refinement_comparison_count": len(refinements),
        "supplied_oracle_numerical_records_compared": 0,
        "reference_error_estimate": provenance["error_estimate"],
        "reference_caution": "observed refinement agreement is not a rigorous interval error bound",
        "arithmetic": "Python Decimal, 100 decimal digits, saved decimal inputs",
        "inputs_sha256": {label(path): digest(path) for path in paths},
        "script_sha256": digest(Path(__file__)),
        "comparisons": comparisons, "refinements": refinements, "assembly_checks": assembly_checks,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("status", "mode", "reference_comparison_count", "independent_refinement_comparison_count")}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
