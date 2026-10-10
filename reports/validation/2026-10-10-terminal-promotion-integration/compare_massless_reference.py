#!/usr/bin/env python3
"""Compare completed native massless predictions to separately generated references.

Reads saved artifacts only. Requires four sample or five Laurent profiles, with
independent precision/order/start/grid changes. No production solver is imported.
"""
import argparse
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/"tools/finite_density"))
from certificate_validation import load_certificates, validate_sample
from compare_massive_reference import (
    ROOT, ZERO, CONFIGURATIONS, complex_decimal, difference, digest, read,
)

TARGETS = ("scalar", "raised_original_numerator")
CUTS = ((), (0,), (1,), (0, 1))
SECTORS = ("vacuum", "single_cut_0", "single_cut_1", "double_cut", "total")
SETTINGS = ("digits", "series_order", "occupied_start_scale", "epsilon_grid_denominator")


def require(value, message):
    if not value:
        raise ValueError(message)


def label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def value(text):
    return (ZERO, ZERO) if text == "0" else complex_decimal(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("sample", "laurent"))
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--additional-predictions", type=Path, action="append", default=[])
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    laurent = args.mode == "laurent"
    dirs = [args.predictions.resolve()] + [p.resolve() for p in args.additional_predictions]
    require(len(set(dirs)) == len(dirs), "duplicate prediction directory")
    refpath = args.reference.resolve()
    reference = read(refpath)
    require(reference["status"] == "independent_reference_generated", "reference generation incomplete")
    require(reference["normalization"] == "raw_euclidean_dDp_over_2pi_power_D_per_loop", "reference normalization mismatch")
    require(reference["oracle_records_read"] == 0 and reference["native_predictions_read"] == 0, "reference must be independent")
    definition_path = ROOT / reference["source"]["input"]
    require(digest(definition_path) == reference["source"]["input_sha256"], "reference input digest changed")
    definition = read(definition_path)
    require(definition["name"] == "massless_two_loop_sunset", "unexpected reference definition")
    input_paths = [d / "input.json" for d in dirs]
    require(all(read(p) == definition for p in input_paths), "prediction definition differs from independent reference")
    configs = tuple((*c, 1000) for c in CONFIGURATIONS) + ((28, 80, 12, 2000),) if laurent else tuple(CONFIGURATIONS)
    paths = []
    duplicate_paths = []
    for config in configs:
        digits, order, start = config[:3]
        suffix = f"-grid-{config[3]}" if laurent else ""
        name = f"prediction-{digits}-{order}-{start}{suffix}.json"
        candidates = [d / name for d in dirs if (d / name).is_file()]
        require(candidates, f"missing completed native profile {config}")
        require(all(digest(p) == digest(candidates[0]) for p in candidates), "conflicting duplicated profile")
        paths.append(candidates[0]); duplicate_paths.extend(candidates[1:])
    inputs = [refpath, definition_path] + input_paths + paths + duplicate_paths
    output = args.output.resolve()
    require(output not in inputs, "output cannot replace an input artifact")
    certificate_arrays=[]
    for directory in dirs:
        certificates, certificate_path=load_certificates(directory,definition,CUTS)
        certificate_arrays.append(certificates); inputs.append(certificate_path)
    require(all(c==certificate_arrays[0]for c in certificate_arrays),"physical-zero certificates changed across resumed runs")
    inputs.append(Path(__file__).with_name("certificate_validation.py"))
    previous = None
    guard = None
    dimension = None
    comparisons = []
    refinements = []
    assembly = []
    for profile_index, (config, path) in enumerate(zip(configs, paths)):
        record = read(path)
        require(record["full_amplitude"] is True, "complete amplitude required")
        require(record["normalization"] == "unscaled Euclidean amplitude", "native normalization mismatch")
        require(record["independent_reference_comparisons"] == 0, "native outputs must precede comparison")
        require(tuple(record[k] for k in SETTINGS[:3]) == config[:3], "native configuration mismatch")
        require(type(record["guard_digits"]) is int and record["guard_digits"] >= 0, "invalid guard precision")
        guard = record["guard_digits"] if guard is None else guard
        require(record["guard_digits"] == guard, "guard precision changed between profiles")
        values = {}
        expected = {}
        if laurent:
            require(record["epsilon_grid_denominator"] == config[3], "epsilon-grid mismatch")
            require(len(record["expansions"]) == len(TARGETS), "target count mismatch")
            require(reference["laurent"]["orders"] == [-2, -1, 0], "reference coefficient orders mismatch")
            require(reference["laurent"]["chemical_potential"] == definition["chemical_potentials"][0], "reference chemical potential mismatch")
            for i, expansion in enumerate(record["expansions"]):
                require(type(expansion["verified_digits"]) is int and expansion["verified_digits"] >= config[0], "shared Laurent fitter did not verify requested digits")
                require(set(expansion["coefficients"]) == {"-2", "-1", "0"}, "requested coefficient set mismatch")
                for j, power in enumerate((-2, -1, 0)):
                    key = (TARGETS[i], "total", power)
                    values[key] = value(expansion["coefficients"][str(power)])
                    expected[key] = value(reference["laurent"]["total"][i][j])
        else:
            d = Fraction(4) - 2 * Fraction(record["epsilon"])
            dimension = d if dimension is None else dimension
            require(d == dimension, "dimension changed between profiles")
            candidates = [s for s in reference["samples"] if Fraction(s["dimension"]) == d and Fraction(s["chemical_potential"]) == Fraction(definition["chemical_potentials"][0])]
            require(len(candidates) == 1, "no unique independent reference for sample")
            sample = candidates[0]
            require(tuple(tuple(c["cut_slots"]) for c in record["contributions"]) == CUTS, "complete physical cut set mismatch")
            occupied = record["occupied_reports"]
            require(tuple(tuple(c["cut_slots"]) for c in occupied) == CUTS[1:], "occupied provenance does not cover all cuts")
            validate_sample(record, certificate_arrays[dirs.index(path.parent)], definition, CUTS,
                            read(path.parent/"configuration.json")["source_options"])
            for sector, entries in zip(SECTORS, [c["values"] for c in record["contributions"]] + [record["values"]]):
                require(len(entries) == len(TARGETS), "target count mismatch")
                for i, target in enumerate(TARGETS):
                    key = (target, sector, None)
                    values[key] = value(entries[i]); expected[key] = value(sample[sector][i])
            for target in TARGETS:
                summed = tuple(sum((values[(target, sector, None)][part] for sector in SECTORS[:-1]), ZERO) for part in range(2))
                assembly.append({"target": target, "configuration": config, **difference(summed, values[(target, "total", None)])})
        for key, actual in values.items():
            target, sector, power = key
            comparisons.append({"target": target, "sector": sector, "power": power, "configuration": config, "prediction_path": label(path), "prediction_real": str(actual[0]), "prediction_imaginary": str(actual[1]), "reference_real": str(expected[key][0]), **difference(actual, expected[key]), "imaginary_zero_passed": abs(actual[1]) <= Decimal("1e-25")})
            if previous is not None:
                refinements.append({"target": target, "sector": sector, "power": power, "from_configuration": configs[profile_index-1], "to_configuration": config, "varied_setting": SETTINGS[profile_index-1], **difference(actual, previous[key])})
        previous = values
    passed = all(x["passed"] and x["imaginary_zero_passed"] for x in comparisons) and all(x["passed"] for x in refinements + assembly)
    report = {"schema_version": 1, "status": "passed" if passed else "failed", "mode": args.mode, "scope": "complete massless two-loop sunset: separate certified physical singleton zeros plus native AMF double cut; no four-loop acceptance", "normalization": "unscaled Euclidean amplitude", "dimension": str(dimension) if dimension else "4-2*epsilon", "configurations": configs, "guard_digits": guard, "source_policy_provenance": "coherent runtime source snapshot; certified singleton physical zeros are distinct from the double-cut native source policy and AMF closure", "relative_tolerance": "1e-12", "small_magnitude_threshold": "1e-20", "small_absolute_tolerance": "1e-25", "imaginary_zero_absolute_tolerance": "1e-25", "reference_comparison_count": len(comparisons), "independent_refinement_comparison_count": len(refinements), "supplied_oracle_numerical_records_compared": 0, "reference_uncertainty": reference["uncertainty"], "arithmetic": "Python Decimal, 100 digits, saved decimal inputs", "inputs_sha256": {label(p): digest(p) for p in inputs}, "script_sha256": digest(Path(__file__)), "comparisons": comparisons, "refinements": refinements, "assembly_checks": assembly}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "mode", "reference_comparison_count", "independent_refinement_comparison_count")}))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
