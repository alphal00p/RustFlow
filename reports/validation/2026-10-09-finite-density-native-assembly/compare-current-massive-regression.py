"""Compare one saved runtime profile; this is not a new refinement study."""
import argparse
import json
import sys
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/finite_density"))
from compare_complete_massive_reference import (
    NORMALIZATION, SECTORS, TARGETS, compare, component_keys, require,
)
from compare_massive_reference import ZERO, complex_decimal, digest, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report_dir = Path(__file__).resolve().parent
    directory = args.predictions.resolve()
    current_path = directory / "prediction-28-80-12.json"
    input_path = directory / "input.json"
    reference_path = report_dir / "independent-reference/quadrature-64-60.json"
    provenance_path = report_dir / "independent-reference/independent-reference.json"
    previous_path = report_dir / "massless-blocks-massive-regression/prediction-28-80-12.json"
    previous_input = previous_path.parent / "input.json"
    paths = [current_path, input_path, reference_path, provenance_path,
             previous_path, previous_input, Path(__file__).resolve()]
    require(args.output.resolve() not in paths, "output must not replace input")
    current, reference, provenance, previous = map(
        read, [current_path, reference_path, provenance_path, previous_path]
    )
    require(read(input_path) == read(previous_input) == provenance["definition"],
            "definitions differ")
    require(provenance["status"] == "refinement_passed", "reference incomplete")
    require((reference["order"], reference["digits"]) == (64, 60), "reference settings")
    require(current["epsilon"] == provenance["epsilon"] == "4/5", "dimension mismatch")
    require(current["normalization"] == provenance["normalization"] == NORMALIZATION,
            "normalization mismatch")
    require(current["full_amplitude"] is True, "full amplitude required")
    require(current["independent_reference_comparisons"] == 0, "predictions first")
    require([current[k] for k in ("digits", "series_order", "occupied_start_scale")]
            == [28, 80, 12], "unexpected profile")
    require(current["guard_digits"] == 40, "unexpected guard precision")
    require([c["cut_slots"] for c in current["contributions"]]
            == [[], [0], [1], [0, 1]], "incomplete or reordered cut set")
    now = [c["values"] for c in current["contributions"]] + [current["values"]]
    before = previous["contributions"] + [previous["assembled"]]
    require(all(len(v) == len(TARGETS) for v in now + before), "target count")
    comparisons, regressions, assembly = [], [], []
    for sector, values, old_values in zip(SECTORS, now, before):
        for target, text, old_text in zip(TARGETS, values, old_values):
            value = complex_decimal(text)
            keys = component_keys(target, sector)
            expected = sum((Decimal(reference["contributions"][k]) for k in keys), ZERO)
            comparisons.append({"sector": sector, "target": target,
                                "reference_component_keys": keys,
                                "prediction": text, "reference_real": str(expected),
                                **compare(value, (expected, ZERO), "sample"),
                                "imaginary_zero_passed": abs(value[1]) <= Decimal("1e-25")})
            regressions.append({"sector": sector, "target": target,
                                "identical_decimal_string": text == old_text,
                                **compare(value, complex_decimal(old_text), "sample")})
    for index, target in enumerate(TARGETS):
        total = tuple(sum((complex_decimal(row[index])[part] for row in now[:-1]), ZERO)
                      for part in range(2))
        assembly.append({"target": target,
                         **compare(total, complex_decimal(now[-1][index]), "sample")})
    passed = (all(c["passed"] and c["imaginary_zero_passed"] for c in comparisons)
              and all(c["passed"] for c in regressions + assembly))
    result = {"schema": 1, "status": "passed" if passed else "failed",
              "scope": "one current-source D12/5 profile, both original targets, vacuum/all cuts/full amplitude; independent reference and same-profile historical regression",
              "full_amplitude": True, "configuration": [28, 80, 12], "guard_digits": 40,
              "reference_comparison_count": len(comparisons),
              "historical_same_profile_comparison_count": len(regressions),
              "independent_refinement_comparison_count": 0,
              "supplied_oracle_numerical_records_compared": 0,
              "reference_error_estimate": provenance["error_estimate"],
              "criterion": "relative <=1e-12 above reference magnitude1e-20; absolute <=1e-25 otherwise; imaginary absolute <=1e-25",
              "inputs_sha256": {str(p.relative_to(ROOT)): digest(p) for p in paths},
              "reference_comparisons": comparisons, "historical_regressions": regressions,
              "assembly_checks": assembly}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("status", "reference_comparison_count",
                                           "historical_same_profile_comparison_count")}))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
