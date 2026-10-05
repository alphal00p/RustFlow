#!/usr/bin/env python3
"""Compare a completed native physical-source Laurent fit with the plugin record.

This program is deliberately separate from boundary generation. It opens the
independent reference only after checking that the native run succeeded. All
acceptance inequalities use exact fractions of the serialized decimal numbers;
Decimal square roots are used only to print diagnostics. No epsilon sample is
compared with a truncated Laurent series and no numerical seed is produced.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / "fixtures/gg-hg/full-systems.json"
PHYSICAL_MAP_SOURCES = {
    "Planar_EW1": "a6bd80f1edeb056f89ee1613561be3066bb3b3ea1fe7f4ab65b928891fa74380",
    "NP_EW1": "4a65e20f6244405152a8d87290daca63bb3a80aca7f63e024efa54cb9d8dcccd",
}
NORMALIZATION = (
    "arXiv:2112.07578v1 canonical F vector; mV^2=1; "
    "exp(2*EulerGamma*epsilon) relative to ordinary two-loop integrals; "
    "physical +i0 and the recorded individual root germs; "
    "no form-factor -1/(mV^2)^2/(4*pi)^4 prefactor"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    require(isinstance(value, str), "numerical coefficients must be exact decimal strings")
    return Fraction(value)


def norm_squared(value):
    return rational(value["real"]) ** 2 + rational(value["imaginary"]) ** 2


def difference_squared(left, right):
    return sum(
        (rational(left[part]) - rational(right[part])) ** 2
        for part in ("real", "imaginary")
    )


def diagnostic_sqrt(value):
    with localcontext() as context:
        context.prec = 100
        return str((Decimal(value.numerator) / Decimal(value.denominator)).sqrt())


def read_json(path):
    data = Path(path).read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def scoped_symbols(values, namespace):
    result = {}
    for symbol, value in values.items():
        prefix = f"{namespace}::{{}}::"
        require(symbol.startswith(prefix), "native coordinate or root has an unexpected namespace")
        name = symbol[len(prefix):]
        require(name and "::" not in name and name not in result, "ambiguous scoped symbol")
        result[name] = value
    return result


def validate_native(native):
    require(native.get("schema") == "rustflow-native-gg-hg-boundary-run-v1", "unsupported native report schema")
    require(native.get("status") == "success", "native generation has not completed successfully; references were not opened")
    require(native.get("verified_bank_saved") is True, "native report has no saved verified boundary")
    require(native.get("family") in PHYSICAL_MAP_SOURCES, "unsupported native physical family")
    require(native.get("physical_map_source_sha256") == PHYSICAL_MAP_SOURCES[native["family"]], "native physical map source differs from the certified canonical map")
    require(native.get("mathematical_source") == "https://arxiv.org/abs/2112.07578v1", "native report uses another canonical normalization")
    boundary = native["boundary"]
    require(boundary["identity"] == native["canonical_identity"], "native boundary identity differs from its declared system")
    provenance = boundary["provenance"]
    require("Native automatic auxiliary-mass boundary;" in provenance and "no numerical reference seed" in provenance, "native AMF provenance is missing")
    require("extracted-map " in provenance, "native report lacks extracted physical-map provenance")
    require(boundary["leading_epsilon_power"] == 0 and boundary["last_epsilon_power"] == 4, "comparison requires the complete canonical Laurent range epsilon^0 through epsilon^4")
    require(boundary["verified_digits"] > 0, "native boundary has no verified precision")
    return boundary


def compare_files(native_path, reference_path=DEFAULT_REFERENCE, digits=20):
    native, native_hash = read_json(native_path)
    boundary = validate_native(native)
    require(isinstance(digits, int) and digits > 0, "comparison digits must be positive")

    # Never move reference access above the native-success gate.
    reference, reference_hash = read_json(reference_path)
    matches = [case for case in reference["cases"] if case["label"] == native["label"]]
    require(len(matches) == 1, "native label must identify exactly one physical SOURCE reference")
    case = matches[0]
    require(case["family"] == native["family"], "reference topology mismatch")
    systems = [system for system in reference["systems"] if system["family"] == case["family"]]
    require(len(systems) == 1, "ambiguous reference canonical system")
    system = systems[0]
    dimension = case["dimension"]
    require(dimension == {"Planar_EW1": 48, "NP_EW1": 61}[case["family"]], "reference does not contain the full physical basis")
    require(dimension == system["dimension"], "reference basis dimension mismatch")
    require(boundary["canonical_basis_indices"] == list(range(1, dimension + 1)), "native canonical basis order mismatch")
    coordinates = scoped_symbols(native["coordinates"], native["namespace"])
    expected_coordinates = dict(zip(system["coordinate_names"], case["start"], strict=True))
    require(coordinates.keys() == expected_coordinates.keys(), "native coordinate set mismatch")
    require(all(rational(coordinates[name]) == rational(value) for name, value in expected_coordinates.items()), "native point differs from the exact reference SOURCE, not destination")
    sheets = {"principal": "Principal", "opposite_principal": "Opposite"}
    expected_germs = {
        root["name"]: sheets[germ]
        for root, germ in zip(system["roots"], case["source_root_germs"], strict=True)
    }
    require(scoped_symbols(native["root_germs"], native["namespace"]) == expected_germs, "native root germs differ from the physical reference")
    cap = min(boundary["verified_digits"], case["source_verified_digits_cap"])
    require(digits <= cap, f"requested {digits} comparison digits exceed recorded native/reference cap {cap}")
    arrays = [boundary["coefficients"], boundary["absolute_errors"], case["source_values"], case["source_errors"]]
    require(all(len(array) == 5 and all(len(row) == dimension for row in array) for array in arrays), "incomplete Laurent coefficient/error array")

    maximum_difference = Fraction(0)
    maximum_scaled = Fraction(0)
    maximum_ratio = Fraction(0)
    zero_allowance_disagreement = False
    failures = []
    tolerance_squared = Fraction(1, 10 ** (2 * digits))
    for power in range(5):
        for index in range(dimension):
            actual, expected = arrays[0][power][index], arrays[2][power][index]
            native_error, reference_error = rational(arrays[1][power][index]), rational(arrays[3][power][index])
            require(native_error >= 0 and reference_error >= 0, "negative component error allowance")
            difference = difference_squared(actual, expected)
            allowance = native_error + reference_error
            scaled = difference / max(Fraction(1), norm_squared(expected))
            maximum_difference = max(maximum_difference, difference)
            maximum_scaled = max(maximum_scaled, scaled)
            if allowance:
                maximum_ratio = max(maximum_ratio, difference / allowance**2)
            elif difference:
                zero_allowance_disagreement = True
            if difference > allowance**2 or scaled > tolerance_squared:
                failures.append({
                    "epsilon_power": power, "canonical_basis_index": index + 1,
                    "native": actual, "reference": expected,
                    "absolute_difference": diagnostic_sqrt(difference),
                    "scaled_difference": diagnostic_sqrt(scaled),
                    "native_absolute_error": arrays[1][power][index],
                    "reference_absolute_error": arrays[3][power][index],
                    "within_combined_recorded_errors": difference <= allowance**2,
                    "within_requested_mixed_accuracy": scaled <= tolerance_squared,
                })
    return {
        "schema": "rustflow-native-gg-hg-source-comparison-v1",
        "status": "passed" if not failures else "reference_disagreement",
        "label": case["label"], "family": case["family"],
        "normalization": NORMALIZATION,
        "exact_source_coordinates": expected_coordinates,
        "root_germs": expected_germs,
        "epsilon_powers": [0, 1, 2, 3, 4], "checked_coefficients": 5 * dimension,
        "requested_mixed_comparison_digits": digits,
        "native_verified_digits": boundary["verified_digits"],
        "reference_verified_digits_cap": case["source_verified_digits_cap"],
        "corroborated_mixed_digits": digits if not failures else None,
        "maximum_absolute_difference": diagnostic_sqrt(maximum_difference),
        "maximum_scaled_difference": diagnostic_sqrt(maximum_scaled),
        "maximum_ratio_to_combined_recorded_errors": "infinity" if zero_allowance_disagreement else diagnostic_sqrt(maximum_ratio),
        "failure_count": len(failures), "failures": failures,
        "metric": "Both |native-reference| <= native_error+reference_error and |native-reference| <= 10^-digits*max(1,|reference|); exact rational squared-norm comparisons",
        "accuracy_note": "Independent recorded errors cap the comparison; printed mantissas do not certify additional digits. Diagnostics are rounded Decimal values, acceptance inequalities are exact. Coefficient comparison only, no finite-epsilon truncation or tail assertion.",
        "reference_opened_after_native_success": True,
        "native_report": {"path": str(native_path), "sha256": native_hash, "build": native["build"], "provenance": boundary["provenance"], "identity": boundary["identity"]},
        "reference": {"path": str(reference_path), "sha256": reference_hash, "upstream_provenance": reference["provenance"], "source_sha256": case["source_sha256"], "source_delta": case["source_delta"], "source_error_decimal_exponent": case["source_error_decimal_exponent"], "matrix_sha256": case["matrix_sha256"]},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_report", type=Path)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--digits", type=int, default=20)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    try:
        report = compare_files(arguments.native_report, arguments.reference, arguments.digits)
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.exit(2, f"comparison rejected: {error}\n")
    text = json.dumps(report, indent=2) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = arguments.output.with_suffix(arguments.output.suffix + ".tmp")
        temporary.write_text(text)
        temporary.replace(arguments.output)
    else:
        print(text, end="")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
