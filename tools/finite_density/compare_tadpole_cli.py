#!/usr/bin/env python3
"""Compare saved native CLI samples with independent D=3 tadpole identities."""
import hashlib
import json
import math
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly/native-interfaces"
getcontext().prec = 100


def atan_inverse(denominator):
    x = Decimal(1) / denominator
    x2 = x * x
    term = x
    result = term
    for k in range(1, 1000):
        term = -term * x2
        addition = term / (2 * k + 1)
        result += addition
        if abs(addition) < Decimal("1e-105"):
            return result
    raise ArithmeticError("arctangent series did not converge")


def euler_gamma():
    # H_n-ln(n)-1/(2n)+sum B_(2k)/(2k*n^(2k)); independently
    # constructed rational Bernoulli numbers, beyond the comparison precision.
    bernoulli = [Fraction(1)]
    for degree in range(1, 81):
        bernoulli.append(-sum((math.comb(degree + 1, k) * bernoulli[k] for k in range(degree)), Fraction(0)) / (degree + 1))
    n = Decimal(128)
    result = sum((Decimal(1) / k for k in range(1, 129)), Decimal(0)) - n.ln() - 1 / (2 * n)
    for degree in range(2, 81, 2):
        b = bernoulli[degree]
        result += Decimal(b.numerator) / Decimal(b.denominator) / (degree * n ** degree)
    return result


def main():
    # Machin's identity; no production normalization or special-function owner.
    pi = 16 * atan_inverse(5) - 4 * atan_inverse(239)
    references = {
        "above": [-1 / (4 * pi), Decimal(0)],
        "below": [-1 / (8 * pi), 1 / (4 * pi)],
    }
    comparisons = []
    hashes = {}
    for case, expected in references.items():
        path = REPORT / f"tadpole-{case}-resources.log"
        record = json.loads(path.read_text())
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        assert record["epsilon"] == "1/2" and record["dimension_base"] == 4
        assert record["normalization"] == "unscaled Euclidean amplitude"
        assert record["full_amplitude"] and record["verified_digits"] is None
        assert [c["cut_slots"] for c in record["contributions"]] == [[], [0]]
        if case == "above":
            assert record["occupied_reports"][0]["construction"] == "compact_polynomial_moments"
        else:
            assert record["empty_support"][0]["cut_slots"] == [0]
            assert all(Decimal(v["real"]) == 0 and Decimal(v["imaginary"]) == 0 for v in record["contributions"][1]["values"])
        for target, reference in enumerate(expected):
            value = record["values"][target]
            delta = Decimal(value["real"]) - reference
            imaginary = Decimal(value["imaginary"])
            absolute = (delta * delta + imaginary * imaginary).sqrt()
            relative = absolute / abs(reference) if reference else None
            small = abs(reference) < Decimal("1e-20")
            tolerance = Decimal("1e-25" if small else "1e-12")
            compared = absolute if small else relative
            comparisons.append({
                "case": case, "target_index": target, "prediction": value,
                "reference": str(reference), "real_difference": str(delta),
                "imaginary_difference": str(imaginary), "absolute_difference": str(absolute),
                "relative_difference": str(relative) if relative is not None else None,
                "criterion": "absolute" if small else "relative",
                "tolerance": str(tolerance), "passed": compared <= tolerance,
            })
    gamma = euler_gamma()
    a = Decimal(1) / 4
    mu = Decimal(1)
    radius = (mu * mu - a).sqrt()
    logarithm = ((mu + radius) / a.sqrt()).ln()
    laurent_reference = [
        {-1: -a / (16 * pi * pi),
         0: a / (16 * pi * pi) * (gamma - 1 - (4 * pi / a).ln())
         - (mu * radius - a * logarithm) / (8 * pi * pi)},
        {-1: 1 / (16 * pi * pi),
         0: ((4 * pi / a).ln() - gamma - 2 * logarithm) / (16 * pi * pi)},
    ]
    laurent_path = REPORT / "tadpole-laurent-resources.log"
    laurent = json.loads(laurent_path.read_text())
    hashes[str(laurent_path.relative_to(ROOT))] = hashlib.sha256(laurent_path.read_bytes()).hexdigest()
    assert laurent["operation"] == "finite-density" and laurent["full_amplitude"]
    assert laurent["normalization"] == "unscaled Euclidean amplitude"
    laurent_comparisons = []
    for target, coefficients in enumerate(laurent_reference):
        expansion = laurent["expansions"][target]
        assert expansion["verified_digits"] >= 12 and expansion["refinements"] >= 1
        for power, reference in coefficients.items():
            value = expansion["coefficients"][str(power)]
            delta = Decimal(value["real"]) - reference
            imaginary = Decimal(value["imaginary"])
            absolute = (delta * delta + imaginary * imaginary).sqrt()
            relative = absolute / abs(reference)
            laurent_comparisons.append({
                "target_index": target, "epsilon_power": power, "prediction": value,
                "reference": str(reference), "real_difference": str(delta),
                "imaginary_difference": str(imaginary), "absolute_difference": str(absolute),
                "relative_difference": str(relative), "criterion": "relative",
                "tolerance": "1e-12", "passed": relative <= Decimal("1e-12"),
            })
    passed = all(c["passed"] for c in comparisons)
    passed = passed and all(c["passed"] for c in laurent_comparisons)
    output = {
        "schema": 1, "status": "passed" if passed else "failed",
        "reference_method": "Independent Euclidean D=3 massive tadpole plus compact occupied shell: I1=-max(m,mu)/(4*pi), I2=1/(8*pi*m) below threshold and zero above; pi from Machin arctangent series at 100 Decimal digits",
        "reference_source": "analytic Gaussian vacuum integral and independent physical-mass derivative of occupied radial integral, not numerical oracle data",
        "oracle_numerical_records_read": 0,
        "scope": "four D=3 target/support-case values and four epsilon=-1,0 Laurent coefficients; complete one-loop CLI execution only, not a multiloop acceptance result",
        "laurent_reference_method": "Independent expansion of Gamma(epsilon-1)*a^(1-epsilon)/(4*pi)^(2-epsilon) minus the D=4 compact occupied radial integral; raised coefficients from -d/da at fixed original numerator and mu. Euler gamma from rational Bernoulli Euler-Maclaurin series, pi from Machin identity.",
        "source_snapshot": "../tangent-source-hashes.json",
        "inputs_sha256": hashes,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "comparisons": comparisons,
        "laurent_comparisons": laurent_comparisons,
    }
    (REPORT / "tadpole-independent-comparison.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"], "sample_comparisons": len(comparisons), "laurent_comparisons": len(laurent_comparisons), "above_raised_absolute_error": comparisons[1]["absolute_difference"], "largest_laurent_relative_difference": str(max(Decimal(c["relative_difference"]) for c in laurent_comparisons))}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
