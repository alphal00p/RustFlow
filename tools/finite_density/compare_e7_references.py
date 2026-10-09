#!/usr/bin/env python3
"""Compare a saved independent E7 reference with supplied I37, never AMF data."""
import ast
import hashlib
import json
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

from compare_tadpole_cli import atan_inverse

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly/independent-e7-reference"
ORACLE = ROOT / "fixtures/finite_density/oracle_results_supplied.json"
REFERENCE = DIRECTORY / "reference.json"
getcontext().prec = 100


def polynomial(expression, constants):
    """Restricted exact Q[pi] evaluator; reject all functions and other symbols."""
    def add(left, right):
        result = left.copy()
        for power, value in right.items():
            result[power] = result.get(power, Fraction(0)) + value
        return {power: value for power, value in result.items() if value}

    def multiply(left, right):
        result = {}
        for lp, lv in left.items():
            for rp, rv in right.items():
                result[lp + rp] = result.get(lp + rp, Fraction(0)) + lv * rv
        return {power: value for power, value in result.items() if value}

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return {0: Fraction(node.value)} if node.value else {}
        if isinstance(node, ast.Name):
            if node.id == "pi":
                return {1: Fraction(1)}
            if node.id in constants:
                return polynomial(constants[node.id], {})
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            sign = -1 if isinstance(node.op, ast.USub) else 1
            return {power: sign * value for power, value in visit(node.operand).items()}
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return add(left, right)
            if isinstance(node.op, ast.Sub):
                return add(left, {power: -value for power, value in right.items()})
            if isinstance(node.op, ast.Mult):
                return multiply(left, right)
            if isinstance(node.op, ast.Div) and set(right) == {0}:
                return {power: value / right[0] for power, value in left.items()}
            if isinstance(node.op, ast.Pow) and set(right) == {0}:
                exponent = right[0]
                assert exponent.denominator == 1 and 0 <= exponent <= 10
                result = {0: Fraction(1)}
                for _ in range(int(exponent)):
                    result = multiply(result, left)
                return result
        raise ValueError(f"unsupported exact reference expression: {ast.dump(node)}")

    return visit(ast.parse(expression.replace("^", "**"), mode="eval").body)


def main():
    helper = ROOT / "tools/finite_density/compare_tadpole_cli.py"
    source_bytes = {str(path.relative_to(ROOT)): path.read_bytes() for path in (ORACLE, REFERENCE, helper)}
    oracle, reference = json.loads(ORACLE.read_text()), json.loads(REFERENCE.read_text())
    assert oracle["schema_version"] == 3
    assert reference["status"] == "reference_checks_passed"
    record = next(record for record in oracle["four_loop_oracles"] if record["id"] == "I37")
    definition = reference["definition"]
    target = record["finite_density_target"]
    assert definition["loops"] == 4 and definition["chemical_potentials"] == ["1"]
    assert definition["loop_charges"] == [[value] for value in record["loop_signatures"]]
    assert definition["numerator_convention"] == "shifted_euclidean"
    assert definition["targets"][0]["numerator"] == target["numerator_polynomial"]
    assert definition["targets"][0]["powers"] == target["denominator_powers"][:7]
    assert all(power == 0 for power in target["denominator_powers"][7:])
    routes = oracle["bases"][record["basis"]]["positive_index_routings"]
    assert [[int(value) for value in edge["routing"]] for edge in definition["edges"]] == routes[:7]
    assert all(Fraction(edge["mass_squared"]) == 0 for edge in definition["edges"])
    assert oracle["conventions"]["normalization"] == "(4*pi)^(-2*L)*(Lambda_bar/2)^(2*L*eps)"
    assert oracle["conventions"]["measure_factor"] == "(exp(symbolica::γ)*Lambda_bar^2/(4*pi))^eps"
    assert not record.get("uncertainties")

    pi = 16 * atan_inverse(5) - 4 * atan_inverse(239)
    comparisons = []
    for power, supplied in record["coefficients"].items():
        derived = reference["analytic_coefficients"][0][power]
        supplied_exact = polynomial(supplied, {"a1": oracle["constants"]["a1"]})
        derived_exact = polynomial(derived, {})
        supplied_decimal = sum((Decimal(value.numerator) / value.denominator * pi ** degree for degree, value in supplied_exact.items()), Decimal(0))
        native_string = reference["native_100_digit_coefficients"][0][power]
        assert native_string.startswith("(") and native_string.endswith("+0i)")
        native_decimal = Decimal(native_string[1:-4])
        delta = native_decimal - supplied_decimal
        relative = abs(delta) / abs(supplied_decimal)
        comparisons.append({
            "epsilon_power": int(power),
            "supplied_expression": supplied,
            "independent_expression": derived,
            "supplied_q_pi_polynomial": {str(k): str(v) for k, v in sorted(supplied_exact.items())},
            "independent_q_pi_polynomial": {str(k): str(v) for k, v in sorted(derived_exact.items())},
            "exact_equality": supplied_exact == derived_exact,
            "saved_independent_native_decimal": native_string,
            "supplied_decimal": str(supplied_decimal),
            "decimal_difference": str(delta),
            "decimal_relative_difference": str(relative),
            "decimal_tolerance": "1e-90",
            "passed": supplied_exact == derived_exact and relative <= Decimal("1e-90"),
        })
    for path in (ORACLE, REFERENCE):
        assert path.read_bytes() == source_bytes[str(path.relative_to(ROOT))]
    passed = all(comparison["passed"] for comparison in comparisons)
    output = {
        "schema": 1,
        "status": "passed" if passed else "failed",
        "kind": "reference_vs_reference_only",
        "supplied_record": "I37",
        "independent_target_index": 0,
        "supplied_reference_records_compared": 1,
        "supplied_coefficients_compared": len(comparisons),
        "amf_predictions_read": 0,
        "amf_oracle_comparisons": 0,
        "native_amf_acceptance": False,
        "method": "Exact rational polynomial identity in Q[pi] after resolving supplied a1; saved independent native 100-digit coefficient values also compared to 100-digit Decimal Machin-pi evaluation. No predictions or generation code are changed.",
        "normalization": reference["normalization"],
        "resolved_supplied_constants": {"a1": oracle["constants"]["a1"]},
        "inputs_sha256": {name: hashlib.sha256(data).hexdigest() for name, data in source_bytes.items()},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "comparisons": comparisons,
    }
    (DIRECTORY / "supplied-I37-reference-comparison.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"], "kind": output["kind"], "exact_coefficient_identities": len(comparisons), "largest_decimal_relative_difference": str(max(Decimal(item["decimal_relative_difference"]) for item in comparisons)), "amf_predictions_read": 0}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
