#!/usr/bin/env python3
"""Exact, validation-only factorized compact Ward checks; no native/oracle input."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUT = ROOT / "examples/finite_density/massless_three_loop_chain.json"


def rank(rows):
    a = [list(row) for row in rows]
    r = 0
    for c in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [v / p for v in a[r]]
        for i in range(len(a)):
            if i != r:
                p = a[i][c]
                a[i] = [u - p * v for u, v in zip(a[i], a[r])]
        r += 1
    return r


def inverse(matrix):
    n = len(matrix)
    a = [list(row) + [Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for c in range(n):
        pivot = next(i for i in range(c, n) if a[i][c])
        a[c], a[pivot] = a[pivot], a[c]
        p = a[c][c]
        a[c] = [v / p for v in a[c]]
        for i in range(n):
            if i != c:
                p = a[i][c]
                a[i] = [u - p * v for u, v in zip(a[i], a[c])]
    return [row[n:] for row in a]


def binomial(a, n):
    result = Q(1)
    for k in range(n):
        result *= (a - k) / (k + 1)
    return result


def moment(dimension, n, upper, energy_power, mu):
    """Raw C_n H_s moment, with the common positive angular factor omitted.

    Fixed-independent-energy mass derivatives give (-1)^(n-1) binom((D-3)/2,n-1).
    Positive H_s is (-1)^(s-1)/(s-1)! times the s-th mu derivative of H_0.
    All chosen samples lie in the convergent origin domain before continuation.
    """
    alpha = Q(dimension - 3, 2)
    beta = dimension + energy_power - 2 * n
    assert beta > 0
    coefficient = (-1) ** (n - 1) * binomial(alpha, n - 1)
    if upper == 0:
        return coefficient * mu ** beta / beta
    derivative = Q(1)
    for k in range(upper - 1):
        derivative *= Q(beta - 1 - k, k + 1)
    return coefficient * (-1) ** (upper - 1) * derivative * mu ** (beta - upper)


def main():
    output = HERE / "exact-checks.json"
    if output.exists():
        raise SystemExit("refusing to overwrite completed evidence")
    samples = []
    for d in [7, 9]:
        for s in range(4):
            for other_n in [1, 2, 3]:
                for other_s in range(3):
                    for other_energy in range(3):
                        mu, other_mu = Q(2, 3), Q(5, 4)
                        rest = moment(d, other_n, other_s, other_energy, other_mu)
                        lhs = (d - 2) * moment(d, 1, s, 0, mu)
                        lhs += (-1 if s == 0 else s) * moment(d, 1, s + 1, 1, mu)
                        assert lhs * rest == 0
                        samples.append({"D": d, "own_upper": s, "other_cut": other_n,
                            "other_upper": other_s, "other_energy_power": other_energy,
                            "other_moment": str(rest), "residual": "0"})
    assert any(row["other_cut"] > 1 and row["other_upper"] > 0 and row["other_moment"] != "0" for row in samples)
    forbidden = {}
    for name, n, r in [("own_raised_C2", 2, 0), ("own_energy_numerator", 1, 1)]:
        value = 5 * moment(7, n, 0, r, Q(2, 3)) - moment(7, n, 1, r + 1, Q(2, 3))
        assert value != 0
        forbidden[name] = str(value)
    definition = json.loads(INPUT.read_text())
    cuts = [0, 3]
    rows = [[Q(v) for v in edge["routing"]] for edge in definition["edges"]]
    routing = [rows[slot] for slot in cuts]
    for i in range(definition["loops"]):
        candidate = [Q(i == j) for j in range(definition["loops"])]
        if rank(routing + [candidate]) > len(routing):
            routing.append(candidate)
    inv = inverse(routing)
    transformed = [[sum((row[i] * inv[i][j] for i in range(len(inv))), Q(0))
                    for j in range(len(inv))] for row in rows]
    guards = []
    for compact, slot in enumerate(cuts):
        dependent = [i for i, row in enumerate(transformed) if i not in cuts and row[compact]]
        guards.append({"compact_loop": compact, "cut_slot": slot,
            "ordinary_physical_base_indices_required_zero": dependent,
            "completion_guard": "also every actual completion with incident Gram/E dependence; not inferred here"})
    output.write_text(json.dumps({
        "scope": "Exact source-identity assessment only; no native closure or numerical prediction",
        "production_admission_changed": False, "native_or_oracle_values_read": 0,
        "checks": len(samples), "samples": samples, "forbidden_case_residuals": forbidden,
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cuts": cuts, "routing": [[str(v) for v in row] for row in routing],
        "transformed_physical_routings": [[str(v) for v in row] for row in transformed],
        "necessary_physical_guards": guards,
    }, indent=2) + "\n")
    print(json.dumps({"checks": len(samples), "status": "passed", "guards": guards}))


if __name__ == "__main__":
    main()
