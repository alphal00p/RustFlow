#!/usr/bin/env python3
"""Exact radial distribution checks, not a native closure/evaluation test."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json


def falling(x, n):
    value = Q(1)
    for i in range(n):
        value *= x - i
    return value


def moment(D, mu, n=1, r=0, s=0):
    # Raw native Cn measure, with common sphere factor A_(D-1)/2 removed.
    alpha = Q(D - 3, 2)
    prefactor = (-1) ** (n - 1) * falling(alpha, n - 1) / factorial(n - 1)
    beta = D + r - 2 * n
    if s == 0:
        assert beta != 0
        return prefactor * mu**beta / beta
    return prefactor * (-1) ** (s - 1) * falling(beta - 1, s - 1) / factorial(s - 1) * mu**(beta - s)


checks = []
for D in (7, 9, 11):
    for mu in (Q(1), Q(3, 2)):
        for H in range(1, 6):
            value = ((D - 2) * moment(D, mu) - mu * moment(D, mu, s=1)) if H == 1 else (
                (D - 1 - H) * moment(D, mu, s=H - 1)
                + (H - 1) * mu * moment(D, mu, s=H)
            )
            assert value == 0
            checks.append({"D": D, "mu": str(mu), "upper_index": H, "residual": str(value)})

negative = []
for D in (7, 9, 11):
    for mu in (Q(1), Q(3, 2)):
        for reason, n, r in (("raised cut excluded", 2, 0), ("energy numerator excluded", 1, 2)):
            residual = (D - 2) * moment(D, mu, n=n, r=r) - mu * moment(D, mu, n=n, r=r, s=1)
            assert residual != 0
            negative.append({"reason": reason, "D": D, "mu": str(mu), "n": n, "r": r, "residual": str(residual)})
for mu in (Q(1), Q(3, 2)):
    m = Q(1, 2)
    # At D3: raw C1 bulk=(mu-m), upper=1 after removing A2/2.
    residual = (mu - m) - mu
    assert residual != 0
    negative.append({"reason": "positive cut mass excluded", "D": 3, "mu": str(mu), "mass": str(m), "residual": str(residual)})

here = Path(__file__).resolve().parent
result = {
    "scope": "Exact radial distribution checks of the proposed polynomial source guards; no production/native flow admission",
    "normalization": "raw native Cn; common A_(D-1)/2 divided out",
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "proof_sha256": hashlib.sha256((here / 'normalized-polynomial-c1-proof.md').read_bytes()).hexdigest(),
    "valid_checks": checks,
    "excluded_domain_counterexamples": negative,
    "passed": len(checks),
    "counterexamples_verified": len(negative),
}
(here / "normalized-polynomial-c1-checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"passed": len(checks), "counterexamples_verified": len(negative)}))
