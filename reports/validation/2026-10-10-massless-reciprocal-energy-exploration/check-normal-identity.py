#!/usr/bin/env python3
"""Exact finite radial normal-IBP coefficient checks; no production predictions."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def binomial(value, order):
    result = Fraction(1)
    for k in range(order):
        result *= (value - k) / (k + 1)
    return result


def main():
    directory = Path(__file__).resolve().parent
    input_path = directory / "normal-identity-inputs.json"
    source = input_path.read_bytes()
    config = json.loads(source)
    records = []
    for dimension in config["spatial_dimensions"]:
        d = Fraction(dimension)
        for n in config["cut_orders"]:
            for p in config["inverse_energy_powers"]:
                a = d / 2 - 1
                beta = d - p - 1 - 2 * n
                assert beta != 0, "true-pole samples are not finite coefficient checks"
                raised = -n * binomial(a, n) / beta
                bulk = (p + 1) * binomial(a, n - 1) / (2 * beta)
                upper = binomial(a, n - 1) / 2
                assert raised + bulk + upper == 0
                records.append({
                    "d": str(d), "cut_order": n, "inverse_energy_power": p,
                    "beta": str(beta), "bare_origin_convergent": beta > 0,
                    "raised_coefficient": str(raised),
                    "bulk_coefficient": str(bulk),
                    "upper_surface_coefficient": str(upper),
                    "residual": "0",
                })
    result = {
        "status": "pass",
        "scope": "exact radial normal-source coefficient identity; no graph admission, native replay or numerical-flow claim",
        "input_sha256": hashlib.sha256(source).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": len(records),
        "records": records,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
