#!/usr/bin/env python3
"""Independent exact Cartesian polynomial and radial boost checks."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            out[monomial] = out.get(monomial, Q(0)) + coefficient
    return {m: c for m, c in out.items() if c}


def scale(polynomial, coefficient):
    return {m: c * coefficient for m, c in polynomial.items() if c * coefficient}


def mul(a, b):
    out = {}
    for x, c in a.items():
        for y, d in b.items():
            z = tuple(i + j for i, j in zip(x, y))
            out[z] = out.get(z, Q(0)) + c * d
    return {m: c for m, c in out.items() if c}


def derivative(polynomial, axis):
    out = {}
    for m, c in polynomial.items():
        if m[axis]:
            power = list(m)
            power[axis] -= 1
            out[tuple(power)] = c * m[axis]
    return out


def binomial(a, n):
    value = Q(1)
    for j in range(n):
        value *= (a - j) / (j + 1)
    return value


def main():
    cartesian = []
    gram_checks = energy_checks = 0
    for dimension in range(2, 6):
        for loops in range(1, 5):
            variables = loops * dimension
            vectors = []
            for i in range(loops):
                vector = []
                for mu in range(dimension):
                    powers = [0] * variables
                    powers[i * dimension + mu] = 1
                    vector.append({tuple(powers): Q(1)})
                vectors.append(vector)

            def dot(a, b):
                return add(*(scale(mul(a[mu], b[mu]), 1 if mu == 0 else -1)
                             for mu in range(dimension)))

            routings = [[Q(int(i == a)) for i in range(loops)] for a in range(loops)]
            routings.append([Q(i + 1, 2) * (-1) ** i for i in range(loops)])
            for routing in routings:
                routed = [add(*(scale(vectors[i][mu], routing[i]) for i in range(loops)))
                          for mu in range(dimension)]
                fields = []
                for i in range(loops):
                    contraction = dot(routed, vectors[i])
                    fields.append([add(mul(vectors[i][0], routed[mu]),
                                       scale(contraction, -1) if mu == 0 else {})
                                   for mu in range(dimension)])

                def action(polynomial):
                    return add(*(mul(fields[i][mu], derivative(polynomial, i * dimension + mu))
                                 for i in range(loops) for mu in range(dimension)))

                divergence = add(*(derivative(fields[i][mu], i * dimension + mu)
                                   for i in range(loops) for mu in range(dimension)))
                assert divergence == scale(routed[0], dimension - 1)
                for i in range(loops):
                    for j in range(i, loops):
                        assert not action(dot(vectors[i], vectors[j]))
                        gram_checks += 1
                    expected = add(mul(vectors[i][0], routed[0]), scale(dot(routed, vectors[i]), -1))
                    assert action(vectors[i][0]) == expected
                    energy_checks += 1
                cartesian.append({"D": dimension, "loops": loops,
                                  "routing": [str(r) for r in routing],
                                  "all_gram_actions_zero": True, "divergence_exact": True})

    radial = []
    for d in [Q(19, 3), Q(37, 3), Q(61, 3)]:
        def bulk(n, power):
            if n <= 0:
                return Q(0)
            return (-1) ** (n - 1) * binomial(d / 2 - 1, n - 1) / (d + power + 1 - 2 * n)

        def surface(n, power):
            return Q(0) if n <= 0 else (-1) ** (n - 1) * binomial(d / 2 - 1, n - 1)

        for n in range(1, 6):
            for r in range(6):
                terms = [(d + r) * bulk(n, r + 1), -r * bulk(n - 1, r - 1),
                         -surface(n, r + 2), surface(n - 1, r)]
                assert sum(terms) == 0
                radial.append({"d": str(d), "cut_order": n, "energy_power": r,
                               "coefficients": [str(t) for t in terms], "residual": "0"})
    positive_upper = []
    for d in [7, 9, 11]:
        for mu in [Q(1), Q(3, 2)]:
            def upper_moment(n, s, power):
                if n <= 0:
                    return Q(0)
                degree = d + power - 2 * n
                # falling(degree,s-1)/(s-1)! = binom(degree,s-1).
                return (Q((-1) ** (n + s - 2)) * binomial(Q(d, 2) - 1, n - 1)
                        * binomial(Q(degree), s - 1) * mu ** (degree - s + 1))

            for n in range(1, 6):
                for r in range(6):
                    for s in range(1, 5):
                        terms = [(d + r) * upper_moment(n, s, r + 1),
                                 -r * upper_moment(n - 1, s, r - 1),
                                 s * upper_moment(n, s + 1, r + 2),
                                 -s * upper_moment(n - 1, s + 1, r)]
                        assert sum(terms) == 0
                        positive_upper.append({"d": d, "mu": str(mu), "cut_order": n,
                                               "energy_power": r, "upper_index": s,
                                               "coefficients": [str(t) for t in terms],
                                               "residual": "0"})
    print(json.dumps({"status": "pass", "scope": "exact Cartesian and radial identities only; no native closure or physical amplitude run",
                      "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "cartesian_fields": len(cartesian), "gram_actions": gram_checks,
                      "energy_actions": energy_checks, "radial_cases": len(radial),
                      "positive_upper_cases": len(positive_upper),
                      "cartesian": cartesian, "radial": radial,
                      "positive_upper": positive_upper}, indent=2))


if __name__ == "__main__":
    main()
