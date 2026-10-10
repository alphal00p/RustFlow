#!/usr/bin/env python3
"""Independent convergent Schwinger/simplex quadrature; no native predictions.

Only Python's Decimal/Fraction arithmetic is used. Gauss-Legendre nodes are
refined at the requested precision. Gamma uses shifted Stirling with exact
Bernoulli numbers and the positive-real next-term remainder bound.
"""
import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import time


def bernoulli(n):
    a = []
    for m in range(n + 1):
        a.append(Q(1, m + 1))
        for j in range(m, 0, -1):
            a[j - 1] = j * (a[j - 1] - a[j])
    return a[0]


def atan_inverse(q, tolerance):
    x = D(1) / q
    power = x
    total = power
    for k in range(1, 10000):
        power *= -x * x
        term = power / (2 * k + 1)
        total += term
        if abs(term) < tolerance:
            return total
    raise ValueError('pi series exhausted')


def gamma_positive(z, digits):
    tolerance = D(10) ** (-digits - 10)
    pi = 16 * atan_inverse(D(5), tolerance) - 4 * atan_inverse(D(239), tolerance)
    w = z + 128
    value = (w - D('0.5')) * w.ln() - w + (2 * pi).ln() / 2
    previous = D('Infinity')
    for k in range(1, 100):
        b = bernoulli(2 * k)
        term = (D(b.numerator) / b.denominator) / (2 * k * (2 * k - 1) * w ** (2 * k - 1))
        if abs(term) >= previous:
            raise ValueError('Stirling remainder stopped decreasing')
        if abs(term) < tolerance:
            remainder = abs(term)
            break
        value += term
        previous = abs(term)
    else:
        raise ValueError('Gamma work limit exhausted')
    for j in range(128):
        value -= (z + j).ln()
    return value.exp(), remainder, k


def gauss_legendre(n, digits):
    tolerance = D(10) ** (-digits - 7)
    nodes = []
    for i in range(1, n // 2 + 1):
        x = D(str(math.cos(math.pi * (i - .25) / (n + .5))))
        for iteration in range(40):
            previous, current = D(1), x
            for k in range(2, n + 1):
                previous, current = current, ((2 * k - 1) * x * current - (k - 1) * previous) / k
            derivative = n * (x * current - previous) / (x * x - 1)
            delta = current / derivative
            x -= delta
            if abs(delta) < tolerance:
                break
        else:
            raise ValueError('Legendre Newton work limit exhausted')
        # Refresh the derivative at the final root before computing its weight.
        previous, current = D(1), x
        for k in range(2, n + 1):
            previous, current = current, ((2 * k - 1) * x * current - (k - 1) * previous) / k
        derivative = n * (x * current - previous) / (x * x - 1)
        weight = 1 / ((1 - x * x) * derivative * derivative)
        nodes.extend([((1 - x) / 2, weight), ((1 + x) / 2, weight)])
    assert len(nodes) == n
    assert abs(sum(w for _, w in nodes) - 1) < D(10) ** (-digits)
    return sorted(nodes)


def evaluate(digits, order, map_power):
    started = time.monotonic()
    with localcontext() as context:
        context.prec = digits + 18
        nodes = gauss_legendre(order, digits)
        dimension = D(12) / 5
        # Primary max-alpha sectors, ordered x/y chart, then t=s^map_power.
        # Integral = 6p integral s^(p*(2-D/2)-1)
        # (1+s^p(1+v))^(D-3) (1+v+s^p*v)^(-D/2) ds dv.
        exponent = map_power * (2 - dimension / 2) - 1
        assert exponent == int(exponent) and exponent >= 0
        total = D(0)
        for s, ws in nodes:
            t = s ** map_power
            prefactor = ws * s ** int(exponent)
            inner = D(0)
            for v, wv in nodes:
                first = ((1 + t * (1 + v)).ln() * (dimension - 3)).exp()
                second = ((1 + v + t * v).ln() * (-dimension / 2)).exp()
                inner += wv * first * second
            total += prefactor * inner
        simplex = 6 * map_power * total
        gamma, gamma_log_error, gamma_terms = gamma_positive(3 - dimension, digits)
        euclidean = gamma * simplex
        return {'decimal_digits': digits, 'nodes_per_axis': order, 'map_power': map_power,
                'function_evaluations': order * order, 'simplex_integral': str(simplex),
                'gamma_factor': str(gamma), 'gamma_log_truncation_bound': str(gamma_log_error),
                'gamma_stirling_terms_checked': gamma_terms,
                'euclidean_pi_normalized': str(euclidean), 'native_minkowski': str(-euclidean),
                'wall_seconds': time.monotonic() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('preserve existing reference outputs')
    profiles = [(35, 32, 5), (50, 64, 5), (65, 96, 5), (65, 96, 10)]
    rows = []
    for profile in profiles:
        row = evaluate(*profile)
        rows.append(row)
        print(json.dumps({'profile': profile, 'wall_seconds': row['wall_seconds']}), flush=True)
    with localcontext() as context:
        context.prec = 80
        selected = D(rows[-2]['native_minkowski'])
        changes = [abs(D(row['native_minkowski']) - selected) / abs(selected) for row in rows]
        assert max(changes[1:]) < D('1e-20'), 'quadrature/order/map refinement insufficient'
        assert changes[0] < D('1e-12'), 'coarse quadrature sanity check failed'
        report = {'scope': 'Independent convergent reference only; no native predictions or oracle records read.',
                  'dimension': '12/5', 'epsilon': '4/5', 'masses_squared': ['1', '1', '1'],
                  'powers': [1, 1, 1], 'normalization': 'Two loop measures d^Dk/(i*pi^(D/2)); native sign (-1)^3 relative to positive Euclidean pi-normalized integral.',
                  'prescription': 'positive unit masses; convergent real Schwinger integral',
                  'profiles': rows, 'relative_changes_from_selected': [str(x) for x in changes],
                  'reference_native': str(selected), 'reference_absolute_imaginary_part': '0',
                  'accepted_comparison_relative_tolerance': '1e-10',
                  'precision_claim': '10 decimal digits, empirically verified by order, arithmetic precision and coordinate-map refinement. Gamma truncation separately bounded; quadrature error is not an interval enclosure.',
                  'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'oracle_records_read': 0, 'native_prediction_files_read': 0}
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
