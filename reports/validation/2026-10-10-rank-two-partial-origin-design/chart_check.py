#!/usr/bin/env python3
"""Exact algebra-only replay of complete rank<=2 chart covers (no native admission)."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
OPS = 0
OP_LIMIT = 2_000_000
MAX_TERMS = 20_000


def budget(n=1):
    global OPS
    OPS += n
    if OPS > OP_LIMIT:
        raise RuntimeError('exact arithmetic budget exceeded; no partial report')


class Poly:
    def __init__(self, n, terms=None):
        self.n = n
        self.t = {tuple(k): Q(v) for k, v in (terms or {}).items() if v}
        assert all(len(k) == n and all(e >= 0 for e in k) for k in self.t)
        if len(self.t) > MAX_TERMS:
            raise RuntimeError('polynomial term budget exceeded; no partial report')

    @classmethod
    def one(cls, n):
        return cls(n, {(0,) * n: 1})

    @classmethod
    def var(cls, n, j):
        e = [0] * n
        e[j] = 1
        return cls(n, {tuple(e): 1})

    def scale(self, c):
        budget(len(self.t))
        return Poly(self.n, {m: c * v for m, v in self.t.items()})

    def __add__(self, b):
        assert self.n == b.n
        out = dict(self.t)
        for m, c in b.t.items():
            budget()
            out[m] = out.get(m, Q(0)) + c
        return Poly(self.n, out)

    def __sub__(self, b):
        return self + b.scale(-1)

    def __mul__(self, b):
        assert self.n == b.n
        out = {}
        for a, av in self.t.items():
            for c, cv in b.t.items():
                budget()
                m = tuple(x + y for x, y in zip(a, c))
                out[m] = out.get(m, Q(0)) + av * cv
        return Poly(self.n, out)

    def __eq__(self, b):
        return self.n == b.n and self.t == b.t

    def divide_x(self):
        assert all(m[0] >= 1 for m in self.t), 'nondivisible chart polynomial'
        return Poly(self.n, {(m[0] - 1,) + m[1:]: c for m, c in self.t.items()})

    def record(self):
        return [[list(m), str(c)] for m, c in sorted(self.t.items())]


def sum_poly(n, terms):
    out = Poly(n)
    for term in terms:
        out = out + term
    return out


def det2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def rank(rows):
    nonzero = [a for a in rows if any(a)]
    if not nonzero:
        return 0
    if len(nonzero[0]) == 1:
        return 1
    return 2 if any(det2(a, b) for a, b in combinations(nonzero, 2)) else 1


def chart(rows, ext, primary, pivot):
    s = len(rows)
    n = s - 1
    alpha = [None] * s
    alpha[primary] = Poly.one(n)
    if len(rows[0]) == 1:
        j = 0
        for e in range(s):
            if e != primary:
                alpha[e] = Poly.var(n, j)
                j += 1
        A = sum_poly(n, [alpha[e].scale(rows[e][0] ** 2) for e in range(s)])
        assert all(c >= 0 for c in A.t.values())
        assert A.t.get((0,) * n) == rows[primary][0] ** 2 > 0
        return {'primary': primary, 'pivot': None, 'unit_lower_bound': str(rows[primary][0] ** 2),
                'U_terms': len(A.t), 'identity_sha256': digest(A.record())}
    assert pivot is not None and det2(rows[primary], rows[pivot])
    x = Poly.var(n, 0)
    alpha[pivot] = x
    j = 1
    transverse = 0
    for e in range(s):
        if det2(rows[primary], rows[e]):
            transverse += 1
        if e in (primary, pivot):
            continue
        alpha[e] = Poly.var(n, j)
        if det2(rows[primary], rows[e]):
            alpha[e] = x * alpha[e]
        j += 1
    assert j == n
    A00 = sum_poly(n, [alpha[e].scale(rows[e][0] ** 2) for e in range(s)])
    A01 = sum_poly(n, [alpha[e].scale(rows[e][0] * rows[e][1]) for e in range(s)])
    A11 = sum_poly(n, [alpha[e].scale(rows[e][1] ** 2) for e in range(s)])
    U = A00 * A11 - A01 * A01
    cb = sum_poly(n, [(alpha[e] * alpha[f]).scale(det2(rows[e], rows[f]) ** 2)
                      for e, f in combinations(range(s), 2)])
    assert U == cb, 'Cauchy--Binet U identity'
    U0 = U.divide_x()
    lower = det2(rows[primary], rows[pivot]) ** 2
    assert all(c >= 0 for c in U0.t.values())
    assert U0.t.get((0,) * n) == lower > 0
    k = len(ext[0])
    B = [[sum_poly(n, [alpha[e].scale(rows[e][r] * ext[e][c]) for e in range(s)])
          for c in range(k)] for r in range(2)]
    M = [[A11 * B[0][c] - A01 * B[1][c] for c in range(k)],
         [A00 * B[1][c] - A01 * B[0][c] for c in range(k)]]
    M0 = [[p.divide_x() for p in row] for row in M]
    V0 = []
    for a in range(k):
        vr = []
        for b in range(k):
            C = sum_poly(n, [alpha[e].scale(ext[e][a] * ext[e][b]) for e in range(s)])
            V = U * C - B[0][a] * M[0][b] - B[1][a] * M[1][b]
            gram = Poly(n)
            for es in combinations(range(s), 3):
                da = det3(*[rows[e] + [ext[e][a]] for e in es])
                db = det3(*[rows[e] + [ext[e][b]] for e in es])
                gram = gram + (alpha[es[0]] * alpha[es[1]] * alpha[es[2]]).scale(da * db)
            assert V == gram, 'full external Gram coefficient Cauchy--Binet identity'
            vr.append(V.divide_x())
        V0.append(vr)
    # Polynomial cofactor numerators and U=x*unit imply at most x^-1 per covariance.
    payload = {'U0': U0.record(), 'M0': [[p.record() for p in row] for row in M0],
               'V0': [[p.record() for p in row] for row in V0]}
    return {'primary': primary, 'pivot': pivot, 'transverse_rows': transverse,
            'jacobian_x_power': transverse - 1, 'unit_lower_bound': str(lower),
            'U_terms': len(U.t), 'mean_terms': sum(len(p.t) for row in M for p in row),
            'V_terms': sum(len(p.t) for row in V0 for p in row),
            'identity_sha256': digest(payload)}


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def Qrows(values):
    return [[Q(x) for x in row] for row in values]


CASES = [
    ('parallel_pair', [[1, 0], [1, 0], [0, 1], [-1, 1]], [[1, 0], [0, 0], [0, 0], [0, 0]]),
    ('rational_shear', [['1/2', '3/2'], [1, 3], [2, -1], ['-3/2', '5/2']],
     [[1, -1], [0, 0], ['2/3', 1], [-1, 2]]),
    ('three_parallel_classes', [[1, 0], [-2, 0], [0, 1], [0, -3], [1, 1], [-2, -2]],
     [[1, 0], [0, 0], [0, 1], [0, 0], [1, -1], ['1/2', '3/4']]),
    ('generic_rank_two', [[1, 2], [2, -1], [3, 1], [-1, 3], [2, 3]],
     [[1, 0], [0, 1], [1, 1], [-1, 2], [0, 0]]),
    ('rank_one_nonunit', [[1], [-2], ['3/2'], ['-1/3']],
     [[1, 0], [0, 1], [1, -1], ['2/3', '1/2']]),
]


def main():
    records = []
    support_count = 0
    full_count = 0
    chart_count = 0
    for name, rv, cv in CASES:
        rows, ext = Qrows(rv), Qrows(cv)
        h = len(rows[0])
        support_records = []
        for mask in range(1 << len(rows)):
            support_count += 1
            active = [e for e in range(len(rows)) if mask & (1 << e)]
            r = [rows[e] for e in active]
            c = [ext[e] for e in active]
            actual_rank = rank(r)
            rec = {'active': active, 'rank': actual_rank}
            if actual_rank == h:
                full_count += 1
                choices = [(i, None) for i in range(len(r))] if h == 1 else [
                    (i, j) for i in range(len(r)) for j in range(len(r)) if det2(r[i], r[j])]
                rec['charts'] = [chart(r, c, i, j) for i, j in choices]
                chart_count += len(choices)
                rec['cover'] = 'all maximum-primary / maximum-nonparallel ordered pairs' if h == 2 else 'all maximum-primary gauges'
            support_records.append(rec)
        records.append({'name': name, 'rows': [[str(x) for x in a] for a in rows],
                        'external_rows': [[str(x) for x in a] for a in ext],
                        'supports': support_records})
    # Every finite integer J gets a proof neighborhood; no universal fixed D is claimed.
    jet_checks = []
    for P in (1, 2, 4, 9, 17):
        for R in (0, 1, 4, 11):
            for J in (0, 1, 3, 19):
                D = 2 * (P + R + J + 2) + Q(1, 3)
                N = (D / 2 + R).__floor__() + 1
                margin = D - P - J - N
                assert D / 2 + R < N < D / 2 + R + 2
                assert margin > 0
                jet_checks.append({'P': P, 'R': R, 'J_including_lower_contacts': J,
                                   'D': str(D), 'fixed_N': N, 'lambda_minus_J_minus_N': str(margin)})
    inventory = Qrows([[1, 0], [1, 0], [0, 1], [-1, 1]])
    active, shifted = [1, 2, 3], [0, 2]
    actual_rank = rank([inventory[e] for e in active])
    shifted_rank = rank([inventory[e] for e in active if e in shifted])
    assert actual_rank == 2 and shifted_rank == 1
    negative = {'complete_rows': ['t+q', 't', 'l', 'l-t'], 'shifted_positions': shifted,
                'active_positions': active, 'actual_virtual_rank': actual_rank,
                'incorrect_shifted_only_rank': shifted_rank,
                'correct_class': 'full-rank shifted-positive virtual vacuum (V=0), not finite-eta zero',
                'no_integral_value_or_native_zero_claim': True}
    out = {'scope': 'exact rational chart algebra and finite-jet inequalities only; no native/HEPKit replay, physical-contour certificate, production permit, or numerical result',
           'complete_supports': support_count, 'full_rank_supports': full_count,
           'complete_charts': chart_count, 'finite_jet_checks': len(jet_checks),
           'arithmetic_operations': OPS, 'arithmetic_budget': OP_LIMIT,
           'max_terms_per_polynomial': MAX_TERMS, 'cases': records,
           'finite_jet_witnesses': jet_checks, 'shifted_only_rank_negative': negative,
           'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / 'chart-checks.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: out[k] for k in ('complete_supports', 'full_rank_supports', 'complete_charts', 'finite_jet_checks', 'arithmetic_operations')}))


if __name__ == '__main__':
    main()
