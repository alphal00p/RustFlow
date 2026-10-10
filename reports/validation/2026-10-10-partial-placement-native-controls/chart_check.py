#!/usr/bin/env python3
"""Exact complete Hepp-chart germ audit; no native or physical-flow admission."""
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def transform(monomial, permutation):
    # alpha[p0]=1, alpha[pj]=t1*...*tj, j>=1.
    return tuple(sum(monomial[permutation[j]] for j in range(i + 1, len(permutation)))
                 for i in range(len(permutation) - 1))


def audit(name, names, u_terms, v_terms, shifted):
    charts = []
    for perm in itertools.permutations(range(len(names))):
        u = [transform(m, perm) for m in u_terms]
        v = [transform(m, perm) for m in v_terms]
        sm = tuple(int(i == shifted) for i in range(len(names)))
        s = transform(sm, perm)
        valuation = tuple(min(m[i] for m in u) for i in range(len(names) - 1))
        unit = [tuple(x - y for x, y in zip(m, valuation)) for m in u]
        assert (0,) * len(valuation) in unit, (name, perm, 'U not a unit')
        ratio = [tuple(x - y - z for x, y, z in zip(m, valuation, s)) for m in v]
        assert all(e >= 0 for m in ratio for e in m), (name, perm, 'ratio not smooth', ratio)
        charts.append({'ordering': [names[i] for i in perm], 'U_valuation': valuation,
                       'S_valuation': s, 'U_positive_unit': unit, 'V_over_US_numerator': ratio,
                       'jacobian_powers': list(range(len(names) - 2, -1, -1)),
                       'U_unit_lower_bound': 1})
    return {'name': name, 'parameters': names, 'U_terms': u_terms, 'V_terms': v_terms,
            'shifted_parameter': names[shifted], 'chart_count': len(charts), 'charts': charts}


def main():
    out = HERE / 'chart-checks.json'
    assert not out.exists(), 'preserve previous reports'
    singleton = audit('cut3_shift0', list('abcd'),
                      [(1,0,1,0),(1,0,0,1),(0,1,1,0),(0,1,0,1),(0,0,1,1)],
                      [(1,1,1,0),(1,1,0,1),(1,0,1,1)], 0)
    bubble_x = audit('double_shift1_and2', list('xy'), [(1,0),(0,1)], [(1,1)], 0)
    bubble_y = audit('double_shift2_and4', list('xy'), [(1,0),(0,1)], [(1,1)], 1)
    report = {'schema': 'partial-placement-complete-germ-charts-v1',
              'scope': 'Independent exact Gaussian full-positive-support chart proof only. No native HEPKit replay, deleted-support inheritance, native discovery, or flow admission.',
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'families': [singleton, bubble_x, bubble_y],
              'complete_chart_count': 28,
              'result': 'Every U is a monomial times a positive unit, and V/(U*S) is a polynomial divided by that positive unit on every closed chart cube.'}
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'charts': 28, 'result': 'PASS'}))


if __name__ == '__main__':
    main()
