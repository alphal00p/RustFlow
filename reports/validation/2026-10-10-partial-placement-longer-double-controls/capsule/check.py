#!/usr/bin/env python3
"""Exhaustive active-support proof replay using exact rational Gaussian rows."""
import hashlib
import importlib.util
import itertools
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
OLD = HERE.parent / '2026-10-10-partial-placement-native-controls/chart_check.py'
spec = importlib.util.spec_from_file_location('charts', OLD)
charts = importlib.util.module_from_spec(spec); spec.loader.exec_module(charts)


def rank(rows, width):
    r = [list(map(Q, row)) for row in rows]
    pivot = 0
    for c in range(width):
        p = next((i for i in range(pivot, len(r)) if r[i][c]), None)
        if p is None:
            continue
        r[pivot], r[p] = r[p], r[pivot]
        scale = r[pivot][c]; r[pivot] = [v / scale for v in r[pivot]]
        for i in range(len(r)):
            if i != pivot:
                scale = r[i][c]; r[i] = [a - scale*b for a, b in zip(r[i], r[pivot])]
        pivot += 1
    return pivot


def restrict(terms, active):
    return [tuple(t[i] for i in active) for t in terms
            if all(not t[i] or i in active for i in range(len(t)))]


def family(name, symbols, slots, rows, external, u, v, shifted, full_rank):
    supports = []
    for flags in itertools.product([False, True], repeat=len(symbols)):
        active = [i for i, flag in enumerate(flags) if flag]
        routed = [rows[i] for i in active]
        r = rank(routed, full_rank)
        record = {'positive_slots': [slots[i] for i in active],
                  'nonpositive_slots': [slots[i] for i in range(len(symbols)) if i not in active],
                  'virtual_rank': r}
        if r < full_rank:
            witnesses = [(1, 0), (0, 1), (1, 1)] if full_rank == 2 else [(1,)]
            null = next(w for w in witnesses if all(sum(a*b for a,b in zip(row,w)) == 0 for row in routed))
            record.update(classification='unrestricted polynomial virtual direction zero', null_vector=list(null))
        elif shifted not in active:
            # Exact affine translation x -> x + c*q makes all positive
            # denominators homogeneous massless and independent of compact q.
            candidates = [(0,0)] if full_rank == 2 else [(0,), (1,), (-1,)]
            translation = next(c for c in candidates if all(external[i] + sum(a*b for a,b in zip(rows[i],c)) == 0 for i in active))
            record.update(classification='translated homogeneous massless virtual vacuum zero',
                          virtual_translation=list(translation),
                          positive_virtual_masses=['0']*len(active),
                          absent_shifted_factor_is_polynomial=True)
        else:
            uu, vv = restrict(u, active), restrict(v, active)
            assert uu and all(any(t) or full_rank == 0 for t in uu)
            checked = charts.audit(name, [symbols[i] for i in active], uu, vv, active.index(shifted))
            record.update(classification='full-rank shifted-positive holomorphic invariant germ',
                          U_terms=uu, V_terms=vv, charts=checked['charts'],
                          chart_count=checked['chart_count'], V_is_zero=not vv)
        supports.append(record)
    return {'name': name, 'all_virtual_physical_slots': slots, 'shifted_virtual_slot': slots[shifted],
            'routed_virtual_rows': rows, 'external_coefficients': external, 'supports': supports}


def main():
    # q=P1-P3, t=P3, l=P2; rows ordered slots0,2,1,4.
    singleton = family('single_cut3_shift0', list('abcd'), [0,2,1,4],
                       [(1,0),(1,0),(0,1),(-1,1)], [1,0,0,0],
                       [(1,0,1,0),(1,0,0,1),(0,1,1,0),(0,1,0,1),(0,0,1,1)],
                       [(1,1,1,0),(1,1,0,1),(1,0,1,1)], 0, 2)
    double_x = family('double_cut03_shift12', list('xy'), [1,4], [(1,),(1,)], [0,-1],
                      [(1,0),(0,1)], [(1,1)], 0, 1)
    double_y = family('double_cut03_shift24', list('xy'), [1,4], [(1,),(1,)], [0,-1],
                      [(1,0),(0,1)], [(1,1)], 1, 1)
    families = [singleton,double_x,double_y]
    summary=[]
    for f in families:
        labels=[s['classification'] for s in f['supports']]
        summary.append({'name':f['name'],'supports':len(labels),
                        'rank_deficient_zeros':sum(x.startswith('unrestricted') for x in labels),
                        'homogeneous_vacuum_zeros':sum(x.startswith('translated') for x in labels),
                        'germs':sum(x.startswith('full-rank') for x in labels),
                        'complete_germ_charts':sum(s.get('chart_count',0) for s in f['supports'])})
    result={'schema':'partial-placement-active-support-audit-v1',
            'scope':'Independent exact support/routing/chart proof. No HEPKit replay, production permit, native rule discovery, period formula, or eta endpoint claim.',
            'integer_indices':'Positive propagators define the actual support. Every nonpositive physical/completion index contributes a finite polynomial numerator; it is not analytically continued from a positive-support Schwinger denominator here.',
            'families':families,'summary':summary,
            'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD,Path(__file__)]}}
    out=HERE/'support-checks.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':
    main()
