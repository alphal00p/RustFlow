#!/usr/bin/env python3
"""Audit saved proof-source images with exact rational polynomial arithmetic."""
import ast
import collections
import fractions
import hashlib
import json
from pathlib import Path

base = Path(__file__).resolve().parent
F = fractions.Fraction
def add(a, b):
    out = dict(a)
    for k, v in b.items(): out[k] = out.get(k, F(0)) + v
    return {k: v for k, v in out.items() if v}
def scale(a, c): return {k: v*c for k, v in a.items() if v*c}
def mul(a, b):
    out = {}
    for m, c in a.items():
        for n, d in b.items():
            powers = collections.Counter(dict(m)); powers.update(dict(n))
            k = tuple(sorted(powers.items()))
            out[k] = out.get(k, F(0)) + c*d
    return {k: v for k, v in out.items() if v}
def poly(text, fixed=None):
    fixed = fixed or {}
    def visit(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, int): return {(): F(n.value)} if n.value else {}
        if isinstance(n, ast.Name):
            if n.id in fixed: return {(): F(fixed[n.id])} if fixed[n.id] else {}
            return {((n.id, 1),): F(1)}
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return scale(visit(n.operand), -1)
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                assert isinstance(n.right, ast.Constant) and isinstance(n.right.value, int) and n.right.value >= 0
                a, out = visit(n.left), {(): F(1)}
                for _ in range(n.right.value): out = mul(out, a)
                return out
            a, b = visit(n.left), visit(n.right)
            if isinstance(n.op, ast.Add): return add(a, b)
            if isinstance(n.op, ast.Sub): return add(a, scale(b, -1))
            if isinstance(n.op, ast.Mult): return mul(a, b)
            if isinstance(n.op, ast.Div):
                assert len(b) == 1 and () in b and b[()]
                return scale(a, 1/b[()])
        raise AssertionError(ast.dump(n))
    return visit(ast.parse(text.replace('^', '**'), mode='eval').body)
def serial(p): return [[list(m), str(c)] for m, c in sorted(p.items())]
def condition_class(text):
    p = poly(text); assert p
    lead = p[min(p)]
    return json.dumps(serial(scale(p, 1/lead)), sort_keys=True)
def inside(label, domain):
    return all((lo is None or x >= lo) and (hi is None or x <= hi) for x, (lo, hi) in zip(label, domain))
def image_audit(data):
    source_rows = data['original_source_rows']
    out = []
    for ordinal, probe in enumerate(data['probes']):
        positives = []
        for ref in probe['rule_sources']:
            source = source_rows[ref['original_ordinal']]
            assert source['source_id'] == ref['source_id']
            seed = ref['seed']; assert inside(seed, source['domain'])
            fixed = {f'a_{i}': x for i, x in enumerate(seed)}
            for condition in source['conditions']: assert poly(condition, fixed)
            row = {}
            for term in source['terms']:
                label = tuple(seed[i] + shift if symbolic else shift for i, (symbolic, shift) in enumerate(term['label']))
                row[label] = add(row.get(label, {}), poly(term['coefficient'], fixed))
            for label, coefficient in sorted(row.items()):
                if not coefficient or label[6] <= 0: continue
                zero = None
                if label[0] <= 0: zero = 'required-cut nonpositive image'
                elif any(inside(label, domain) for domain in data['zero_domains']): zero = 'explicit certified zero-domain image'
                positives.append({'source_id': source['source_id'], 'seed': seed, 'indices': list(label),
                                  'exact_specialized_coefficient': serial(coefficient), 'native_zero_reason': zero})
        out.append({'point_ordinal': ordinal, 'target': probe['target'],
                    'positive_occupied_energy_seed_count': sum(r['seed'][6] > 0 for r in probe['rule_sources']),
                    'positive_occupied_energy_rhs': [t for t in probe['rhs'] if t['indices'][6] > 0],
                    'nonzero_coefficient_positive_source_images': positives,
                    'positive_source_images_after_known_zeros': sum(x['native_zero_reason'] is None for x in positives)})
    return out

data = {name: json.loads((base/(name+'.json')).read_text()) for name in ['polynomial62','reciprocal56']}
for value in data.values():
    assert value['roles'] == [1,0,0,0,0,0,0,0,0,2,2,0]
audits = {name: image_audit(value) for name, value in data.items()}
(base/'proof-energy-images.json').write_text(json.dumps({
    'scope': 'Exact coefficient specialization of the source rows named by every saved native rule proof. This audits proof-source image labels and final one-step RHS, not every transient internal elimination row. Required-cut and explicit certified zero images are reported separately. No additional integral relations are installed.',
    'arms': audits}, indent=2)+'\n')
comparisons = []
for i, (a,b) in enumerate(zip(data['polynomial62']['probes'],data['reciprocal56']['probes'])):
    assert a['target'] == b['target']
    assert a['persistence_roundtrip_replay'] and b['persistence_roundtrip_replay']
    comparisons.append({'ordinal':i,'target':a['target'],
                        'polynomial_status':a['status'],'reciprocal_status':b['status'],
                        'polynomial_rhs_terms':len(a['rhs']),'reciprocal_rhs_terms':len(b['rhs']),
                        'rhs_identical':a['rhs']==b['rhs'],
                        'raw_conditions_identical':a['conditions']==b['conditions'],
                        'conditions_same_up_to_nonzero_rational_scaling_and_repetition':set(map(condition_class,a['conditions']))==set(map(condition_class,b['conditions'])),
                        'polynomial_conditions':a['conditions'],'reciprocal_conditions':b['conditions'],
                        'different_rhs':None if a['rhs']==b['rhs'] else {'polynomial62':a['rhs'],'reciprocal56':b['rhs']}})
totals={}
for name, value in data.items():
    resource=json.loads((base/(name+'-resources.json')).read_text())
    audit=audits[name]
    totals[name]={'source_count':value['source_count'],'points':value['points'],
                  'applied':sum(x['status'].startswith('Applied') for x in value['probes']),
                  'unresolved':sum(x['status'].startswith('Unresolved') for x in value['probes']),
                  'rhs_terms':sum(len(x['rhs']) for x in value['probes']),
                  'proof_source_rows':sum(len(x['rule_sources']) for x in value['probes']),
                  'positive_energy_seeds':sum(x['positive_occupied_energy_seed_count'] for x in audit),
                  'positive_energy_rhs_terms':sum(len(x['positive_occupied_energy_rhs']) for x in audit),
                  'positive_energy_proof_images_before_zeros':sum(len(x['nonzero_coefficient_positive_source_images']) for x in audit),
                  'positive_energy_proof_images_after_known_zeros':sum(x['positive_source_images_after_known_zeros'] for x in audit),
                  'wall_seconds':resource['wall_seconds'],'peak_child_rss_kib':resource['peak_child_rss_kib']}
summary={'scope':'Deterministic 19-label sample of a later physical frontier; one native point domain, depth3, sample_seed0, max_domains1 per point. Fresh source contexts, no imported historical rules. No full closure or reciprocal necessity claim.',
         'totals':totals,'rhs_identical_points':sum(x['rhs_identical'] for x in comparisons),
         'conditions_equivalent_points':sum(x['conditions_same_up_to_nonzero_rational_scaling_and_repetition'] for x in comparisons),
         'comparisons':comparisons}
(base/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='comparisons'},indent=2))
