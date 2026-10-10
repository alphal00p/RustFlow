#!/usr/bin/env python3
"""Compare bounded native coverage/presentation only, never integral values."""
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
MODES = ['single-all', 'single-slot0', 'double-all', 'double-slots12', 'double-slots24']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(mode, report):
    probes = report['probes']
    assert len(probes) <= 16
    assert all(p['roundtrip'] for p in probes)
    return {'mode': mode, 'cuts': report['cuts'], 'shifted': report['shifted'],
            'sources': report['source_count'], 'points': len(probes),
            'applied': sum(p['application']['status'].startswith('Applied') for p in probes),
            'zero': sum(p['application']['status'] == 'Zero' for p in probes),
            'unresolved': sum(p['application']['status'].startswith('Unresolved') for p in probes),
            'rhs_terms': sum(len(p['application']['rhs']) for p in probes),
            'conditions': sum(len(p['application']['conditions']) for p in probes),
            'point_fallbacks': sum(p['point_fallback'] for p in probes),
            'discovery_gaps': sum(len(p['gaps']) for p in probes),
            'resources': json.loads((BASE / f'{mode}-resources.json').read_text()),
            'result_sha256': sha(BASE / mode / 'result.json')}


reports = {}
for mode in MODES:
    binding = json.loads((BASE / f'{mode}-binding.json').read_text())
    assert binding['exit_code'] == 0 and binding['inputs_unchanged']
    report = json.loads((BASE / mode / 'result.json').read_text())
    assert report['proof_binding_sha256'] == binding['proof_binding_sha256']
    assert report['limits'] == {'depth': 3, 'ray_domains': 1, 'point_domains': 1, 'max_selected_points': 16}
    reports[mode] = report

pairs = []
for baseline, alternative in [('single-all', 'single-slot0'), ('double-all', 'double-slots12'), ('double-all', 'double-slots24')]:
    a, b = reports[baseline], reports[alternative]
    for key in ['selected', 'roots', 'all_initial_derivatives', 'zero_domains', 'admitted_domain', 'source_ids', 'source_count', 'coordinates', 'inverse_routing', 'physical_arity', 'proof_binding_sha256']:
        assert a[key] == b[key], (baseline, alternative, key)
    assert a['source_identity'] != b['source_identity']
    changes = []
    for p, q in zip(a['probes'], b['probes']):
        assert p['point'] == q['point']
        if p['application'] != q['application']:
            changes.append({'point': p['point'], 'baseline': p['application'], 'alternative': q['application']})
    pairs.append({'baseline': baseline, 'alternative': alternative,
                  'identical_labels_domains_source_ids_and_limits': True,
                  'separate_deformation_bound_contexts': True, 'changed_outputs': changes})

out = BASE / 'summary.json'
assert not out.exists()
summary = {'scope': 'Bounded fresh native source-only comparison. Different masks define different finite-eta integrals; changed RHS is not an equality/accuracy test. Both members use restricted full-positive lower-origin zeros and omit raw Ward/free-virtual zeros. No numerical, derivative-closure or production-admission claim.',
           'rows': [summarize(m, reports[m]) for m in MODES], 'pairs': pairs,
           'roundtrip_gates': sum(len(r['probes']) for r in reports.values()),
           'inputs_sha256': {str(p): sha(p) for p in [BASE / 'compare.py', BASE / 'build-binding.json', BASE / 'chart-checks.json', BASE / 'origin-proof.md']}}
out.write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps([{'mode': r['mode'], 'sources': r['sources'], 'points': r['points'], 'applied': r['applied'], 'zero': r['zero'], 'unresolved': r['unresolved'], 'rhs_terms': r['rhs_terms']} for r in summary['rows']]))
