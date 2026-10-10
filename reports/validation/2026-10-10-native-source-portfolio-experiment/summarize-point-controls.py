#!/usr/bin/env python3
"""Summarize preserved, oracle-free native point controls without rerunning them."""
import collections
import hashlib
import json
from pathlib import Path

base = Path(__file__).resolve().parent
output = base / 'original-19-point-summary.json'
assert not output.exists(), 'Preserve historical summaries.'
read = lambda name: json.loads((base / name).read_text())
comparison = read('original-19-portfolio-comparison.json')
proposal = read('proposal-original-19-points.json')
focused = read('focused-controls.json')
assert focused['passed'] == 19 and focused['failed'] == 0
arms = []
for index, label in enumerate(['strict protected powers', 'nonincreasing cut, index-dependent']):
    trials = [s['portfolio']['trials'][index] for row in proposal['probes']
              for s in row['search_stats'] if s['portfolio'] is not None]
    outcomes = collections.Counter(t['completion'] or 'not run after earlier exact zero' for t in trials)
    arms.append({'label': label, 'outcomes': dict(sorted(outcomes.items())),
                 **{key: sum(t[key] for t in trials) for key in
                    ['attempted_rows', 'accepted_rows', 'guard_rejected_rows', 'empty_rows',
                     'seeds', 'independent_rows', 'exact_trace_rows', 'exact_trace_terms']}})
files = ['original-19-portfolio-comparison.json', 'focused-controls.json',
         'focused-controls-binding.json', 'focused-controls-resources.json',
         'isolated-source-map.json', 'isolated-native.patch', 'native-build/build-binding.json',
         'portfolio-probe.rs', 'compare-portfolio.py', 'condition-polynomials.py', 'points.json',
         'summarize-point-controls.py']
resources = {}
cross = {}
for tag in ['baseline', 'proposal']:
    resources[tag] = read(tag + '-original-19-points-resources.json')
    r = read('cross-schema-' + tag + '-resources.json')
    log = (base / ('cross-schema-' + tag + '-resources.log')).read_text()
    assert r['exit_code'] == 0 and '"cross_schema_rejected":true' in log
    cross[tag] = {'exit_code': 0, 'incompatible_guarded_program_rejected': True}
    files += [tag + '-original-19-points' + suffix for suffix in
              ['.json', '-binding.json', '-resources.json', '-resources.log', '-resources.provenance.json']]
    files += ['cross-schema-' + tag + '-resources' + suffix for suffix in
              ['.json', '.log', '.provenance.json']]
report = {
    'schema': 1, 'status': 'isolated_point_controls_passed',
    'scope': 'Fresh native discovery from the complete original 52-source corpus, exact replay and persistence roundtrip. No historical rules imported; no production modification, active closure or physical amplitude acceptance.',
    'focused_controls_passed': focused['passed'],
    'comparison': {k: comparison[k] for k in ['points', 'source_count', 'baseline_applied',
        'proposal_applied', 'baseline_rhs_terms', 'proposal_rhs_terms', 'strict_improvements',
        'condition_gate_pass_points']},
    'improved_points': [{k: row[k] for k in ['target', 'baseline_rhs_terms', 'proposal_rhs_terms',
        'baseline_conditions', 'proposal_conditions']} for row in comparison['rows'] if row['strict_improvement']],
    'optional_arm_work': arms,
    'cross_schema_controls': cross,
    'resources': resources,
    'cost_interpretation': 'Optional searches increase the measured local point-control time. These small controls establish two shorter certified rules, not a net closure or runtime advantage.',
    'limits_per_optional_arm': {'depth': 3, 'attempted_rows': 2048,
        'accepted_rows': 512, 'exact_trace_rows': 64, 'exact_trace_terms': 16384},
    'condition_scope': 'Every admitted alternative condition is a nonzero constant or an exact nonzero rational multiple of an original sealed baseline condition; original conditions are retained unchanged.',
    'failure_semantics': 'Budget, exhaustion and alternative certification rejection retain baseline. Fatal arithmetic, sample and exact-lift failures propagate. Ordinary and ray paths are disabled.',
    'artifact_sha256': {name: hashlib.sha256((base / name).read_bytes()).hexdigest() for name in files},
    'production_modified': False, 'closure_or_numerical_acceptance': False,
}
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report['comparison']))
