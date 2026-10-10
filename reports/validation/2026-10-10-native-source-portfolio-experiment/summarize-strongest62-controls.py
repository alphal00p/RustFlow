#!/usr/bin/env python3
"""Bind matched portfolio controls to all source guards and the upstream corpus."""
import collections
import hashlib
import json
from pathlib import Path

base = Path(__file__).resolve().parent
root = base.parents[2]
upstream = root / 'reports/validation/2026-10-10-polynomial-raw-ward-attribution/native-source-probe/strongest62-boost-ward'
output = base / 'strongest62-point-summary.json'
assert not output.exists(), 'Preserve historical summaries.'
read = lambda name: json.loads((base / name).read_text())
b = read('baseline-strongest62-19-points.json')
q = read('proposal-strongest62-19-points.json')
old = read('baseline-original-19-points.json')
comparison = read('strongest62-portfolio-comparison.json')
n = json.loads((upstream / 'result.json').read_text())
keys = ['imported_source_schema', 'original_measure_id', 'roles',
        'index_variable_positions', 'zero_domains', 'original_source_rows']
assert all(b[key] == q[key] for key in keys)
assert b['imported_source_schema'] == 'rustred.guarded-source-program.v2'
assert b['source_count'] == q['source_count'] == 62
assert b['original_source_rows'][:52] == old['original_source_rows']
assert n['source_count'] == 62 and len(b['original_source_rows'][52:]) == 10
upstream_additions = {row['id']: row for key in ['additional_sources',
    'shifted_polynomial_global_boost_sources', 'shifted_temporal_U_sources'] for row in n[key]}
for actual in b['original_source_rows'][52:]:
    expected = upstream_additions[actual['source_id']]
    assert actual['source_id'] == expected['id']
    assert actual['domain'] == expected['domain'] and actual['conditions'] == expected['conditions']
    assert [{'coefficient': t['coefficient'], 'powers': t['label']}
            for t in actual['terms']] == [
                {'coefficient': t['coefficient'], 'powers': t['powers'] if 'powers' in t
                 else [[True, shift] for shift in t['shift']]} for t in expected['terms']]
    assert all(p[0] for t in actual['terms'] for p in t['label'])
upstream_points = {tuple(row['point']): row for row in n['probes']}
for row in b['probes']:
    known = upstream_points[tuple(row['target'])]
    assert row['rhs'] == known['terms'] and row['conditions'] == known['conditions']
assert comparison['points'] == comparison['baseline_applied'] == comparison['proposal_applied'] == 19
assert comparison['condition_gate_pass_points'] == 19
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
canonical_hash = lambda value: hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
arms = []
for index, label in enumerate(['strict protected powers', 'nonincreasing cut, index-dependent']):
    trials = [s['portfolio']['trials'][index] for row in q['probes']
              for s in row['search_stats'] if s['portfolio'] is not None]
    outcomes = collections.Counter(t['completion'] or 'not run after earlier exact zero' for t in trials)
    arms.append({'label': label, 'outcomes': dict(sorted(outcomes.items())),
                 **{key: sum(t[key] for t in trials) for key in
                    ['attempted_rows', 'accepted_rows', 'guard_rejected_rows', 'empty_rows',
                     'seeds', 'independent_rows', 'exact_trace_rows', 'exact_trace_terms']}})
files = [base / name for name in ['strongest62-portfolio-comparison.json',
    'baseline-original-19-points.json', 'points.json', 'generic-source-probe.rs',
    'run-generic-point-control.py', 'summarize-strongest62-controls.py',
    'compare-portfolio.py', 'condition-polynomials.py', 'focused-controls.json',
    'native-build/build-binding.json']]
resources = {}
for tag in ['baseline', 'proposal']:
    binding = read(tag + '-strongest62-19-points-binding.json')
    assert binding['exit_code'] == 0 and binding['post_run_inputs_unchanged']
    assert binding['source_corpus_sha256'] == sha(upstream / 'program.bin')
    assert binding['output_sha256'] == sha(base / (tag + '-strongest62-19-points.json'))
    resources[tag] = read(tag + '-strongest62-19-points-resources.json')
    files += [base / (tag + '-generic-probe-build.json')]
    files += [base / (tag + '-strongest62-19-points' + suffix) for suffix in
              ['.json', '-binding.json', '-resources.json', '-resources.log', '-resources.provenance.json']]
files += [upstream / 'program.bin', upstream / 'result.json']
report = {
    'schema': 1, 'status': 'isolated_point_controls_passed',
    'scope': 'Fresh native discovery and proof replay on the unchanged complete 62-source polynomial corpus. This is a local search experiment, not a closed system or physical numerical result.',
    'comparison': {key: comparison[key] for key in ['points', 'source_count', 'baseline_applied',
        'proposal_applied', 'baseline_rhs_terms', 'proposal_rhs_terms', 'strict_improvements', 'condition_gate_pass_points']},
    'source_binding_checks': {'all_complete_context_metadata_equal': True,
        'original_52_rows_byte_equal_as_decoded_json': True,
        'all_10_added_rows_match_upstream_terms_domains_and_conditions': True,
        'all_19_baseline_maps_match_upstream_native_result': True,
        'source_guards_and_zero_domains_unchanged_between_arms': True,
        'historical_rules_imported': False},
    'complete_context_canonical_json_sha256': {key: canonical_hash(b[key]) for key in keys},
    'improved_points': [{key: row[key] for key in ['target', 'baseline_rhs_terms', 'proposal_rhs_terms',
        'baseline_conditions', 'proposal_conditions']} for row in comparison['rows'] if row['strict_improvement']],
    'optional_arm_work': arms, 'resources': resources,
    'interpretation': 'The two certified improvements reduce 143 RHS terms to 127. Optional point search is slower in this control (about 0.59 to 1.55 seconds); only a separately bounded active-closure pilot can test overall benefit.',
    'limits_per_optional_arm': {'depth': 3, 'attempted_rows': 2048, 'accepted_rows': 512,
        'exact_trace_rows': 64, 'exact_trace_terms': 16384},
    'artifact_sha256': {str(path.relative_to(root)): sha(path) for path in files},
    'production_modified': False, 'closure_or_numerical_acceptance': False,
}
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report['comparison']))
