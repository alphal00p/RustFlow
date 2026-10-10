#!/usr/bin/env python3
"""Compare complete saved native active closures; never evaluate periods."""
import collections
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path

base = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('conditions', base/'condition-polynomials-v2.py')
conditions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(conditions)
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
metadata = lambda path: {'path':str(path.relative_to(base)), 'bytes':path.stat().st_size, 'sha256':sha(path)}
results = {arm:json.loads((base/'active-resume'/f'polynomial62-{arm}'/'result.json').read_text()) for arm in ['baseline','proposal']}
b, p = results['baseline'], results['proposal']
assert b['status'] == p['status'] == 'closed'
assert b['original_targets'] == p['original_targets']
assert b['closed']['basis'] == p['closed']['basis']
assert b['closed']['target_then_derivative_rows'] == p['closed']['target_then_derivative_rows']
assert len(b['original_targets']) == 2 and sum(map(len,b['original_targets'])) == 3
assert len(b['closed']['basis']) == 7 and len(b['closed']['target_then_derivative_rows']) == 9
assert b['measure_id'] == p['measure_id']
limits = ['max_frontier','max_requested','max_rounds','max_rules','depth','source_count',
          'zero_attempts_per_round','ray_domains_per_requested_point','point_domains_per_fallback',
          'shifted_ordinary_slots']
assert all(b[key] == p[key] for key in limits)
groups, summaries, guard_evaluations = {}, {}, {}
for arm, result in results.items():
    directory = base/'active-resume'/f'polynomial62-{arm}'
    parent = base/'active-pilot'/f'polynomial62-{arm}'
    assert result['conditions_encountered'] == result['closed']['conditions']
    assert len(result['conditions_encountered']) == len(set(result['conditions_encountered']))
    groups[arm] = collections.defaultdict(list)
    guard_evaluations[arm] = []
    for expression in result['closed']['conditions']:
        polynomial = conditions.parse(expression)
        assert polynomial, 'zero retained guard'
        canonical = conditions.canonical(polynomial)
        groups[arm][canonical].append(expression)
        evaluated = {}
        for monomial, coefficient in polynomial.items():
            powers = dict(monomial)
            assert set(powers) <= {'epsilon','eta'}
            coefficient *= Fraction(-5,4)**powers.get('epsilon',0)
            eta_power = powers.get('eta',0)
            evaluated[eta_power] = evaluated.get(eta_power,Fraction(0)) + coefficient
        evaluated = {k:v for k,v in evaluated.items() if v}
        assert len(evaluated) == 1, 'guard is not proved nonzero on positive eta by this audit'
        power, coefficient = next(iter(evaluated.items()))
        assert coefficient and power >= 0
        guard_evaluations[arm].append({'condition':expression,'epsilon':'-5/4',
                                       'coefficient':str(coefficient),'eta_power':power,
                                       'nonzero_for_every_eta_positive':True})
    resource = json.loads((directory/'resources.json').read_text())
    parent_resource = json.loads((parent/'resources.json').read_text())
    binding = json.loads((directory/'run-binding.json').read_text())
    assert resource['exit_code'] == binding['exit_code'] == 0
    assert parent_resource['exit_code'] == 124
    assert not binding['original_corpus_rules_imported']
    assert binding['parent_discovered_rules_restored_with_native_replay']
    cumulative = resource['wall_seconds'] + parent_resource['wall_seconds']
    assert abs(cumulative-binding['cumulative_wall_seconds']) < 1e-8
    for name, digest in binding['frozen_parent_and_build_artifact_sha256'].items():
        path = Path(name)
        assert path.exists() and sha(path) == digest, name
    summaries[arm] = {
        'status':result['status'], 'basis_size':len(result['closed']['basis']),
        'completed_rounds':len(result['history']), 'native_rules':result['history'][-1]['rules'],
        'frontier_history':[x['frontier'] for x in result['history']],
        'max_frontier_observed':max(x['frontier'] for x in result['history']),
        'raw_retained_conditions':len(result['closed']['conditions']),
        'rational_associate_classes':len(groups[arm]),
        'source_zero_rules_added':sum(x['zero_rules'] for x in result['history']),
        'parent_timeout_exit_code':124,'parent_wall_seconds':parent_resource['wall_seconds'],
        'resume_exit_code':0,'resume_wall_seconds':resource['wall_seconds'],
        'cumulative_process_wall_seconds':cumulative,
        'parent_peak_child_rss_kib':parent_resource['peak_child_rss_kib'],
        'resume_peak_child_rss_kib':resource['peak_child_rss_kib'],
        'inputs':[metadata(directory/x) for x in ['result.json','closed.bin','run-binding.json','resources.json']],
    }

def serial_class(canonical):
    return [{'powers':list(monomial),'coefficient':str(coefficient)} for monomial,coefficient in canonical]
def class_records(keys):
    return [{'canonical_polynomial':serial_class(key),
             'baseline_raw_conditions':groups['baseline'].get(key,[]),
             'proposal_raw_conditions':groups['proposal'].get(key,[])} for key in sorted(keys)]
bs, ps = set(groups['baseline']), set(groups['proposal'])
condition_report = {
    'scope':'Only exact polynomial equality up to nonzero rational scaling and repetition. No factorization, product splitting, root-set equivalence or condition deletion. Every original condition is preserved.',
    'shared_class_count':len(bs & ps),'baseline_only_class_count':len(bs-ps),'proposal_only_class_count':len(ps-bs),
    'shared_classes':class_records(bs & ps),'baseline_only_classes':class_records(bs-ps),
    'proposal_only_classes':class_records(ps-bs),'exact_epsilon_minus_5_over_4_checks':guard_evaluations,
    'sample_scope':'At epsilon=-5/4 every stored guard is a nonzero rational times a nonnegative integer power of eta. This proves guard admissibility for eta>0 only; it evaluates no integral and does not license eta=0.'}
(base/'active-closure-condition-comparison.json').write_text(json.dumps(condition_report,indent=2)+'\n')
summary = {
    'status':'both_native_active_closures_passed',
    'scope':'One occupied singleton channel of the three-loop massless chain; both original weighted output sums and derivatives of every member of the discovered seven-element spanning basis pass the final native audit. No complete finite-density amplitude, basis minimality or numerical period claim.',
    'shared_limits':{key:b[key] for key in limits},'arms':summaries,
    'same_exact_ordered_basis':True,'same_exact_original_target_rows':True,'same_exact_basis_derivative_rows':True,
    'same_retained_condition_classes':bs==ps,
    'original_targets':b['original_targets'],'basis':b['closed']['basis'],
    'target_then_derivative_rows':b['closed']['target_then_derivative_rows'],
    'native_audit_contract':'The saved runner rebinds exactly the final basis terminals through native replay, reduces the actual two original weighted sums and all seven actual derivative sums, and asserts no failures and no unresolved leaves before writing closed.bin.',
    'conditions':{'shared_classes':len(bs&ps),'baseline_only_classes':len(bs-ps),'proposal_only_classes':len(ps-bs),
                  'both_sets_nonzero_at_epsilon_minus_5_over_4_and_eta_positive':True},
    'timing_scope':'Parents and continuations ran on an overlapping shared host. Cumulative values include timeout work and repeated unfinished-round work; exclude compilation and Nix startup. No performance superiority claim.',
    'production_native_changed':False,'numerical_periods_evaluated':False,
    'analysis_inputs':[metadata(Path(__file__)),metadata(base/'condition-polynomials-v2.py'),metadata(base/'active-resume/pilot.rs')],
}
(base/'active-closure-comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'status':summary['status'],'arms':summaries,'conditions':summary['conditions']},indent=2))
