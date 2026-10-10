#!/usr/bin/env python3
"""Record passed integration/regression gates and both failed three-loop runs."""
import hashlib
import json
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIX = 'polynomial-closure-v2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = BASE / (PREFIX + '-checkpoint-summary.json')
    assert not output.exists()
    read = lambda path: json.loads(path.read_text())
    gate_path = BASE / (PREFIX + '-root-gates.json')
    native_path = BASE / (PREFIX + '-native-gate-reuse.json')
    lower_path = BASE / (PREFIX + '-lower-loop-regressions.json')
    snapshot = BASE / (PREFIX + '-source-hashes.json')
    build = BASE / (PREFIX + '-build-provenance.json')
    gates, native, lower = map(read, (gate_path, native_path, lower_path))
    assert gates['status'] == lower['status'] == 'passed'
    assert sum(row['passed'] for row in gates['gates']) == 105
    assert native['reused_passed_tests'] == 84 and native['new_native_test_invocations'] == 0
    failures = []
    bound = [gate_path, native_path, lower_path, snapshot, build, Path(__file__)]
    for suffix, passes in (('three-loop-fixed', 1), ('three-loop-fixed-guards3', 3)):
        stem = PREFIX + '-' + suffix
        directory = BASE / stem
        binding_path = BASE / (stem + '-binding.json')
        resources_path = BASE / (stem + '-resources.json')
        failure_path = directory / 'failure.json'
        binding, resources, failure = map(read, (binding_path, resources_path, failure_path))
        assert binding['exit_code'] == resources['exit_code'] == 101
        assert binding['source_snapshot_sha256'] == sha(snapshot)
        assert binding['build_provenance_sha256'] == sha(build)
        assert binding['post_run_source_and_executable_unchanged']
        assert binding['settings']['RUSTFLOW_WEIGHTED_GUARD_PASSES'] == str(passes)
        assert not failure['numerical_acceptance']
        assert not list(directory.glob('prediction-*.json'))
        rounds = sorted((directory / 'native-closure/cut-0').glob('round-*-provisional.json'))
        frontier = [len(read(path)['frontier']) for path in rounds]
        assert frontier == [7, 15, 105, 6]
        assert 'rounds=4, provisional=6, unresolved_terms=3' in failure['error']
        assert 'ConditionVanished { rule: 4218, condition: 0 }' in failure['error']
        refine_path = directory / 'native-closure/cut-0/conditional-point-refinements.json'
        refinements = read(refine_path)['refinements']
        assert len(refinements) == 1 and refinements[0]['pass'] == 1
        log_path = BASE / (stem + '-resources.log')
        times = re.findall(r'test result: FAILED\..*?finished in ([0-9.]+)s', log_path.read_text())
        assert len(times) == 1
        failures.append({
            'attempt': suffix, 'status': 'failed_preparation', 'guard_pass_budget': passes,
            'exit_code': 101, 'harness_seconds': float(times[0]),
            'wall_seconds': resources['wall_seconds'],
            'peak_child_rss_kib': resources['peak_child_rss_kib'],
            'attempted_cut': [0], 'completed_rounds': 4, 'frontier_progression': frontier,
            'unresolved_terms': 3, 'failure': failure['error'],
            'refinements_performed': refinements,
            'full_preparation_complete': False, 'predictions_saved': 0,
            'numerical_acceptance': False,
        })
        bound += [binding_path, resources_path, failure_path, refine_path, log_path, *rounds]
    assert failures[0]['failure'] == failures[1]['failure']
    report = {
        'schema': 1,
        'status': 'source_integration_and_lower_loop_regressions_passed_three_loop_blocked',
        'source_policy': 'polynomial-closure-v1', 'unchanged_native_rustred': True,
        'new_root_tests_passed': 105, 'native_tests_reused_by_exact_identity': 84,
        'new_native_test_invocations': 0,
        'actual_factory_source_equivalence': {'ordered_sources': 62, 'historical_rules_imported': False},
        'lower_loop_regressions': lower['gates'],
        'three_loop_attempts': failures,
        'three_loop_numerical_acceptance': False, 'four_loop_numerical_acceptance': False,
        'supplied_oracle_numerical_records_compared': 0,
        'scope': 'Optional exact source policy integration with independently compared two-loop regressions. Both full three-loop requests stop during first-singleton preparation at the same guarded exceptional child; neither saves numerical predictions. The source-only standalone closure is separate evidence and is not a full-amplitude prediction.',
        'bindings_sha256': {str(path.relative_to(ROOT)): sha(path) for path in bound},
    }
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'three_loop_attempts': len(failures)}))


if __name__ == '__main__':
    main()
