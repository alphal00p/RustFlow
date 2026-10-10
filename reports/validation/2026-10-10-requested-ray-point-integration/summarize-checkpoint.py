#!/usr/bin/env python3
"""Summarize completed gates, failed numerical requests and one native closure."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIX = 'requested-ray-point'
read = lambda path: json.loads(path.read_text())
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = BASE / (PREFIX + '-checkpoint-summary.json')
    assert not output.exists()
    gate_path = BASE / (PREFIX + '-root-gates.json')
    native_path = BASE / (PREFIX + '-native-gate-reuse.json')
    lower_path = BASE / (PREFIX + '-lower-loop-regressions.json')
    snapshot = BASE / (PREFIX + '-source-hashes.json')
    build = BASE / (PREFIX + '-build-provenance.json')
    verified = BASE / (PREFIX + '-final-source-verification.json')
    gates, native, lower = map(read, (gate_path, native_path, lower_path))
    assert gates['status'] == lower['status'] == 'passed'
    assert sum(row['passed'] for row in gates['gates']) == 115
    assert native['reused_passed_tests'] == 84 and native['new_native_test_invocations'] == 0
    attempts = []
    bound = [gate_path, native_path, lower_path, snapshot, build, verified, Path(__file__)]
    expected = [('three-loop-fixed', [3], 3, 1474, True),
                ('three-loop-double', [0, 3], 4, 2344, True),
                ('three-loop-double-residual-control', [0, 3], 7, 7346, False),
                ('three-loop-single3-residual-control', [3], 5, 2, False)]
    for suffix, cuts, count, unresolved, schedule in expected:
        stem = PREFIX + '-' + suffix
        directory = BASE / stem
        binding_path = BASE / (stem + '-binding.json')
        resources_path = BASE / (stem + '-resources.json')
        failure_path = directory / 'failure.json'
        config_path = directory / 'configuration.json'
        binding, resources, failure, config = map(read, (binding_path, resources_path, failure_path, config_path))
        assert binding['exit_code'] == resources['exit_code'] == 101
        assert binding['source_snapshot_sha256'] == sha(snapshot)
        assert binding['build_provenance_sha256'] == sha(build)
        assert binding['post_run_source_and_executable_unchanged']
        assert binding['settings']['RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT'] == str(schedule).lower()
        assert config['source_options']['policy'] == 'polynomial-closure-v1'
        assert config['independent_reference_comparisons'] == 0
        assert not failure['numerical_acceptance']
        assert not list(directory.glob('prediction-*.json'))
        corpus = directory / 'native-closure'
        if suffix == 'three-loop-fixed':
            corpus = corpus / 'cut-3'
        rounds = sorted(corpus.glob('round-*-provisional.json'))
        assert len(rounds) == count
        frontier = [len(read(path)['frontier']) for path in rounds]
        assert 'rounds=' + str(count) + ', provisional=' + str(frontier[-1]) + ', unresolved_terms=' + str(unresolved) in failure['error']
        log_path = BASE / (stem + '-resources.log')
        times = re.findall(r'test result: FAILED\..*?finished in ([0-9.]+)s', log_path.read_text())
        assert len(times) == 1
        attempts.append({'attempt': suffix, 'status': 'failed_preparation',
            'requested_ray_point': schedule, 'full_amplitude_requested': suffix == 'three-loop-fixed',
            'attempted_cut_at_failure': cuts, 'completed_rounds_at_failure': count,
            'frontier_progression': frontier, 'unresolved_terms': unresolved,
            'failure_stage': failure['stage'], 'failure': failure['error'],
            'exit_code': 101, 'harness_seconds': float(times[0]),
            'wall_seconds': resources['wall_seconds'], 'peak_child_rss_kib': resources['peak_child_rss_kib'],
            'predictions_saved': 0, 'numerical_acceptance': False})
        bound += [binding_path, resources_path, failure_path, config_path, directory / 'input.json', log_path, *rounds]

    # Integrity of the completed cut-0 native proof, independently of a failed
    # whole-amplitude preparation. This does not evaluate a physical period.
    tool_record_path = BASE / 'native-digest-build.json'
    tool_record = read(tool_record_path)
    exe = Path(tool_record['executable']['path'])
    assert sha(exe) == tool_record['executable']['sha256']
    def b3(path):
        words = subprocess.check_output([str(exe), str(path)], text=True).split()
        assert len(words) == 2 and int(words[1]) == path.stat().st_size
        return words[0]
    folder = BASE / (PREFIX + '-three-loop-fixed/native-closure/cut-0')
    closed = folder / 'round-015-closed.json'
    row = read(closed)
    assert row['status'] == 'closed' and row['round'] == 15 and len(row['frontier']) == 7
    assert row['discovery_schedule'] == 'requested-ray-point-v1'
    proof = closed.with_suffix('.bin')
    assert proof.stat().st_size == row['program_bytes'] and b3(proof) == row['program_blake3']
    link = row['active_state']['requested_discovery']
    assert re.fullmatch(r'requested-discovery-[0-9]+[.]json', link['path'])
    transaction_path = folder / link['path']
    transaction = read(transaction_path)
    assert b3(transaction_path) == link['blake3']
    assert transaction['schema'] == 1 and transaction['status'] == 'complete'
    assert transaction['schedule'] == row['discovery_schedule']
    assert transaction['source_measure_id'] == row['measure_id']
    for key in ('physical_arity', 'storage_capacity'):
        assert transaction[key] == row[key]
    assert transaction['policy'] == row['requested_ray_point']
    assert transaction['empty_terminals'] is True and transaction['completed_search_implies_coverage'] is False
    provisional = folder / 'round-015-provisional.json'
    p = read(provisional)
    provisional_proof = provisional.with_suffix('.bin')
    assert p['active_state'] == row['active_state'] and p['round'] == row['round']
    assert provisional_proof.stat().st_size == p['program_bytes'] == transaction['native_program']['bytes']
    assert b3(provisional_proof) == p['program_blake3'] == transaction['native_program']['blake3']
    cut0_rounds = sorted(folder.glob('round-*-provisional.json'))
    assert len(cut0_rounds) == 16
    cut0 = {'cut_slots': [0], 'status': 'native_connection_closed', 'completed_rounds': 16,
        'basis_size': 7, 'native_rule_count': row['native_rule_count'],
        'frontier_progression': [len(read(path)['frontier']) for path in cut0_rounds],
        'closed_program_blake3': row['program_blake3'], 'closed_program_bytes': row['program_bytes'],
        'closed_program_sha256': sha(proof), 'transaction_and_provisional_digests_verified': True,
        'raw_protected_paths': [str(path.relative_to(BASE)) for path in (closed, proof, provisional, provisional_proof, transaction_path)],
        'scope': 'One occupied native connection only. Full preparation fails later; no three-loop numerical value is accepted.'}
    bound += [tool_record_path, closed, proof, provisional, provisional_proof, transaction_path, *cut0_rounds]
    report = {'schema': 1,
        'status': 'requested_schedule_integration_and_lower_loop_regressions_passed_three_loop_blocked',
        'source_policy': 'polynomial-closure-v1', 'unchanged_native_rustred': True,
        'new_root_tests_passed': 115, 'native_tests_reused_by_exact_identity': 84,
        'new_native_test_invocations': 0,
        'actual_factory_source_equivalence': {'ordered_sources': 62, 'historical_rules_imported': False},
        'lower_loop_regressions': lower['gates'], 'three_loop_attempts': attempts,
        'completed_cut0_connection': cut0,
        'three_loop_laurent_run_performed': False, 'three_loop_numerical_acceptance': False,
        'four_loop_numerical_acceptance': False, 'supplied_oracle_numerical_records_compared': 0,
        'scope': 'Passed default-branch two-loop regressions and optional schedule integration tests. Four three-loop requests fail before any prediction. A completed cut-0 native connection does not certify a full amplitude. Timings overlap other work; no performance comparison is claimed. Reference precision is empirical refinement evidence, not a rigorous bound.',
        'bindings_sha256': {str(path.relative_to(ROOT)): sha(path) for path in sorted(set(bound))}}
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'three_loop_failures': len(attempts), 'cut0_closed_basis': 7}))

if __name__ == '__main__':
    main()
