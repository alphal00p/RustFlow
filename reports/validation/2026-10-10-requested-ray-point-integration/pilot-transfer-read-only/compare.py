#!/usr/bin/env python3
"""Compare published production checkpoints with frozen source-only pilot JSON.

No native program is decoded, imported or modified. No numerical period or
reference values are read. Reported equality concerns scheduling and labels.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
P = ROOT / 'reports/validation/2026-10-10-native-source-portfolio-experiment'
K = ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(path):
    if not path.exists():
        path = path.with_suffix(path.suffix + '.gz')
    raw = path.read_bytes()
    payload = gzip.decompress(raw) if path.suffix == '.gz' else raw
    return json.loads(payload), {
        'path': str(path.relative_to(ROOT)),
        'bytes': len(raw), 'sha256': sha(raw),
        'decoded_json_sha256': sha(payload),
    }


def points(rows):
    return {tuple(row) for row in rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--last-round', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert 0 <= args.last_round < 16
    assert not args.output.exists(), 'Do not overwrite prior comparison evidence'
    production = K / 'requested-ray-point-three-loop-fixed/native-closure/cut-0'
    results = []
    for number in range(args.last_round + 1):
        pilot_dir = P / ('active-pilot' if number < 8 else 'active-resume') / 'polynomial62-baseline'
        old, old_binding = read(pilot_dir / f'round-{number:03}.json')
        new, new_binding = read(production / f'round-{number:03}-provisional.json')
        assert old['round'] == new['round'] == number
        state = new['active_state']
        transaction_path = production / state['requested_discovery']['path']
        transaction, transaction_binding = read(transaction_path)
        assert transaction['status'] == 'complete'
        assert transaction['schedule'] == 'requested-ray-point-v1'
        old_probes = {tuple(row['target']): row for row in old['discovery']}
        new_probes = {tuple(row['integral']): row for row in transaction['transactions']
                      if 'ray_probe' in row}
        differences = []
        for point in sorted(old_probes.keys() | new_probes.keys()):
            before, after = old_probes.get(point), new_probes.get(point)
            if before is None or after is None:
                differences.append({'point': point, 'field': 'search_present',
                                    'pilot': before is not None, 'production': after is not None})
                continue
            for old_field, new_field in [('ray_status', 'ray_probe'),
                                         ('point_status', 'point_probe')]:
                a, b = before.get(old_field), after.get(new_field, {}).get('status')
                if a != b:
                    differences.append({'point': point, 'field': old_field,
                                        'pilot': a, 'production': b})
        checks = {
            'frontier_labels': points(old['frontier']) == points(new['frontier']),
            'native_rule_count': old['rule_count'] == new['native_rule_count'],
            'submitted_labels': points(old['requested']) == points(state['submitted_requests']),
            'historical_labels': points(old['attempted_exact_points']) == points(state['historical_requests']),
            'next_needed_labels': points(old['next_needed']) == points(state['next_needed']),
            'searched_labels_and_application_statuses': not differences,
            'no_deferred_transactions': not transaction['deferred_requests'],
        }
        results.append({
            'round': number, 'checks': checks, 'all_checks_match': all(checks.values()),
            'counts': {'frontier': [len(old['frontier']), len(new['frontier'])],
                       'rules': [old['rule_count'], new['native_rule_count']],
                       'submitted': [len(old['requested']), len(state['submitted_requests'])],
                       'searched': [len(old_probes), len(new_probes)],
                       'allocated_domains': transaction['allocated_domains'],
                       'deferred': len(transaction['deferred_requests'])},
            'application_differences': differences,
            'inputs': [old_binding, new_binding, transaction_binding],
            'production_checkpoint_program_binding': {
                'blake3': new['program_blake3'], 'bytes': new['program_bytes']},
            'production_transaction_checkpoint_binding': state['requested_discovery'],
        })
    report = {
        'scope': 'Read-only scheduling/label comparison. Source equivalence is established by the separate production 62-row gate. Conditions/RHS coefficients are not asserted equal here. No native proof or numerical value was imported.',
        'script_sha256': sha(Path(__file__).read_bytes()),
        'last_completed_round': args.last_round,
        'all_compared_checks_match': all(row['all_checks_match'] for row in results),
        'rounds': results,
    }
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'rounds': len(results), 'all_compared_checks_match': report['all_compared_checks_match'],
                      'output': str(args.output), 'sha256': sha(args.output.read_bytes())}))


if __name__ == '__main__':
    main()
