#!/usr/bin/env python3
"""Seal the completed routing diagnostic without changing its failed result."""
import gzip
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
K = ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration'
read = lambda p: json.loads(p.read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    assert not (HERE/'artifact-manifest.json').exists()
    binding = read(HERE/'run-binding.json')
    resources = read(HERE/'resources.json')
    failure = read(HERE/'selected-cut/failure.json')
    audit = read(HERE/'routing-audit.json')
    assert binding['exit_code'] == resources['exit_code'] == 101
    assert binding['post_run_source_and_executable_unchanged']
    assert not failure['numerical_acceptance']
    assert not list((HERE/'selected-cut').glob('prediction-*.json'))
    assert sha(HERE/'run.py') == binding['launcher_sha256']
    assert sha(HERE/'rerouted-input.json') == binding['rerouted_input_sha256']
    assert read(HERE/'rerouted-input.json') == read(HERE/'selected-cut/input.json')
    assert sha(ROOT/binding['original_input']) == binding['original_sha256']
    assert sha(K/'requested-ray-point-build-provenance.json') == binding['build_provenance_sha256']
    assert sha(K/'requested-ray-point-source-hashes.json') == binding['source_snapshot_sha256']
    for name, digest in audit['inputs_sha256'].items(): assert sha(ROOT/name) == digest
    corpus = HERE/'selected-cut/native-closure'
    rounds = sorted(corpus.glob('round-*-provisional.json'))
    frontier = [len(read(p)['frontier']) for p in rounds]
    assert frontier == [20, 124, 603, 984, 1025]
    assert 'rounds=5, provisional=1025, unresolved_terms=1054' in failure['error']
    log = HERE/'resources.log'
    times = re.findall(r'test result: FAILED\..*?finished in ([0-9.]+)s', log.read_text())
    assert len(times) == 1
    baseline_path = K/'requested-ray-point-three-loop-fixed-binding.json'
    baseline = read(baseline_path)
    changes = {key: {'before': baseline['settings'].get(key), 'after': binding['settings'].get(key)}
        for key in sorted(set(baseline['settings']) | set(binding['settings']))
        if baseline['settings'].get(key) != binding['settings'].get(key)}
    bound = [HERE/'run-binding.json', HERE/'resources.json', HERE/'resources.provenance.json', log,
        HERE/'routing-audit.json', HERE/'check_routing.py', HERE/'run.py', HERE/'rerouted-input.json',
        HERE/'selected-cut/input.json', HERE/'selected-cut/configuration.json', HERE/'selected-cut/sector.json',
        HERE/'selected-cut/failure.json', *rounds, baseline_path,
        K/'requested-ray-point-build-provenance.json', K/'requested-ray-point-source-hashes.json']
    summary = {'status': 'bounded_native_closure_failed', 'exact_routing_audit': 'passed',
        'original_cut': [3], 'rerouted_cut': [0], 'original_raised_slot': 0, 'rerouted_raised_slot': 3,
        'determinant': -1, 'absolute_jacobian': 1, 'additional_cut_or_wick_phase': False,
        'frontier_progression': frontier, 'completed_rounds': 5, 'unresolved_terms': 1054,
        'failure': failure, 'exit_code': 101, 'harness_seconds': float(times[0]),
        'wall_seconds': resources['wall_seconds'], 'peak_child_rss_kib': resources['peak_child_rss_kib'],
        'settings_changes_from_full_K': changes, 'selected_test': 'runtime_graph_occupied_flow',
        'timeout_seconds': 600, 'original_full_timeout_seconds': 3600,
        'source_or_reducer_changes': False, 'historical_rules_imported': False,
        'numerical_predictions': 0, 'full_amplitude_acceptance': False,
        'scope': 'Exact graph/input equivalence and a failed selected-cut native preparation. Routing, physical-edge order and completion coordinates change together. No isolated-completion, numerical-equivalence or speed conclusion.',
        'bindings_sha256': {str(p.relative_to(ROOT)): sha(p) for p in bound}}
    (HERE/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    preserve = set()
    for number in (0, 4):
        path = corpus/f'round-{number:03}-provisional.json'
        preserve.update((path, path.with_suffix('.bin')))
        preserve.add(corpus/read(path)['active_state']['requested_discovery']['path'])
    preserved = {str(p.relative_to(HERE)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(preserve)}
    entries = []
    for path in sorted([p for p in corpus.rglob('*') if p.is_file() and p not in preserve] + [log]):
        raw = path.read_bytes(); target = path.with_name(path.name+'.gz'); assert not target.exists()
        encoded = gzip.compress(raw, compresslevel=9, mtime=0)
        target.write_bytes(encoded)
        assert gzip.decompress(target.read_bytes()) == raw
        entries.append({'original_path': str(path.relative_to(HERE)), 'archive_path': str(target.relative_to(HERE)),
            'original_sha256': hashlib.sha256(raw).hexdigest(), 'archive_sha256': sha(target),
            'original_bytes': len(raw), 'archive_bytes': len(encoded)})
        path.unlink()
    archive = {'scope': 'Completed proof payloads and resource log; readers released. Original/rerouted inputs, bindings, failure and round000/004 proof pairs plus their transactions remain raw.',
        'gzip_mtime': 0, 'entries': entries, 'preserved_raw': preserved, 'restore_checks_passed': True,
        'file_count': len(entries), 'original_bytes': sum(r['original_bytes'] for r in entries),
        'archive_bytes': sum(r['archive_bytes'] for r in entries)}
    (HERE/'archive-map.json').write_text(json.dumps(archive, indent=2)+'\n')
    mapping = {r['original_path']:r for r in entries}
    for name, digest in summary['bindings_sha256'].items():
        path = ROOT/name
        raw = path.read_bytes() if path.exists() else gzip.decompress((HERE/mapping[str(path.relative_to(HERE))]['archive_path']).read_bytes())
        assert hashlib.sha256(raw).hexdigest() == digest
    for name, row in preserved.items(): assert sha(HERE/name) == row['sha256']
    files = {str(p.relative_to(HERE)): {'sha256':sha(p), 'bytes':p.stat().st_size}
        for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
    manifest = {'scope': summary['scope'], 'status': summary['status'], 'files':files,
        'file_count':len(files), 'all_historical_bindings_and_archives_verified':True}
    (HERE/'artifact-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':'frozen','files':len(files),'archives':len(entries),
        'raw_bytes':archive['original_bytes'],'gzip_bytes':archive['archive_bytes'],
        'manifest_sha256':sha(HERE/'artifact-manifest.json')}))

if __name__ == '__main__': main()
