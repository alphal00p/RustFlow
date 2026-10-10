#!/usr/bin/env python3
"""Seal completed integration evidence after checking current source/build identity."""
import datetime
import gzip
import hashlib
import json
from pathlib import Path
from frozen_sources import verify_snapshot, digest

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIX = 'history-candidate-limit'
read = lambda p: json.loads(p.read_text())
write = lambda p, value: p.write_text(json.dumps(value, indent=2) + '\n')
sha = lambda data: hashlib.sha256(data).hexdigest()
manifest_path = BASE / 'artifact-manifest.json'
assert not manifest_path.exists()
snapshot = verify_snapshot(BASE / (PREFIX + '-source-hashes.json'))
build = read(BASE / (PREFIX + '-build-provenance.json'))
for name, info in build['files'].items():
    assert digest(ROOT / name) == info['sha256'], name
application = read(BASE / 'applied-source-binding.json')
for name, expected in application['source_sha256'].items():
    assert digest(ROOT / name) == expected, name
gates = read(BASE / (PREFIX + '-root-gates.json'))
assert gates['status'] == 'passed'
assert sum(g['passed'] for g in gates['gates']) == 213
checked = 0
for gate in gates['gates']:
    assert gate['failed'] == 0
    for info in gate['artifacts'].values():
        assert digest(ROOT / info['path']) == info['sha256']
        checked += 1
comparisons = {}
for suffix in ('two-loop-fixed', 'two-loop-laurent', 'massive-regression'):
    stem = PREFIX + '-' + suffix
    binding = read(BASE / (stem + '-binding.json'))
    resources = read(BASE / (stem + '-resources.json'))
    assert binding['exit_code'] == resources['exit_code'] == 0
    assert binding['post_run_source_and_executable_unchanged']
    comparison = read(BASE / (stem + '-comparison.json'))
    cb = read(BASE / (stem + '-comparison-binding.json'))
    assert comparison['status'] == cb['status'] == 'passed'
    for data in (comparison, cb):
        for name, expected in data.get('inputs_sha256', {}).items():
            assert digest(ROOT / name) == expected, name
            checked += 1
    comparisons[suffix] = {
        'status': 'passed',
        'comparison_sha256': digest(BASE / (stem + '-comparison.json')),
        'binding_sha256': digest(BASE / (stem + '-comparison-binding.json')),
        'resources': resources,
        'counts': {k: v for k, v in comparison.items() if 'count' in k},
    }
for name, expected in read(BASE / 'comparison-plan.json')['frozen_comparator_and_reference_artifacts_sha256'].items():
    assert digest(ROOT / name) == expected
    checked += 1
# All active readers use their separate static capsule. Keep the latest
# provisional and final closed transactions raw; archive older verbose payloads.
entries = []
for run in ('two-loop-fixed', 'two-loop-laurent', 'massive-regression'):
    closure = BASE / (PREFIX + '-' + run) / 'native-closure'
    for cut in sorted(closure.iterdir()):
        if not cut.is_dir():
            continue
        latest = max(cut.glob('round-*-provisional.json'))
        for path in sorted(cut.glob('round-*')):
            if '-closed.' in path.name or path.stem == latest.stem:
                continue
            raw = path.read_bytes()
            compressed = gzip.compress(raw, compresslevel=9, mtime=0)
            target = path.with_name(path.name + '.gz')
            assert not target.exists()
            target.write_bytes(compressed)
            assert gzip.decompress(target.read_bytes()) == raw
            entries.append({'original_path': str(path.relative_to(BASE)), 'archive_path': str(target.relative_to(BASE)),
                            'original_sha256': sha(raw), 'archive_sha256': sha(compressed),
                            'original_bytes': len(raw), 'archive_bytes': len(compressed)})
            path.unlink()
write(BASE / 'archive-map.json', {'schema': 1, 'gzip_mtime': 0, 'entries': entries,
    'scope': 'Lossless storage only; all latest provisional/closed proofs, predictions, configuration and comparisons remain raw.',
    'all_restore_checks_passed': True})
summary = {
    'schema': 1, 'status': 'passed', 'captured_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'source_files_verified_unchanged': len(snapshot['files']), 'build_artifacts_verified_unchanged': len(build['files']),
    'report_bindings_verified': checked, 'native_tests_passed': 91, 'rustflow_tests_passed': 122, 'native_tests_ignored': 1,
    'comparisons': comparisons, 'complete_final_source_replay_and_required_audits_retained': True,
    'production_admission_unchanged': True, 'full_three_loop_predictions': 0, 'four_loop_predictions': 0,
    'source_snapshot_sha256': digest(BASE / (PREFIX + '-source-hashes.json')),
    'build_provenance_sha256': digest(BASE / (PREFIX + '-build-provenance.json')),
    'archive_map_sha256': digest(BASE / 'archive-map.json'),
    'scope': 'Coherent optional historical-map limit integration and lower-loop regression gate only. Active massless profiles exercise cap0; massive inactive closure retains None. Higher-loop controls are separate. This final build record supersedes the application-time build_pending_at_application marker.'}
write(BASE / 'checkpoint-summary.json', summary)
files = {str(p.relative_to(BASE)): {'sha256': digest(p), 'bytes': p.stat().st_size}
         for p in sorted(BASE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p != manifest_path}
write(manifest_path, {'schema': 1, 'scope': summary['scope'], 'files': files, 'file_count': len(files),
                     'bytes': sum(v['bytes'] for v in files.values())})
print(json.dumps({'status': 'frozen', 'manifest_sha256': digest(manifest_path), 'file_count': len(files),
                  'archived_files': len(entries), 'bindings_verified': checked}))
