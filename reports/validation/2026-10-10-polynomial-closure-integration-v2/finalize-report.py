#!/usr/bin/env python3
"""Verify lossless archives and frozen report bindings, then seal this attempt."""
import gzip
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIX = 'polynomial-closure-v2'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    output = BASE / 'artifact-manifest.json'
    verification = BASE / 'archive-and-binding-verification.json'
    assert not output.exists() and not verification.exists()
    read = lambda path: json.loads(path.read_text())
    archive_path = BASE / (PREFIX + '-native-archives.json')
    archives = read(archive_path)
    mapping = {row['original_path']: row for row in archives['entries']}
    for row in archives['entries']:
        compressed = (BASE / row['archive_path']).read_bytes()
        assert digest(compressed) == row['archive_sha256']
        raw = gzip.decompress(compressed)
        assert digest(raw) == row['original_sha256'] and len(raw) == row['original_bytes']
        assert not (BASE / row['original_path']).exists()
    for name, row in archives['preserved_raw'].items():
        assert digest((BASE / name).read_bytes()) == row['sha256']

    def resolve(path):
        if path.exists():
            return path.read_bytes()
        name = str(path.relative_to(BASE))
        return gzip.decompress((BASE / mapping[name]['archive_path']).read_bytes())

    checked = 0
    summary_names = [PREFIX + '-checkpoint-summary.json', PREFIX + '-domain-retry-addendum.json',
                     PREFIX + '-per-residual-control-addendum.json']
    for name in summary_names:
        for file, expected in read(BASE / name)['bindings_sha256'].items():
            assert digest(resolve(ROOT / file)) == expected
            checked += 1
    lower = read(BASE / (PREFIX + '-lower-loop-regressions.json'))
    for row in lower['gates']:
        for info in row['bindings'].values():
            assert digest(resolve(ROOT / info['path'])) == info['sha256']
            checked += 1
    for path in BASE.glob('*-comparison-binding.json'):
        data = read(path)
        assert data['status'] == 'passed'
        for file, expected in data['inputs_sha256'].items():
            assert digest(resolve(ROOT / file)) == expected
            checked += 1
    diagnostic_manifests = {}
    for name in ('conditional-chain-diagnostic', 'conditional-chain-domains16384-diagnostic'):
        path = BASE / name / 'manifest.json'
        data = read(path)
        for row in data['files']:
            assert digest((ROOT / row['path']).read_bytes()) == row['sha256']
        diagnostic_manifests[name] = {'sha256': digest(path.read_bytes()),
                                      'verified_files': len(data['files'])}
    result = {'status': 'passed', 'lossless_archives_verified': len(mapping),
              'preserved_raw_files_verified': len(archives['preserved_raw']),
              'historical_report_bindings_verified': checked,
              'diagnostic_manifests': diagnostic_manifests,
              'scope': 'Report and archive integrity only. Final source/build verification was recorded before releasing that freeze; this check does not claim that later source edits belong to the same build.'}
    verification.write_text(json.dumps(result, indent=2) + '\n')
    files = {str(path.relative_to(BASE)): {'sha256': digest(path.read_bytes()),
             'bytes': path.stat().st_size} for path in sorted(BASE.rglob('*'))
             if path.is_file() and '__pycache__' not in path.parts and path != output}
    manifest = {'schema': 1,
        'scope': 'Frozen optional polynomial-closure-v1 source integration and lower-loop regression attempt, including all four failed three-loop controls. No three-loop or four-loop numerical acceptance; no supplied oracle numerical comparison.',
        'new_root_tests_passed': 105, 'native_tests_reused_by_exact_identity': 84,
        'new_native_test_invocations': 0, 'successful_lower_loop_runs': 3,
        'failed_three_loop_runs': 4, 'three_loop_predictions': 0,
        'archive_map_sha256': digest(archive_path.read_bytes()),
        'verification_sha256': digest(verification.read_bytes()),
        'files': files, 'file_count': len(files),
        'bytes': sum(row['bytes'] for row in files.values())}
    output.write_text(json.dumps(manifest, indent=2) + '\n')
    assert all(digest((BASE / name).read_bytes()) == row['sha256'] for name, row in files.items())
    print(json.dumps({'status': 'frozen', 'file_count': len(files),
                      'manifest_sha256': digest(output.read_bytes()),
                      'verified_archives': len(mapping)}))


if __name__ == '__main__':
    main()
