#!/usr/bin/env python3
"""Verify lossless archives and immutable evidence, then seal this attempt."""
import gzip
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIX = 'requested-ray-point'
read = lambda path: json.loads(path.read_text())
digest = lambda data: hashlib.sha256(data).hexdigest()


def main():
    output = BASE / 'artifact-manifest.json'
    verification = BASE / 'archive-and-binding-verification.json'
    assert not output.exists() and not verification.exists()
    archive_path = BASE / (PREFIX + '-native-archives.json')
    archives = read(archive_path)
    mapping = {row['original_path']: row for row in archives['entries']}
    assert len(mapping) == len(archives['entries'])
    for row in archives['entries']:
        compressed = (BASE / row['archive_path']).read_bytes()
        assert digest(compressed) == row['archive_sha256'] and len(compressed) == row['archive_bytes']
        raw = gzip.decompress(compressed)
        assert digest(raw) == row['original_sha256'] and len(raw) == row['original_bytes']
        assert not (BASE / row['original_path']).exists()
    for name, row in archives['preserved_raw'].items():
        data = (BASE / name).read_bytes()
        assert digest(data) == row['sha256'] and len(data) == row['bytes']

    def resolve(path):
        if path.exists():
            return path.read_bytes()
        name = str(path.relative_to(BASE))
        return gzip.decompress((BASE / mapping[name]['archive_path']).read_bytes())

    checked = 0
    summary = read(BASE / (PREFIX + '-checkpoint-summary.json'))
    for name, expected in summary['bindings_sha256'].items():
        assert digest(resolve(ROOT / name)) == expected
        checked += 1
    lower = read(BASE / (PREFIX + '-lower-loop-regressions.json'))
    for row in lower['gates']:
        for info in row['bindings'].values():
            assert digest(resolve(ROOT / info['path'])) == info['sha256']
            checked += 1
    for path in BASE.glob('*-comparison-binding.json'):
        data = read(path)
        assert data['status'] == 'passed'
        for name, expected in data['inputs_sha256'].items():
            assert digest(resolve(ROOT / name)) == expected
            checked += 1
    gates = read(BASE / (PREFIX + '-root-gates.json'))
    for gate in gates['gates']:
        for info in gate['artifacts'].values():
            assert digest(resolve(ROOT / info['path'])) == info['sha256']
            checked += 1
    for name, expected in read(BASE / 'comparison-plan.json')['frozen_comparator_and_reference_artifacts_sha256'].items():
        assert digest((ROOT / name).read_bytes()) == expected
        checked += 1

    nested_manifests = {}
    paths = [BASE / 'pilot-transfer-read-only/manifest.json',
             ROOT / 'reports/validation/2026-10-10-native-ray-point-fallback-fixture/manifest.json']
    for path in paths:
        data = read(path)
        for row in data['files']:
            raw = resolve(ROOT / row['path'])
            assert digest(raw) == row['sha256'] and len(raw) == row['bytes']
        nested_manifests[str(path.relative_to(ROOT))] = {
            'sha256': digest(path.read_bytes()), 'verified_files': len(data['files'])}
    # Reporting prose may be updated after the frozen build. Its separate hash
    # is not presented as a changed compilation input.
    docs = ROOT / 'docs/finite-density-plan.md'
    external = {str(docs.relative_to(ROOT)): {'sha256': digest(docs.read_bytes()),
        'scope': 'Post-run reporting text; separate from the frozen source/build identity.'}}
    result = {'status': 'passed', 'lossless_archives_verified': len(mapping),
        'preserved_raw_files_verified': len(archives['preserved_raw']),
        'historical_report_bindings_verified': checked, 'nested_manifests': nested_manifests,
        'external_reporting_bindings': external,
        'scope': 'Archive/report integrity only. The final source/build verification was captured before releasing that freeze. No later source edit is attributed to the recorded build.'}
    verification.write_text(json.dumps(result, indent=2) + '\n')
    files = {str(path.relative_to(BASE)): {'sha256': digest(path.read_bytes()), 'bytes': path.stat().st_size}
        for path in sorted(BASE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts and path != output}
    manifest = {'schema': 1,
        'scope': 'Optional requested-ray-point native discovery integration and lower-loop default-branch regression evidence, with all four failed three-loop numerical requests preserved. One cut-0 native connection closed; no three-loop or four-loop amplitude is accepted.',
        'new_root_tests_passed': 115, 'native_tests_reused_by_exact_identity': 84,
        'new_native_test_invocations': 0, 'successful_lower_loop_runs': 3,
        'failed_three_loop_runs': 4, 'three_loop_predictions': 0,
        'supplied_oracle_numerical_records_compared': 0,
        'archive_map_sha256': digest(archive_path.read_bytes()),
        'verification_sha256': digest(verification.read_bytes()),
        'files': files, 'file_count': len(files), 'bytes': sum(row['bytes'] for row in files.values())}
    output.write_text(json.dumps(manifest, indent=2) + '\n')
    assert all(digest((BASE / name).read_bytes()) == row['sha256'] for name, row in files.items())
    print(json.dumps({'status': 'frozen', 'file_count': len(files),
        'manifest_sha256': digest(output.read_bytes()), 'verified_archives': len(mapping),
        'verified_historical_bindings': checked}))

if __name__ == '__main__':
    main()
