#!/usr/bin/env python3
"""Freeze the isolated native prototype evidence; no new compilation or test run."""
import gzip
import hashlib
import json
from pathlib import Path
import re

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, record):
    assert not path.exists(), path
    path.write_text(json.dumps(record, indent=2) + '\n')


def main():
    assert not (HERE / 'artifact-manifest.json').exists()
    builds = {kind: json.loads((HERE / f'{kind}-build.json').read_bytes())
              for kind in ['library', 'tests']}
    for kind, build in builds.items():
        assert build['exit_code'] == 0
        assert build['production_sources_unchanged'] and build['isolated_sources_unchanged']
        assert digest(HERE / 'isolated-native.patch') == build['patch_sha256']
        assert digest(Path(build['output']['path'])) == build['output']['sha256']
        for name, expected in build['base_sources'].items():
            assert digest(ROOT / 'vendor/rustred/crates/rustred-core' / name) == expected
        for name, expected in build['isolated_sources'].items():
            assert digest(Path('/tmp/rustred-verified-program-reuse-20261010') / name) == expected
        for item in build['externs'].values():
            assert digest(Path(item['path'])) == item['sha256']
    assert builds['library']['base_sources'] == builds['tests']['base_sources']
    assert builds['library']['isolated_sources'] == builds['tests']['isolated_sources']
    resources = json.loads((HERE / 'native-tests-resources.json').read_bytes())
    assert resources['exit_code'] == 0
    assert builds['tests']['output']['path'] in resources['command']
    found = re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;.*finished in ([0-9.]+)s',
                       (HERE / 'native-tests-resources.log').read_text())
    assert found == [('91', '0', '1', '0.10')]
    inner_path = HERE / 'reuse-profile/artifact-manifest.json'
    inner = json.loads(inner_path.read_bytes())
    for name, row in inner['files'].items():
        path = inner_path.parent / name
        assert digest(path) == row['sha256'] and path.stat().st_size == row['bytes']
    archives = json.loads((inner_path.parent / 'archive-map.json').read_bytes())
    for row in archives['entries']:
        path = inner_path.parent / row['archive_path']
        assert digest(path) == row['archive_sha256']
        data = gzip.decompress(path.read_bytes())
        assert hashlib.sha256(data).hexdigest() == row['sha256'] and len(data) == row['bytes']
    privacy = json.loads((inner_path.parent / 'privacy-controls/result.json').read_bytes())
    assert privacy['result'] == 'PASS' and len(privacy['cases']) == 6
    profile = json.loads((inner_path.parent / 'summary.json').read_bytes())
    assert profile['status'] == 'PASS' and all(row['all_encoded_comparisons_exact'] for row in profile['profiles'])
    record = {
        'status': 'frozen isolated prototype passed',
        'scope': 'Native immutable verified reuse only. No production integration, new closure, numerical period or three-/four-loop acceptance.',
        'patch_sha256': digest(HERE / 'isolated-native.patch'),
        'builds': {kind: {'record_sha256': digest(HERE / f'{kind}-build.json'),
                         'wall_seconds': build['wall_seconds'], 'output': build['output']}
                   for kind, build in builds.items()},
        'native_guarded_tests': {'passed': 91, 'failed': 0, 'ignored': 1,
                                 'harness_seconds': 0.10, 'wall_seconds': resources['wall_seconds'],
                                 'peak_child_rss_kib': resources['peak_child_rss_kib'],
                                 'resources_sha256': digest(HERE / 'native-tests-resources.json'),
                                 'log_sha256': digest(HERE / 'native-tests-resources.log'),
                                 'launch_sha256': digest(HERE / 'native-tests-resources.provenance.json')},
        'external_privacy_controls': {'passed': 6, 'positive_compile': 1, 'expected_compile_failures': 5},
        'reuse_profile': {'manifest_sha256': digest(inner_path), 'files_verified': len(inner['files']),
                          'large_encoded_restorations_verified': len(archives['entries']),
                          'summary_sha256': digest(inner_path.parent / 'summary.json')},
        'design_audit_sha256': digest(HERE / 'design-audit.md'),
        'production_and_isolated_sources_verified_unchanged': True,
        'native_proof_schema': 'rustred.guarded-source-program.v2',
        'cold_decode_still_replays': True, 'production_admission_changed': False,
    }
    save(HERE / 'summary.json', record)
    files = {str(p.relative_to(HERE)): {'sha256': digest(p), 'bytes': p.stat().st_size}
             for p in sorted(HERE.rglob('*')) if p.is_file() and p != HERE / 'artifact-manifest.json'}
    save(HERE / 'artifact-manifest.json', {'status': 'frozen', 'scope': record['scope'],
                                         'files': files, 'file_count': len(files)})
    print(json.dumps({'status': 'frozen', 'files': len(files),
                      'manifest_sha256': digest(HERE / 'artifact-manifest.json')}))


if __name__ == '__main__':
    main()
