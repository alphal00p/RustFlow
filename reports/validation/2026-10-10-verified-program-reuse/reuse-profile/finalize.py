#!/usr/bin/env python3
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    assert not path.exists(), path
    path.write_text(json.dumps(value, indent=2) + '\n')


def main():
    summaries = []
    archived = []
    build = json.loads((HERE / 'build-binding.json').read_bytes())
    assert build['exit_code'] == 0
    assert digest(Path(build['executable']['path'])) == build['executable']['sha256']
    native_path = Path(build['native_build']['path'])
    assert digest(native_path) == build['native_build']['sha256']
    native = json.loads(native_path.read_bytes())
    assert all(digest(ROOT / 'vendor/rustred/crates/rustred-core' / name) == expected
               for name, expected in native['base_sources'].items())
    for mode in ['single3', 'double']:
        directory = HERE / mode
        binding = json.loads((directory / 'run-binding.json').read_bytes())
        resources = json.loads((directory / 'resources.json').read_bytes())
        result = json.loads((directory / 'profile/result.json').read_bytes())
        assert resources['exit_code'] == binding['exit_code'] == 0
        assert binding['inputs_unchanged'] and binding['libraries_unchanged'] and binding['production_sources_unchanged']
        assert all(digest(Path(path)) == expected for path, expected in binding['inputs'].items())
        assert result['initial_loaded_proof_replay']
        assert result['verified_union_equals_loaded_encoded_bytes'] and result['verified_union_equals_replayed_encoded_bytes']
        phases = {row['phase']: row['wall_seconds'] for row in result['timings']}
        proof = directory / 'profile/unchanged-program.bin'
        raw = proof.read_bytes()
        archive = proof.with_suffix('.bin.gz')
        assert not archive.exists()
        packed = gzip.compress(raw, compresslevel=9, mtime=0)
        archive.write_bytes(packed)
        assert gzip.decompress(archive.read_bytes()) == raw
        raw_sha = hashlib.sha256(raw).hexdigest()
        archived.append({'logical_path': str(proof.relative_to(HERE)), 'sha256': raw_sha,
                         'bytes': len(raw), 'archive_path': str(archive.relative_to(HERE)),
                         'archive_sha256': digest(archive), 'archive_bytes': len(packed), 'gzip_mtime': 0})
        proof.unlink()
        summaries.append({'mode': mode, 'result': 'PASS', 'rules': result['rule_count'],
                          'sources': result['source_count'], 'phases_wall_seconds': phases,
                          'resource_wall_seconds': resources['wall_seconds'],
                          'peak_child_rss_kib': resources['peak_child_rss_kib'],
                          'verified_encoded_program_sha256': raw_sha,
                          'matches_original_input_bytes': raw_sha == digest(Path(next(p for p in binding['inputs'] if p.endswith('/input.bin')))),
                          'all_encoded_comparisons_exact': True,
                          'initial_decode_and_later_union_independently_replay': True,
                          'result_sha256': digest(directory / 'profile/result.json'),
                          'binding_sha256': digest(directory / 'run-binding.json')})
    privacy = json.loads((HERE / 'privacy-controls/result.json').read_bytes())
    assert privacy['result'] == 'PASS' and len(privacy['cases']) == 6
    save(HERE / 'archive-map.json', {'scope': 'Lossless storage for unchanged encoded output programs only. Original frozen input corpora stay raw.', 'entries': archived})
    save(HERE / 'summary.json', {
        'status': 'PASS', 'scope': 'Isolated prototype mechanism and privacy validation only; no production integration, closure, numerical result or end-to-end speed claim.',
        'source_build': {'path': str(native_path), 'sha256': digest(native_path)},
        'harness_build': {'path': str(HERE / 'build-binding.json'), 'sha256': digest(HERE / 'build-binding.json')},
        'profiles': summaries, 'external_privacy_controls_passed': 6,
        'privacy_result_sha256': digest(HERE / 'privacy-controls/result.json'),
        'production_sources_unchanged': True,
        'timing_caveat': 'Single operation samples on a shared host. Tiny verified timings are clock-sensitive. Encoding, construction and loading are separately recorded; whole closure/application cost was not measured.',
        'replay_scope': 'Large programs were loaded with full native replay and then independently replayed once more by union_replayed. Verified union was byte-identical between them, so a third large decode was not run. Root-owned small unit tests separately decode verified output programs.'
    })
    files = {str(p.relative_to(HERE)): {'sha256': digest(p), 'bytes': p.stat().st_size}
             for p in sorted(HERE.rglob('*')) if p.is_file() and p.name != 'artifact-manifest.json'}
    save(HERE / 'artifact-manifest.json', {'status': 'frozen', 'scope': 'Isolated prototype profile and external privacy evidence.', 'files': files, 'file_count': len(files)})
    print(json.dumps({'status': 'frozen', 'files': len(files), 'manifest_sha256': digest(HERE / 'artifact-manifest.json')}))


if __name__ == '__main__':
    main()
