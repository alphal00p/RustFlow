#!/usr/bin/env python3
"""Losslessly archive selected completed runs, preserving explicit probe inputs."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
PREFIX = 'polynomial-closure-v2'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readers-released', action='store_true', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('runs', nargs='+')
    args = parser.parse_args()
    output = BASE / args.manifest
    assert not output.exists()
    protected = set()
    candidates = []
    for suffix in args.runs:
        assert '/' not in suffix and suffix.startswith(('two-loop-', 'massive-', 'three-loop-'))
        stem = PREFIX + '-' + suffix
        binding = json.loads((BASE / (stem + '-binding.json')).read_text())
        resources = json.loads((BASE / (stem + '-resources.json')).read_text())
        assert isinstance(binding['exit_code'], int)
        assert resources['exit_code'] == binding['exit_code']
        assert binding['post_run_source_and_executable_unchanged']
        directory = BASE / stem
        for path in (directory / 'native-closure').rglob('*'):
            if not path.is_file() or path.suffix == '.gz':
                continue
            # Successful final proofs remain directly available to comparison
            # bindings. Failed round checkpoints are reserved for peer probes.
            keep = '-closed.' in path.name
            if suffix.startswith('three-loop-fixed'):
                keep |= path.name in ('round-002-provisional.json', 'round-002-provisional.bin',
                                      'round-003-provisional.json', 'round-003-provisional.bin',
                                      'conditional-point-refinements.json')
            if keep:
                protected.add(path)
            else:
                candidates.append(path)
        log = BASE / (stem + '-resources.log')
        assert log.exists()
        candidates.append(log)
    assert len(candidates) == len(set(candidates))
    preserved = {str(path.relative_to(BASE)): {'sha256': sha(path.read_bytes()),
                  'bytes': path.stat().st_size} for path in sorted(protected)}
    entries = []
    for path in sorted(candidates):
        data = path.read_bytes()
        target = path.with_name(path.name + '.gz')
        assert not target.exists()
        compressed = gzip.compress(data, compresslevel=9, mtime=0)
        assert gzip.decompress(compressed) == data
        target.write_bytes(compressed)
        assert gzip.decompress(target.read_bytes()) == data
        entries.append({'original_path': str(path.relative_to(BASE)),
                        'archive_path': str(target.relative_to(BASE)),
                        'original_sha256': sha(data), 'archive_sha256': sha(compressed),
                        'original_bytes': len(data), 'archive_bytes': len(compressed)})
        path.unlink()
    assert all(sha((BASE / name).read_bytes()) == row['sha256']
               for name, row in preserved.items())
    report = {'schema': 1, 'scope': 'Completed verbose native payloads and logs only. Inputs, configuration, predictions, comparison reports, final closed proofs and reserved failure checkpoints are retained raw. No proof or numerical value is modified.',
              'runs': args.runs, 'gzip_mtime': 0, 'compression_level': 9,
              'reader_release_explicit': True, 'preserved_raw': preserved,
              'entries': entries, 'file_count': len(entries),
              'original_bytes': sum(row['original_bytes'] for row in entries),
              'archive_bytes': sum(row['archive_bytes'] for row in entries),
              'all_restore_checks_passed': True,
              'archive_script_sha256': sha(Path(__file__).read_bytes())}
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('file_count', 'original_bytes', 'archive_bytes')}))


if __name__ == '__main__':
    main()
