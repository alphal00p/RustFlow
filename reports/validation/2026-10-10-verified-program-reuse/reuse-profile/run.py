#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
PRIOR = ROOT / 'reports/validation/2026-10-10-native-program-work-profile'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    mode = sys.argv[1]
    assert mode in ('single3', 'double')
    destination = HERE / mode
    assert not destination.exists(), 'preserve prior attempts'
    destination.mkdir()
    build = json.loads((HERE / 'build-binding.json').read_bytes())
    assert build['exit_code'] == 0
    native_file = Path(build['native_build']['path'])
    assert digest(native_file) == build['native_build']['sha256']
    native = json.loads(native_file.read_bytes())
    for name, expected in native['base_sources'].items():
        assert digest(ROOT / 'vendor/rustred/crates/rustred-core' / name) == expected
    libraries = build['dependencies']
    assert all(digest(Path(v['path'])) == v['sha256'] for v in libraries.values())
    exe = Path(build['executable']['path'])
    assert digest(exe) == build['executable']['sha256']
    corpus = PRIOR / mode / 'input.bin'
    old_manifest = json.loads((PRIOR / 'manifest.json').read_bytes())
    for row in old_manifest['files']:
        path = PRIOR / row['path']
        assert digest(path) == row['sha256']
        assert path.stat().st_size == row['bytes']
    inputs = [corpus, PRIOR / mode / 'input.json', PRIOR / 'manifest.json',
              PRIOR / 'input-bindings.json', HERE / 'probe.rs', HERE / 'build.py',
              HERE / 'build-binding.json', Path(__file__).resolve(), exe, native_file]
    hashes = {str(path): digest(path) for path in inputs}
    command = ['timeout', '600', str(exe), str(corpus), str(destination / 'profile')]
    record = {'scope': 'Same decoded real program, sequential verified then replayed empty union. One sample per phase; no fair end-to-end benchmark or closure claim.',
              'command': command, 'inputs': hashes, 'libraries': libraries,
              'prior_input_mode': mode, 'production_sources_verified_before': True,
              'frozen_input_report_manifest_verified': len(old_manifest['files'])}
    binding = destination / 'run-binding.json'
    binding.write_text(json.dumps(record, indent=2) + '\n')
    run = subprocess.run([sys.executable, str(ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),
                          str(destination / 'resources.json'), *command])
    record['exit_code'] = run.returncode
    record['inputs_unchanged'] = all(digest(Path(path)) == expected for path, expected in hashes.items())
    record['libraries_unchanged'] = all(digest(Path(v['path'])) == v['sha256'] for v in libraries.values())
    record['production_sources_unchanged'] = all(digest(ROOT / 'vendor/rustred/crates/rustred-core' / name) == expected for name, expected in native['base_sources'].items())
    binding.write_text(json.dumps(record, indent=2) + '\n')
    assert record['inputs_unchanged'] and record['libraries_unchanged'] and record['production_sources_unchanged']
    raise SystemExit(run.returncode)


if __name__ == '__main__':
    main()
