#!/usr/bin/env python3
"""Build a validation-only adapter after coherent capture; never invokes Cargo."""
import json
import subprocess
import time
from pathlib import Path
from frozen_sources import BASE, ROOT, digest, verify_snapshot

PREFIX = 'requested-ray-point'
OUTPUT = Path('/tmp/rustflow-requested-ray-point-digest')


def main():
    report_path = BASE / 'native-digest-build.json'
    assert not report_path.exists(), 'preserve the existing utility build'
    snapshot = BASE / (PREFIX + '-source-hashes.json')
    build_path = BASE / (PREFIX + '-build-provenance.json')
    build = json.loads(build_path.read_text())
    verify_snapshot(snapshot)
    assert build['source_snapshot_sha256'] == digest(snapshot)
    libraries = sorted((ROOT / 'target/release/deps').glob('libblake3-*.rlib'))
    assert libraries, 'use the already cached native dependency'
    library = libraries[0]
    source = BASE / 'native-digest.rs'
    inputs = {str(path): digest(path) for path in (source, library)}
    native_paths = sorted((ROOT / 'target/release/build').glob('blake3-*/out'))
    native_inputs = {str(path): digest(path) for directory in native_paths
                     for path in directory.glob('*.a')}
    command = ['nix', 'develop', '--command', 'rustc', '--edition=2024', '-O',
               '-L', 'dependency=target/release/deps',
               '--extern', 'blake3=' + str(library)]
    for directory in native_paths:
        command += ['-L', 'native=' + str(directory)]
    command += [str(source), '-o', str(OUTPUT)]
    started = time.monotonic()
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    (BASE / 'native-digest-build.log').write_text(result.stdout + result.stderr)
    assert all(digest(Path(path)) == sha for path, sha in (inputs | native_inputs).items())
    verify_snapshot(snapshot)
    report = {'scope': 'Tiny validation-only adapter linked to the cached native BLAKE3 dependency. No Cargo, RustFlow source mutation, native proof change or reference-value calculation.',
              'exit_code': result.returncode, 'wall_seconds': time.monotonic() - started,
              'command': command, 'source_snapshot_sha256': digest(snapshot),
              'build_provenance_sha256': digest(build_path), 'inputs_sha256': inputs,
              'native_archives_sha256': native_inputs,
              'builder_sha256': digest(Path(__file__))}
    if result.returncode == 0:
        report['executable'] = {'path': str(OUTPUT), 'sha256': digest(OUTPUT)}
        report['controls'] = []
        # Standard algorithm check vectors exercise both an empty and short file.
        for name, content, expected in [
            ('empty', b'', 'af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262'),
            ('abc', b'abc', '6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85')]:
            fixture = BASE / ('native-digest-control-' + name + '.bin')
            fixture.write_bytes(content)
            control = subprocess.run([str(OUTPUT), str(fixture)], text=True, capture_output=True)
            actual = control.stdout.strip()
            passed = control.returncode == 0 and actual == expected + ' ' + str(len(content))
            report['controls'].append({'name': name, 'bytes': len(content),
                                       'blake3': expected, 'actual_output': actual,
                                       'exit_code': control.returncode, 'stderr': control.stderr,
                                       'passed': passed})
        if not all(row['passed'] for row in report['controls']):
            report['exit_code'] = 1
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    raise SystemExit(report['exit_code'])


if __name__ == '__main__':
    main()
