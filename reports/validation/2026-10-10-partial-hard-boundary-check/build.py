#!/usr/bin/env python3
"""Compile the report-only harness against an explicitly captured coherent build.

No Cargo invocation. This must be launched only after the root releases the
coherent native library for standalone validation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]


def digest(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coherent-build', type=Path, required=True)
    args = parser.parse_args()
    record = BASE / 'native-build-binding.json'
    assert not record.exists(), 'preserve previous native build attempts'
    build_path = args.coherent_build.resolve()
    build = json.loads(build_path.read_text())
    assert build['checks']['combined_release_build'] == 'passed'
    root_library = ROOT / 'target/release/deps/libsymbolica_amflow-bc386786e4706ace.rlib'
    native_library = ROOT / 'target/release/deps/librustred-392075d0c6995fd1.rlib'
    for library in [root_library, native_library]:
        assert digest(library) == build['files'][str(library.relative_to(ROOT))]['sha256']
    libraries = {
        'symbolica_amflow': root_library,
        'symbolica': ROOT / 'target/release/deps/libsymbolica-02ed301b2d4bd7f2.rlib',
        'serde_json': ROOT / 'target/release/deps/libserde_json-84255895f7dd0f90.rlib',
    }
    source = BASE / 'native.rs'
    output = Path('/tmp/rustflow-partial-hard-boundary-20261010')
    command = ['rustc', '--edition=2024', '--crate-name', 'partial_hard_boundary_gate',
               '-C', 'opt-level=0', '-C', 'debuginfo=0']
    for name, library in libraries.items():
        command += ['--extern', f'{name}={library}']
    command += ['-L', 'dependency=' + str(ROOT / 'target/release/deps')]
    for directory in sorted((ROOT / 'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):
        command += ['-L', 'native=' + str(directory)]
    command += [str(source), '-o', str(output)]
    checked = [source, build_path, native_library, *libraries.values(), Path(__file__)]
    hashes = {str(p): digest(p) for p in checked}
    binding = {'scope': 'Standalone ordinary hard-boundary validation, compiled only from the report harness and coherent existing libraries.',
               'command': command, 'inputs_sha256': hashes,
               'coherent_build': str(build_path), 'compile_identities': build['compile_identities'],
               'cargo_invoked': False, 'reference_files_read': 0}
    environment = dict(os.environ, CARGO_CRATE_NAME='partial_hard_boundary_gate')
    process = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True, text=True)
    (BASE / 'native-build.log').write_text(process.stdout + process.stderr)
    binding['exit_code'] = process.returncode
    assert all(digest(Path(name)) == sha for name, sha in hashes.items())
    if process.returncode == 0:
        binding['executable'] = str(output)
        binding['executable_sha256'] = digest(output)
    record.write_text(json.dumps(binding, indent=2) + '\n')
    print(process.stdout + process.stderr, end='')
    raise SystemExit(process.returncode)


if __name__ == '__main__':
    main()
