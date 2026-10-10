#!/usr/bin/env python3
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
NATIVE_BUILD = HERE.parent / 'library-build.json'
EXE = Path('/tmp/rustflow-verified-reuse-profile-20261010')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert not (HERE / 'build-binding.json').exists()
    native = json.loads(NATIVE_BUILD.read_bytes())
    assert native['exit_code'] == 0
    assert native['production_sources_unchanged'] and native['isolated_sources_unchanged']
    old = json.loads((ROOT / 'reports/validation/2026-10-10-native-program-work-profile/build-binding.json').read_bytes())
    dependencies = old['dependencies'].copy()
    dependencies['rustred'] = native['output']
    for item in dependencies.values():
        assert digest(Path(item['path'])) == item['sha256']
    source = HERE / 'probe.rs'
    command = ['rustc', '--edition=2024', '--crate-name', 'verified_reuse_profile',
               '-C', 'opt-level=0', '-C', 'debuginfo=0']
    for name, item in dependencies.items():
        command += ['--extern', f"{name}={item['path']}"]
    command += ['-L', f'dependency={ROOT}/target/release/deps']
    for path in sorted((ROOT / 'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):
        command += ['-L', f'native={path}']
    command += [str(source), '-o', str(EXE)]
    source_hash = digest(source)
    started = time.monotonic()
    run = subprocess.run(command, capture_output=True, text=True,
                         env={**os.environ, 'CARGO_CRATE_NAME': 'verified_reuse_profile'})
    (HERE / 'build.log').write_text(run.stdout + run.stderr)
    record = {'scope': 'Standalone diagnostic; no Cargo or shared-target outputs.',
              'native_build': {'path': str(NATIVE_BUILD), 'sha256': digest(NATIVE_BUILD)},
              'dependencies': dependencies, 'command': command,
              'source': {'path': str(source), 'sha256': source_hash},
              'compile_environment_override': {'CARGO_CRATE_NAME': 'verified_reuse_profile'},
              'exit_code': run.returncode, 'wall_seconds': time.monotonic() - started,
              'cargo_invoked': False}
    assert digest(source) == source_hash
    assert all(digest(Path(item['path'])) == item['sha256'] for item in dependencies.values())
    if run.returncode == 0:
        record['executable'] = {'path': str(EXE), 'sha256': digest(EXE), 'bytes': EXE.stat().st_size}
    (HERE / 'build-binding.json').write_text(json.dumps(record, indent=2) + '\n')
    print(run.stdout + run.stderr)
    print(json.dumps({'exit_code': run.returncode, 'wall_seconds': record['wall_seconds']}))
    raise SystemExit(run.returncode)


if __name__ == '__main__':
    main()
