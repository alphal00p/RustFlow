"""Prepare or compile one complete isolated RustFlow crate against the memo graph.

Preparation and compilation are separate explicit commands. No Cargo invocation,
live source edit, old RustFlow --extern, or shared target output is permitted.
The overlays must already contain the agreed module declarations and tests.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import tomllib

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = ROOT / 'reports/validation/2026-10-10-terminal-promotion-integration'
SNAPSHOT = BASE / 'terminal-promotion-source-hashes.json'
BUILD = BASE / 'terminal-promotion-build-provenance.json'
GATES = BASE / 'terminal-promotion-root-gates.json'
CAP = Path('/tmp/rustflow-partial-flow-full-crate-20261010')
COMPILER = Path('/nix/store/paaihvr59q1101200mabmk7y4yhw42g7-rustc-wrapper-1.97.1/bin/rustc')
LINKER = Path('/nix/store/w88q44gqd1qg5wmkk7v0h97rpiqvam0l-gcc-wrapper-15.3.0/bin/cc')
PYTHON_LIB = Path('/nix/store/0bq2kpdgy8kimy82pzpjjw787h1sx6a9-python3-3.14.6-env/lib')


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, data):
    assert not path.exists(), f'Preserve prior attempt: {path}'
    path.write_text(json.dumps(data, indent=2) + '\n')


def record(path):
    return {'path': str(path), 'sha256': digest(path), 'bytes': path.stat().st_size}


def inventory(directory):
    result = {}
    for path in sorted(directory.rglob('*')):
        assert not path.is_symlink(), f'Source symlink is not a frozen copy: {path}'
        if path.is_file():
            result[str(path.relative_to(directory))] = digest(path)
    return result


def dependency_graph(build):
    root_name = next(name for name in build['files']
                     if Path(name).name.startswith('libsymbolica_amflow-') and name.endswith('.rlib'))
    assert digest(ROOT / root_name) == build['files'][root_name]['sha256']
    suffix = Path(root_name).stem.removeprefix('libsymbolica_amflow-')
    fp = ROOT / 'target/release/.fingerprint'
    root_fp = fp / f'symbolica-amflow-{suffix}/lib-symbolica_amflow.json'
    graph = json.loads(root_fp.read_text())
    test_record = build['feature_fingerprints']['rustflow']
    test_fp = ROOT / test_record['path']
    assert digest(test_fp) == test_record['sha256']
    test_graph = json.loads(test_fp.read_text())
    assert graph['deps'] == test_graph['deps']
    assert graph['features'] == test_graph['features']
    selected = {}
    evidence = [record(root_fp), record(test_fp)]
    for _, name, _, expected in graph['deps']:
        if name == 'build_script_build':
            continue
        matches = []
        for marker in fp.glob('*/lib-' + name):
            raw = bytes.fromhex(marker.read_text().strip())
            if len(raw) != 8 or int.from_bytes(raw, 'little') != expected:
                continue
            suffix = marker.parent.name.rsplit('-', 1)[1]
            for extension in ['rlib', 'so']:
                library = ROOT / f'target/release/deps/lib{name}-{suffix}.{extension}'
                if library.is_file():
                    matches.append((marker, library))
        assert len(matches) == 1, (name, matches)
        marker, library = matches[0]
        selected[name] = record(library)
        selected[name]['cargo_dependency_fingerprint'] = expected
        evidence += [record(marker), record(marker.with_suffix('.json'))]
    assert 'symbolica_amflow' not in selected
    native_name = next(name for name in build['files']
                       if Path(name).name.startswith('librustred-') and name.endswith('.rlib'))
    assert selected['rustred']['sha256'] == build['files'][native_name]['sha256']
    return graph, selected, evidence


def prepare(overlays, revision):
    cap = CAP.with_name(CAP.name + '-' + revision)
    state = HERE / revision
    assert not cap.exists(), 'Use a fresh capsule; preserve failed preparations.'
    assert not state.exists(), 'Use a fresh revision; preserve prior preparations.'
    snapshot = json.loads(SNAPSHOT.read_text())
    build = json.loads(BUILD.read_text())
    gates = json.loads(GATES.read_text())
    assert build['source_snapshot_sha256'] == digest(SNAPSHOT)
    assert build['checks']['combined_release_build'] == 'passed'
    assert gates['status'] == 'passed'
    assert gates['build_provenance']['sha256'] == digest(BUILD)
    assert gates['source_snapshot']['sha256'] == digest(SNAPSHOT)
    assert all(row['failed'] == 0 and row['status'] == 'passed' for row in gates['gates'])
    graph, selected, fingerprints = dependency_graph(build)
    base_files = dict(snapshot['files'])
    # The build inventory does not include every doc/fixture embedded by lib.rs.
    for prefix in ['docs', 'fixtures']:
        for path in (ROOT / prefix).rglob('*'):
            if path.is_file():
                base_files.setdefault(str(path.relative_to(ROOT)), digest(path))
    for name, expected in base_files.items():
        assert digest(ROOT / name) == expected, f'Memo baseline changed: {name}'
    overlay_files = {}
    origins = []
    for directory in overlays:
        assert directory.is_dir()
        for path in sorted(directory.rglob('*')):
            if not path.is_file():
                continue
            assert not path.is_symlink()
            relative = str(path.relative_to(directory))
            assert relative.startswith(('src/', 'tests/')), relative
            value = digest(path)
            assert relative not in overlay_files or overlay_files[relative][1] == value, relative
            overlay_files[relative] = (path, value)
            origins.append({'destination': relative, **record(path)})
    assert 'src/finite_density/mod.rs' in overlay_files, 'Explicit module integration is required.'
    module_text = overlay_files['src/finite_density/mod.rs'][0].read_text()
    assert 'mod partial_origin;' in module_text and 'mod source_class;' in module_text
    state.mkdir()
    shutil.copyfile(Path(__file__), state / 'build-full-crate.py')
    source = cap / 'source'
    source.mkdir(parents=True)
    for name, expected in base_files.items():
        destination = source / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, destination)
        assert digest(destination) == expected
    for name, (original, expected) in overlay_files.items():
        destination = source / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, destination)
        assert digest(destination) == expected
    sources = inventory(source)
    source_hash = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    (cap / 'deps').mkdir()
    libraries = {}
    for path in sorted((ROOT / 'target/release/deps').iterdir()):
        if path.suffix not in ['.rlib', '.so'] or path.name.startswith('libsymbolica_amflow-'):
            continue
        (cap / 'deps' / path.name).symlink_to(path)
        libraries[str(path)] = digest(path)
    native_dirs = sorted(set(
        p for pattern in ['gmp-mpfr-sys-*/out/lib', 'blake3-*/out', 'psm-*/out']
        for p in (ROOT / 'target/release/build').glob(pattern)
    ))
    native_archives = {str(p): digest(p) for d in native_dirs for p in sorted(d.rglob('*.a'))}
    package = tomllib.loads((source / 'Cargo.toml').read_text())['package']
    identities = dict(build['compile_identities'])
    identities['PORT_SOURCE_DIGEST'] = 'isolated-partial-flow-v1:sha256:' + source_hash
    environment = {**identities, 'CARGO_MANIFEST_DIR': str(source),
                   'CARGO_PKG_NAME': package['name'], 'CARGO_PKG_VERSION': package['version'],
                   'CARGO_CRATE_NAME': 'symbolica_amflow'}
    compiler_version = subprocess.run([str(COMPILER), '-vV'], check=True,
                                      capture_output=True, text=True).stdout
    prepared = {
        'scope': 'One copied full crate; memo native/dependency graph reused; no old RustFlow library in externs.',
        'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'cargo_invoked': False, 'compile_executed': False, 'shared_target_written': False,
        'baseline': [record(SNAPSHOT), record(BUILD), record(GATES)],
        'baseline_files': base_files, 'overlays': origins, 'source_files': sources,
        'source_inventory_sha256': source_hash, 'source_root': str(source), 'capsule': str(cap),
        'features': json.loads(graph['features']), 'direct_libraries': selected,
        'dependency_fingerprints': fingerprints, 'library_files': libraries,
        'native_directories': [str(p) for p in native_dirs], 'native_archives': native_archives,
        'compiler': {**record(COMPILER), 'version': compiler_version}, 'linker': record(LINKER),
        'compile_environment': environment, 'builder_sha256': digest(Path(__file__)),
        'limits': 'Rustc invoked only by a separate build command. Native/CAS work bounds are not proved here.',
    }
    save(state / 'prepared.json', prepared)
    print(json.dumps({'prepared': str(state / 'prepared.json'), 'files': len(sources),
                      'direct_dependencies': len(selected), 'compilation_started': False}))


def verify(prepared):
    assert inventory(Path(prepared['source_root'])) == prepared['source_files']
    for key in ['library_files', 'native_archives']:
        for name, expected in prepared[key].items():
            assert digest(Path(name)) == expected, name
    for item in prepared['dependency_fingerprints'] + prepared['baseline']:
        assert digest(Path(item['path'])) == item['sha256'], item['path']
    for key in ['compiler', 'linker']:
        assert digest(Path(prepared[key]['path'])) == prepared[key]['sha256']
    # Absolute test fixture includes remain source-bound during compilation.
    for name, expected in prepared['baseline_files'].items():
        if name.startswith(('examples/', 'fixtures/')):
            assert digest(ROOT / name) == expected, name


def build(mode, attempt, timeout, revision):
    assert attempt and all(c.isalnum() or c in '-_' for c in attempt)
    state = HERE / revision
    prepared_path = state / 'prepared.json'
    prepared = json.loads(prepared_path.read_text())
    cap = Path(prepared['capsule'])
    assert prepared['builder_sha256'] == digest(Path(__file__)), 'Freeze this launcher before preparation.'
    verify(prepared)
    folder = state / attempt
    folder.mkdir()
    output = cap / (attempt + ('.rmeta' if mode.endswith('check') else '.bin'))
    command = [str(COMPILER), '--edition=2024', '--crate-name', 'symbolica_amflow',
               '-C', 'linker=' + str(LINKER), '-C', 'opt-level=0', '-C', 'debuginfo=0',
               '-C', 'metadata=isolated_partial_flow_' + attempt]
    command += ['--test'] if mode in ['test', 'test-check'] else ['--crate-type', 'rlib']
    if mode.endswith('check'):
        command += ['--emit=metadata']
    for feature in prepared['features']:
        command += ['--cfg', 'feature="' + feature + '"']
    for name, row in prepared['direct_libraries'].items():
        assert name != 'symbolica_amflow'
        command += ['--extern', name + '=' + str(cap / 'deps' / Path(row['path']).name)]
    command += ['-L', 'dependency=' + str(cap / 'deps')]
    for path in prepared['native_directories'] + [str(PYTHON_LIB)]:
        command += ['-L', 'native=' + path]
    command += ['-C', 'link-arg=-Wl,-rpath,' + str(PYTHON_LIB),
                str(Path(prepared['source_root']) / 'src/lib.rs'), '-o', str(output)]
    save(folder / 'launch.json', {'command': command, 'mode': mode,
         'prepared_sha256': digest(prepared_path), 'timeout_seconds': timeout,
         'environment_overrides': prepared['compile_environment']})
    start = time.monotonic()
    with (folder / 'build.log').open('w') as log:
        try:
            result = subprocess.run(command, cwd=prepared['source_root'], stdout=log, stderr=subprocess.STDOUT,
                                    env={**os.environ, **prepared['compile_environment']}, timeout=timeout)
            exit_code = result.returncode
        except subprocess.TimeoutExpired:
            exit_code = 124
    after = None
    try:
        verify(prepared)
    except Exception as error:
        after = repr(error)
    save(folder / 'result.json', {'exit_code': exit_code, 'wall_seconds': time.monotonic() - start,
         'mode': mode, 'compile_only': True, 'tests_executed': False,
         'cargo_invoked': False, 'shared_target_written': False,
         'post_integrity_passed': after is None, 'post_integrity_error': after,
         'prepared_sha256': digest(prepared_path), 'launcher_sha256': digest(Path(__file__)),
         'output': record(output) if output.is_file() else None})
    raise SystemExit(exit_code if exit_code else (1 if after else 0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    preparation = commands.add_parser('prepare')
    preparation.add_argument('--revision', required=True)
    preparation.add_argument('--overlay', type=Path, action='append', required=True)
    compilation = commands.add_parser('build')
    compilation.add_argument('--revision', required=True)
    compilation.add_argument('--mode', choices=['check', 'test-check', 'test'], required=True)
    compilation.add_argument('--attempt', required=True)
    compilation.add_argument('--timeout', type=int, default=1800)
    args = parser.parse_args()
    assert args.revision and all(c.isalnum() or c in '-_' for c in args.revision)
    if args.command == 'prepare':
        prepare([path.resolve() for path in args.overlay], args.revision)
    else:
        assert args.timeout > 0
        build(args.mode, args.attempt, args.timeout, args.revision)
