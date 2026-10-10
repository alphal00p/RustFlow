"""Link the joined validation harness to a captured, gated production library.

No Cargo build or prediction run. Sources and direct libraries are hash-bound;
the executable lives outside Git. Later runs use only the immutable capsule.
"""
from pathlib import Path
import datetime
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = HERE.parent / '2026-10-10-singleton-physical-zero-integration-v2'
PREFIX = 'singleton-physical-zero-v2'
CAP = Path('/tmp/rustflow-three-loop-joint-20261010')

def h(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')

def main():
    assert not CAP.exists() and not (HERE / 'capsule-binding.json').exists()
    spec = importlib.util.spec_from_file_location('frozen_sources', PARENT / 'frozen_sources.py')
    owner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(owner)
    sp = PARENT / (PREFIX + '-source-hashes.json')
    bp = PARENT / (PREFIX + '-build-provenance.json')
    gp = PARENT / (PREFIX + '-root-gates.json')
    snapshot = owner.verify_snapshot(sp)
    build = json.loads(bp.read_text())
    gates = json.loads(gp.read_text())
    assert build['source_snapshot_sha256'] == h(sp)
    assert build['checks']['combined_release_build'] == 'passed'
    assert gates['status'] == 'passed' and gates['build_provenance']['sha256'] == h(bp)
    assert gates['source_snapshot']['sha256'] == h(sp)
    assert all(row['status'] == 'passed' and row['failed'] == 0 for row in gates['gates'])
    harness = json.loads((HERE / 'harness-binding.json').read_text())
    assert harness['original_sha256'] == snapshot['files'][harness['original_path']]
    assert harness['generated_sha256'] == h(HERE / 'joint-runtime.rs')
    assert harness['generator_sha256'] == h(HERE / 'generate-harness.py')
    artifacts = [sp, bp, gp, PARENT / (PREFIX + '-actual62-source-equivalence.json'),
                 PARENT / 'run-resource-command.py']
    for gate in gates['gates']:
        for item in gate['artifacts'].values():
            path = ROOT / item['path']
            assert h(path) == item['sha256']
            artifacts.append(path)
    CAP.mkdir()
    (CAP / 'deps').mkdir()
    (HERE / 'evidence').mkdir(exist_ok=True)
    assert not list((HERE / 'evidence').iterdir()), 'preserve prior evidence capture'
    # Symlinks are used only during this source-frozen link. The resulting
    # statically linked executable, not these development libraries, runs later.
    for source in (ROOT / 'target/release/deps').iterdir():
        if source.suffix in ('.rlib', '.so'):
            (CAP / 'deps' / source.name).symlink_to(source)
    libraries = {}
    for name in ['symbolica_amflow', 'rustred']:
        matches = [p for p in build['files'] if Path(p).name.startswith('lib' + name + '-') and p.endswith('.rlib')]
        assert len(matches) == 1
        source = ROOT / matches[0]
        assert h(source) == build['files'][matches[0]]['sha256']
        libraries[name] = CAP / 'deps' / source.name
    root_suffix = libraries['symbolica_amflow'].stem.removeprefix('libsymbolica_amflow-')
    root_fingerprint = ROOT / 'target/release/.fingerprint' / ('symbolica-amflow-' + root_suffix) / 'lib-symbolica_amflow.json'
    fingerprint = json.loads(root_fingerprint.read_text())
    artifacts.append(root_fingerprint)
    dependency_selection = {}
    for name in ['symbolica', 'serde_json']:
        recorded = [row[3] for row in fingerprint['deps'] if row[1] == name]
        assert len(recorded) == 1
        matches = []
        for marker in (ROOT / 'target/release/.fingerprint').glob(name + '-*/lib-' + name):
            # Cargo stores this dependency fingerprint as eight little-endian
            # bytes in hexadecimal. Bind the exact direct dependency, not age.
            if int.from_bytes(bytes.fromhex(marker.read_text().strip()), 'little') == recorded[0]:
                suffix = marker.parent.name.removeprefix(name + '-')
                library = CAP / 'deps' / ('lib' + name + '-' + suffix + '.rlib')
                if library.is_file():
                    matches.append((marker, library))
        assert len(matches) == 1, (name, matches)
        marker, libraries[name] = matches[0]
        artifacts.extend([marker, marker.with_suffix('.json')])
        dependency_selection[name] = {'root_dependency_fingerprint': recorded[0],
            'fingerprint_marker': str(marker.relative_to(ROOT)), 'fingerprint_sha256': h(marker),
            'metadata_sha256': h(marker.with_suffix('.json'))}
    library_binding = {name: {'path': str(path), 'sha256': h(path), 'bytes': path.stat().st_size}
                       for name, path in libraries.items()}
    source = CAP / 'joint-runtime.rs'
    shutil.copy2(HERE / 'joint-runtime.rs', source)
    compiler = Path(shutil.which('rustc')).resolve()
    compiler_version = subprocess.run([str(compiler), '-vV'], text=True,
                                      capture_output=True, check=True).stdout
    compiler_binding = {'path': str(compiler), 'sha256': h(compiler),
                        'version': compiler_version}
    command = [str(compiler), '--edition=2024', '--crate-name', 'rustflow_joint_validation', '--test',
               '-C', 'opt-level=0', '-C', 'debuginfo=0']
    for name, path in libraries.items():
        command += ['--extern', name + '=' + str(path)]
    command += ['-L', 'dependency=' + str(CAP / 'deps')]
    for path in sorted((ROOT / 'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):
        command += ['-L', 'native=' + str(path)]
    python_lib = '/nix/store/0bq2kpdgy8kimy82pzpjjw787h1sx6a9-python3-3.14.6-env/lib'
    command += ['-L', 'native=' + python_lib, '-C', 'link-arg=-Wl,-rpath,' + python_lib,
                str(source), '-o', str(CAP / 'runtime')]
    started = time.monotonic()
    result = subprocess.run(command, text=True, capture_output=True,
                            env={**os.environ, 'CARGO_MANIFEST_DIR': str(ROOT),
                                 'CARGO_CRATE_NAME': 'rustflow_joint_validation'})
    (HERE / 'build.log').write_text(result.stdout + result.stderr)
    save(HERE / 'build.json', {'command': command, 'exit_code': result.returncode,
         'wall_seconds': time.monotonic() - started, 'libraries': library_binding,
         'compiler': compiler_binding,
         'direct_dependency_selection': dependency_selection,
         'root_library_fingerprint_sha256': h(root_fingerprint),
         'source_sha256': h(source), 'shared_target_modified': False, 'cargo_invoked': False,
         'production_library_build_provenance_sha256': h(bp)})
    print(result.stdout + result.stderr)
    if result.returncode:
        raise SystemExit(result.returncode)
    for row in library_binding.values():
        assert h(Path(row['path'])) == row['sha256']
    input_name = 'examples/finite_density/massless_three_loop_chain.json'
    assert h(ROOT / input_name) == snapshot['files'][input_name]
    shutil.copy2(ROOT / input_name, CAP / 'input.json')
    shutil.copy2(ROOT / input_name, HERE / 'input.json')
    shutil.copy2(PARENT / 'run-resource-command.py', CAP / 'run-resource-command.py')
    os.chmod(CAP / 'runtime', 0o555)
    copies = []
    for path in dict.fromkeys(artifacts):
        copy = HERE / 'evidence' / path.name
        shutil.copy2(path, copy)
        assert h(copy) == h(path)
        copies.append({'source_path': str(path.relative_to(ROOT)), 'copy_path': str(copy.relative_to(HERE)),
                       'sha256': h(copy), 'bytes': copy.stat().st_size})
    archive = HERE / 'frozen-source-assets.tar.gz'
    with archive.open('wb') as raw, gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as gz, tarfile.open(fileobj=gz, mode='w|') as tar:
        for name, expected in sorted(snapshot['files'].items()):
            path = ROOT / name
            assert h(path) == expected
            info = tar.gettarinfo(str(path), arcname=name)
            info.mtime = info.uid = info.gid = 0
            info.uname = info.gname = ''
            info.mode = 0o644
            with path.open('rb') as stream:
                tar.addfile(info, stream)
    with tarfile.open(archive, 'r:gz') as tar:
        restored = {member.name: hashlib.file_digest(tar.extractfile(member), 'sha256').hexdigest() for member in tar}
    assert restored == snapshot['files']
    owner.verify_snapshot(sp)
    ldd = subprocess.run(['ldd', str(CAP / 'runtime')], text=True, capture_output=True)
    assert ldd.returncode == 0
    (HERE / 'runtime-dynamic-libraries.txt').write_text(ldd.stdout + ldd.stderr)
    record = {
        'scope': 'Fresh original complete three-loop input; one production preparation shared by unchanged fixed-D and Laurent phase bodies. No imported rules/checkpoints or reference values.',
        'captured_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'capsule_path': str(CAP), 'source_snapshot_sha256': h(sp), 'build_provenance_sha256': h(bp),
        'root_gates_sha256': h(gp), 'compile_identities': gates['compile_identities'], 'features': build['features'],
        'source_asset_count': len(snapshot['files']),
        'source_archive': {'path': archive.name, 'sha256': h(archive), 'bytes': archive.stat().st_size,
                           'all_restored_member_hashes_match_snapshot': True},
        'evidence_copies': copies, 'harness_binding_sha256': h(HERE / 'harness-binding.json'),
        'standalone_build_sha256': h(HERE / 'build.json'),
        'static_files': [{'path': str(path), 'sha256': h(path), 'bytes': path.stat().st_size}
                         for path in sorted(CAP.iterdir()) if path.is_file()],
        'post_run_policy': 'Captured executable/input/harness/source archive/evidence only; no dependency on later shared source/build state.',
        'large_executable_storage': '/tmp only; do not commit executable', 'prepare_script_sha256': h(Path(__file__)),
    }
    save(HERE / 'capsule-binding.json', record)
    print(json.dumps({'runtime_sha256': h(CAP / 'runtime'), 'capsule_binding_sha256': h(HERE / 'capsule-binding.json')}))

if __name__ == '__main__':
    main()
