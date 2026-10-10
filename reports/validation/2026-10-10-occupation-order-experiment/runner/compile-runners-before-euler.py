"""Build only standalone probes, against an explicitly supplied native library.
Run inside the existing Nix development shell; no Cargo invocation or writes.
"""
import argparse, datetime, hashlib, json, os, pathlib, subprocess, time
parser = argparse.ArgumentParser()
parser.add_argument('--tag', required=True, choices=['baseline', 'degree-first'])
parser.add_argument('--library', required=True)
args = parser.parse_args()
root = pathlib.Path.cwd()
report = root / 'reports/validation/2026-10-10-occupation-order-experiment'
source_dir = report / 'runner'
outputs = pathlib.Path('/tmp/rustflow-occupation-order-runners')
outputs.mkdir(exist_ok=True)
def meta(path):
    path = pathlib.Path(path)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
library = pathlib.Path(args.library).resolve()
if args.tag == 'baseline':
    assert meta(library)['sha256'] == 'ac85472324303a84e909023a9aac6fe77f4993934985519fd88275cae0e1aad0'
dependencies = {
    'rustred': library,
    'symbolica': root / 'target/release/deps/libsymbolica-02ed301b2d4bd7f2.rlib',
    'serde_json': root / 'target/release/deps/libserde_json-84255895f7dd0f90.rlib',
    'bincode': root / 'target/release/deps/libbincode-30fae4fd0bc9d1b5.rlib',
}
results = []
for kind in ['point-witness', 'active-pilot']:
    source = source_dir / (kind + '.rs')
    executable = outputs / (args.tag + '-' + kind)
    command = ['rustc', '--edition=2024', '--crate-name', kind.replace('-', '_'),
               '-C', 'opt-level=0', '-C', 'debuginfo=0']
    for name, path in dependencies.items():
        command += ['--extern', name + '=' + str(path)]
    command += ['-L', 'dependency=' + str(root / 'target/release/deps')]
    for directory in ['gmp-mpfr-sys-2e6d668657be102d', 'gmp-mpfr-sys-c789de89ed67d3db']:
        command += ['-L', 'native=' + str(root / 'target/release/build' / directory / 'out/lib')]
    command += [str(source), '-o', str(executable)]
    binding = {'captured_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'tag': args.tag, 'source': meta(source),
               'dependencies': {name: meta(path) for name, path in dependencies.items()},
               'command': command, 'cargo_invoked': False}
    binding_file = report / (args.tag + '-' + kind + '-build-binding.json')
    binding_file.write_text(json.dumps(binding, indent=2) + '\n')
    started = time.monotonic()
    with (report / (args.tag + '-' + kind + '-build.log')).open('w') as log:
        result = subprocess.run(command, env={**os.environ, 'CARGO_CRATE_NAME': kind.replace('-', '_')},
                                stdout=log, stderr=subprocess.STDOUT)
    item = {'kind': kind, 'returncode': result.returncode,
            'wall_seconds': time.monotonic() - started, 'command': command}
    if result.returncode == 0:
        item['executable'] = meta(executable)
    results.append(item)
    (report / (args.tag + '-runner-builds.json')).write_text(json.dumps(results, indent=2) + '\n')
    if result.returncode:
        raise SystemExit(result.returncode)
    print(args.tag, kind, 'compiled', item['wall_seconds'], flush=True)
