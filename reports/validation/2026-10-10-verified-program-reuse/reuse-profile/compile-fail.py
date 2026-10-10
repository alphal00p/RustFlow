#!/usr/bin/env python3
"""External-crate privacy controls, including a successful immutable-access control."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    out = HERE / 'privacy-controls'
    assert not out.exists(), 'preserve prior attempts'
    out.mkdir()
    native_path = HERE.parent / 'library-build.json'
    native = json.loads(native_path.read_bytes())
    assert native['exit_code'] == 0
    library = Path(native['output']['path'])
    assert digest(library) == native['output']['sha256']
    prefix = 'use rustred::solver::guarded::{GuardedProgram,GuardedSourceSystem};\nuse std::sync::Arc;\n'
    cases = {
        'immutable_access': ('fn inspect(p: &GuardedProgram<1>) { let _ = (p.sources(),p.rules(),p.terminals()); }', None),
        'struct_literal': ('fn forge(s: Arc<GuardedSourceSystem<1>>) -> GuardedProgram<1> { GuardedProgram { sources:s, rules:vec![], terminals:Default::default() } }', 'E0451'),
        'private_field': ('fn mutate(p: &mut GuardedProgram<1>) { p.rules.clear(); }', 'E0616'),
        'rules_getter': ('fn mutate(p: &mut GuardedProgram<1>) { p.rules().reverse(); }', 'E0596'),
        'terminals_getter': ('fn mutate(p: &mut GuardedProgram<1>) { p.terminals().clear(); }', 'E0596'),
        'sources_getter': ('fn mutate(p: &mut GuardedProgram<1>) { let _ = Arc::get_mut(p.sources()); }', 'E0308'),
    }
    results = []
    for name, (body, expected) in cases.items():
        source = out / f'{name}.rs'
        source.write_text(prefix + body + '\n')
        command = ['rustc', '--edition=2024', '--crate-name', f'verified_privacy_{name}',
                   '--crate-type', 'lib', '--emit=metadata', '--error-format=json',
                   '--extern', f'rustred={library}', '-L', f'dependency={ROOT}/target/release/deps',
                   str(source), '-o', f'/tmp/verified-privacy-{name}.rmeta']
        started = time.monotonic()
        run = subprocess.run(command, capture_output=True, text=True,
                             env={**os.environ, 'CARGO_CRATE_NAME': f'verified_privacy_{name}'})
        (out / f'{name}.log').write_text(run.stdout + run.stderr)
        diagnostics = [json.loads(line) for line in run.stderr.splitlines() if line.startswith('{')]
        errors = [d for d in diagnostics if d.get('level') == 'error']
        codes = sorted({d['code']['code'] for d in errors if d.get('code')})
        if expected is None:
            assert run.returncode == 0 and not errors, (name, codes)
        else:
            assert run.returncode != 0 and expected in codes, (name, codes)
            assert set(codes) == {expected}, (name, 'unexpected independent compiler failure', codes)
        results.append({'case': name, 'command': command, 'source_sha256': digest(source),
                        'log_sha256': digest(out / f'{name}.log'), 'exit_code': run.returncode,
                        'expected_error': expected, 'error_codes': codes,
                        'wall_seconds': time.monotonic() - started, 'result': 'PASS'})
    assert digest(library) == native['output']['sha256']
    (out / 'result.json').write_text(json.dumps({'scope': 'External Rust privacy/type checks. Expected failures are successful controls; no runtime proof or numerical claim.',
        'native_build': {'path': str(native_path), 'sha256': digest(native_path)},
        'library': native['output'], 'result': 'PASS', 'cases': results}, indent=2) + '\n')
    print(json.dumps({'result': 'PASS', 'cases': len(results)}))


if __name__ == '__main__':
    main()
