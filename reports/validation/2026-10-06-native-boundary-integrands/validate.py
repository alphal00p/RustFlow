import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

root = Path.cwd()
raw = root / 'target/native-boundary-validation/final'
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
tracked = subprocess.check_output(['git', 'ls-files', '-z']).decode().split('\0')
paths = sorted(p for p in tracked if p and (
    p in ['Cargo.lock', 'Cargo.toml', 'build.rs'] or
    p.startswith(('src/', 'build_support/', 'tests/', 'examples/', 'reference/', 'fixtures/'))
))
def snapshot():
    return {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}
before = snapshot()
(raw/'source-before.json').write_text(json.dumps({
    'commit': commit, 'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'files': before}, indent=2)+'\n')
targets = ['native_boundary_integrand', 'integrand', 'recursive', 'ft', 'linear', 'linear_two_directions', 'end_to_end', 'tensor_vacuum', 'single_mass', 'complex_kinematics', 'partial_cut_boundaries', 'mixed_cut_flow', 'native_cut_interfaces', 'regions', 'supplied_flow', 'hepkit', 'hepkit_cuts', 'native_arity', 'native_family_cache', 'native_graphs', 'native_parametric', 'native_physical', 'native_symmetry', 'reduction_graphs', 'scoped_reductions', 'python_cut_interfaces', 'cli', 'gg_hg_amplitude', 'kinematics', 'kinematic_derivative', 'core', 'cuts', 'massless_phase_space', 'physical_family', 'multimass_physical']
tests = ['cargo', 'test', '--locked', '--release', '--features', 'python_stubgen', '--lib']
for target in targets:
    tests.extend(['--test', target])
tests.extend(['--', '--test-threads=4'])
commands = [
    ('fmt', ['cargo', 'fmt', '--check']),
    ('clippy', ['cargo', 'clippy', '--locked', '--release', '--all-targets',
                '--features', 'python_stubgen', '--', '-D', 'warnings']),
    ('tests', tests),
]
records = []
for name, command in commands:
    start = time.monotonic()
    with (raw/(name+'.log')).open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    records.append({'name':name, 'command':command, 'exit_code':result.returncode,
                    'elapsed_seconds':time.monotonic()-start})
    (raw/'processes.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records[-1]), flush=True)
    if result.returncode:
        break
after = snapshot()
(raw/'source-after.json').write_text(json.dumps({'commit':commit,
    'unchanged':before == after, 'files':after}, indent=2)+'\n')
assert before == after, 'build inputs changed during validation'
raise SystemExit(next((p['exit_code'] for p in records if p['exit_code']), 0))
