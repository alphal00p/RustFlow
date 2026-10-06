import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

root = Path.cwd()
raw = root / 'target/python-sheet-integration'
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
tracked = subprocess.check_output(['git', 'ls-files', '-z']).decode().split('\0')
paths = sorted(p for p in tracked if p and (
    p in ['Cargo.lock', 'Cargo.toml', 'build.rs'] or
    p.startswith(('src/', 'build_support/', 'tests/', 'examples/', 'reference/'))
))
def snapshot():
    return {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}
before = snapshot()
(raw/'source-before.json').write_text(json.dumps({
    'commit': commit, 'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'files': before}, indent=2)+'\n')
targets = ['algebraic_cache', 'algebraic_endpoints', 'algebraic_preparation',
    'canonical_algebraic', 'constrained_singular_endpoints', 'epsilon_algebraic',
    'frobenius', 'prepared_frobenius', 'physical_conditions', 'physical_transport',
    'prescribed_cache', 'supplied_singular_endpoints', 'transport_cache',
    'python_algebraic_transport', 'python_prescribed_transport',
    'python_endpoint_constraints', 'python_cut_interfaces']
tests = ['cargo', 'test', '--locked', '--release', '--features', 'python_stubgen']
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
