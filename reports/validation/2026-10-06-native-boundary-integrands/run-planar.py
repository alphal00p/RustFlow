import datetime, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

root = Path.cwd()
raw = root/'target/native-boundary-validation'
run = raw/'scientific-planar'
run.mkdir(exist_ok=False)
binary = run/'gg_hg_boundaries'
shutil.copy2('/common/dev/amflow/target/release/examples/gg_hg_boundaries', binary)
def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()
frozen = json.loads((raw/'final/source-before.json').read_text())
assert all(sha(root/p) == h for p,h in frozen['files'].items())
command = [str(binary), '--family', 'pl', '--digits', '20', '--guard-digits', '40',
    '--order', '80', '--workers', '16', '--cache-directory', str(run/'state'),
    '--cancel-file', str(run/'cancel')]
record = {'status':'running', 'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'tested_commit':frozen['commit'], 'command':command, 'binary_sha256':sha(binary),
    'cold_numerical_and_exact_cache':True, 'reference_role':'Comparison only after native boundary fit.',
    'wall_cap_seconds':1800, 'cancel_grace_seconds':60}
(run/'process.json').write_text(json.dumps(record,indent=2)+'\n')
with (run/'describe.json').open('w') as out, (run/'describe.stderr').open('w') as err:
    result = subprocess.run(command+['--describe'],stdout=out,stderr=err)
assert result.returncode == 0
started = time.monotonic()
with (run/'stdout.log').open('w') as out, (run/'stderr.log').open('w') as err:
    process = subprocess.Popen(command,stdout=out,stderr=err)
    record['pid'] = process.pid
    (run/'process.json').write_text(json.dumps(record,indent=2)+'\n')
    try:
        code = process.wait(timeout=1800)
    except subprocess.TimeoutExpired:
        record['cooperative_wall_cancel'] = True
        (run/'cancel').touch()
        try:
            code = process.wait(timeout=60)
        except subprocess.TimeoutExpired:
            process.terminate()
            try:
                code = process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.kill()
                code = process.wait()
record.update(status='completed' if code == 0 else 'failed_or_cancelled', exit_code=code,
    elapsed_seconds=time.monotonic()-started,
    finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    source_unchanged=all(sha(root/p) == h for p,h in frozen['files'].items()),
    binary_sha256_after=sha(binary))
(run/'process.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
raise SystemExit(code)
