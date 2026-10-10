from pathlib import Path
import subprocess, time, json, hashlib, os
ROOT=Path('/common/dev/rustflow_fermi')
REPORT=ROOT/'reports/validation/2026-10-10-partial-endpoint-capability-prototype'
BINARY=Path('/tmp/rustflow-partial-endpoint-20261010/prototype-tests')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
build=json.loads((REPORT/'build-01.json').read_text())
assert sha(BINARY)==build['executable']['sha256']
evidence=REPORT/'final-evidence';evidence.mkdir(exist_ok=True)
origin=evidence/'origin';origin.mkdir(exist_ok=True)
endpoint=evidence/'endpoint';endpoint.mkdir(exist_ok=True)
command=[str(BINARY),'finite_density::partial_origin::','--test-threads=1','--nocapture']
start=time.monotonic();result=subprocess.run(command,capture_output=True,text=True,timeout=120,env={**os.environ,'PARTIAL_ORIGIN_REPORT':str(origin),'PARTIAL_ENDPOINT_REPORT':str(endpoint)})
log=REPORT/'final-gates.log';log.write_text(result.stdout+result.stderr)
report={'command':command,'exit_code':result.returncode,'wall_seconds':time.monotonic()-start,'binary_sha256':sha(BINARY),'build_file':'build-01.json','build_sha256':sha(REPORT/'build-01.json'),'runner_sha256':sha(Path(__file__)),'log_sha256':sha(log),'evidence':{str(p.relative_to(REPORT)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(evidence.rglob('*')) if p.is_file()},'production_admission':False,'native_mode_authority':False,'cargo_invoked':False}
(REPORT/'final-gates.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(result.stdout+result.stderr);raise SystemExit(result.returncode)
