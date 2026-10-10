from pathlib import Path
import subprocess,time,json,hashlib
R=Path(__file__).resolve().parent
t=time.monotonic()
with (R/'tests.log').open('wb') as log:
 try:
  process=subprocess.run(['/tmp/rustflow-integrated-angular-aggregate-20261010','--test-threads=1'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
  status=process.returncode
 except subprocess.TimeoutExpired:status=124
record={'exit_code':status,'wall_seconds':time.monotonic()-t,'timeout_seconds':120,'source_sha256':hashlib.sha256((R/'angular.rs').read_bytes()).hexdigest()}
(R/'tests.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record));print((R/'tests.log').read_text())
raise SystemExit(status)
