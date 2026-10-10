from pathlib import Path
import subprocess,time,json,hashlib
R=Path(__file__).resolve().parent
cmd=['/tmp/rustflow-integrated-gram-shift-direct-20261010',str(R/'direct-result.json')]
t=time.monotonic()
record={'command':cmd,'timeout_seconds':180,'scope':'Generic bridge verification, not physical coefficient generation','check_sha256':hashlib.sha256((R/'check-direct.rs').read_bytes()).hexdigest()}
with (R/'direct-run.log').open('wb') as log:
 try:
  p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180)
  record['exit_code']=p.returncode
 except subprocess.TimeoutExpired: record['exit_code']=124
record['wall_seconds']=time.monotonic()-t
if (R/'direct-result.json').exists(): record['result_sha256']=hashlib.sha256((R/'direct-result.json').read_bytes()).hexdigest()
(R/'direct-run.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
raise SystemExit(record['exit_code'])
