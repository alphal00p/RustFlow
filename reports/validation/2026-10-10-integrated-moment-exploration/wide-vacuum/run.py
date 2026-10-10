from pathlib import Path
import subprocess,time,json,hashlib
R=Path(__file__).resolve().parent
cmd=['/tmp/rustflow-integrated-wide-vacuum-20261010',str(R/'result.json')]
t=time.monotonic()
record={'command':cmd,'timeout_seconds':180,'scope':'Native wide-index verification, not physical coefficient generation','check_sha256':hashlib.sha256((R/'check.rs').read_bytes()).hexdigest()}
with (R/'run.log').open('wb') as log:
 try:
  p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180)
  record['exit_code']=p.returncode
 except subprocess.TimeoutExpired: record['exit_code']=124
record['wall_seconds']=time.monotonic()-t
if (R/'result.json').exists(): record['result_sha256']=hashlib.sha256((R/'result.json').read_bytes()).hexdigest()
(R/'run.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
raise SystemExit(record['exit_code'])
