from pathlib import Path
import json,hashlib,subprocess,time
R=Path(__file__).resolve().parent;ROOT=R.parents[2];build=json.loads((R/'census-build.json').read_text());exe=Path(build['executable']['path']);assert hashlib.sha256(exe.read_bytes()).hexdigest()==build['executable']['sha256'];out=R/'census.json';assert not out.exists()
cmd=[str(exe),str(ROOT/'examples/finite_density/massless_three_loop_chain.json'),str(out)];t=time.monotonic();record={'command':cmd,'timeout_seconds':180,'build_sha256':hashlib.sha256((R/'census-build.json').read_bytes()).hexdigest(),'references_loaded':False}
with (R/'census-run.log').open('wb') as log:
 try:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180);record['exit_code']=p.returncode
 except subprocess.TimeoutExpired:record['exit_code']=124
record['wall_seconds']=time.monotonic()-t
if out.exists():record['output']={'path':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'bytes':out.stat().st_size}
(R/'census-run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record));raise SystemExit(record['exit_code'])
