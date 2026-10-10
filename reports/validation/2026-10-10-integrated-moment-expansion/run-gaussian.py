from pathlib import Path
import json,hashlib,subprocess,time,sys
R=Path(__file__).resolve().parent;ROOT=R.parents[2];build=json.loads((R/'gaussian-build.json').read_text());exe=Path(build['executable']['path']);assert hashlib.sha256(exe.read_bytes()).hexdigest()==build['executable']['sha256'];order=int(sys.argv[1]);out=R/f'gaussian-prefix-{order}';assert not out.exists();out.mkdir()
cmd=[str(exe),str(ROOT/'examples/finite_density/massless_three_loop_chain.json'),str(out),str(order)];t=time.monotonic();record={'command':cmd,'timeout_seconds':180,'build_sha256':hashlib.sha256((R/'gaussian-build.json').read_bytes()).hexdigest(),'references_loaded':False}
with (out/'run.log').open('wb') as log:
 try:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180);record['exit_code']=p.returncode
 except subprocess.TimeoutExpired:record['exit_code']=124
record['wall_seconds']=time.monotonic()-t
p=out/'coefficients.json'
if p.exists():record['output']={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
(out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record));raise SystemExit(record['exit_code'])
