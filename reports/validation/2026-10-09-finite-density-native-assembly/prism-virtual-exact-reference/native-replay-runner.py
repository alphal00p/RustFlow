from pathlib import Path
import os,sys,subprocess,time,resource,hashlib,json,signal,datetime
out=Path('/tmp/rustflow-prism-native-replay-2');out.mkdir(exist_ok=True)
exe=Path('/tmp/finite_density_prism_virtual_replay');src=exe.with_suffix('.rs');candidate=Path('/tmp/rustflow-prism-symbolic-f5yyyumi/work/candidate-relations.json')
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
meta={'scope':'validation-only native exact source proof of optional external candidates; no AMF or numerical reference values','sha256':{str(p):digest(p) for p in [exe,src,candidate]},'timeout_seconds':180,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
env=dict(os.environ);env.update(RUSTFLOW_PRISM_VIRTUAL_REPLAY_REPORT=str(out),RUSTFLOW_PRISM_VIRTUAL_CANDIDATE=str(candidate));cmd=[str(exe),'optional_prism_candidates_are_native_exact_source_consequences','--exact','--ignored','--nocapture'];meta['command']=cmd
(out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
start=time.monotonic()
with (out/'run.log').open('w') as log:
 p=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 try:code=p.wait(timeout=180);timed_out=False
 except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait();code=124;timed_out=True
r=resource.getrusage(resource.RUSAGE_CHILDREN);meta.update(exit_code=code,timed_out=timed_out,wall_seconds=time.monotonic()-start,peak_rss_kib=r.ru_maxrss,user_seconds=r.ru_utime,system_seconds=r.ru_stime,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(out/'resources.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta));print((out/'run.log').read_text()[-3000:]);sys.exit(code)
