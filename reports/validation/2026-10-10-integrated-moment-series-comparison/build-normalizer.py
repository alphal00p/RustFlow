"""Standalone adapter compile against the producer's already captured library."""
from pathlib import Path
import hashlib,json,os,subprocess,time,sys
BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
BINDING=ROOT/'reports/validation/2026-10-10-integrated-moment-validation/producer-freeze/generator-build.json'
def bind(p):
 p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def main():
 attempt=sys.argv[1];out=BASE/attempt;out.mkdir()
 old=json.loads(BINDING.read_text());assert old['exit_code']==0
 for r in old['dependencies'].values():assert bind(r['path'])['sha256']==r['sha256']
 command=old['command'][:];command[command.index('--crate-name')+1]='integrated_moment_normalize'
 source=BASE/'source/normalize.rs';exe=Path('/tmp/rustflow-integrated-moment-normalize-20261010')
 command[-3:]=[str(source),'-o',str(exe)]
 env=dict(os.environ);env['CARGO_CRATE_NAME']='integrated_moment_normalize'
 record={'scope':'isolated report adapter; no Cargo/shared target/production edit','command':command,'source':bind(source),'producer_build_binding':bind(BINDING),'dependencies':old['dependencies'],'timeout_seconds':180}
 start=time.monotonic()
 with (out/'compile.log').open('wb') as log:
  try:r=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180);code=r.returncode
  except subprocess.TimeoutExpired:code=124
 record.update(exit_code=code,wall_seconds=time.monotonic()-start)
 if code==0:record['executable']=bind(exe)
 record['source_unchanged']=bind(source)==record['source']
 record['dependencies_unchanged']=all(bind(r['path'])['sha256']==r['sha256'] for r in old['dependencies'].values())
 (out/'binding.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps({k:record[k] for k in ['exit_code','wall_seconds','source_unchanged','dependencies_unchanged']}));raise SystemExit(code)
if __name__=='__main__':main()
