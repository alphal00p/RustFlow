#!/usr/bin/env python3
"""Link the tiny control runner to the frozen isolated native rlib; no Cargo."""
import datetime,hashlib,json,os,subprocess,time
from pathlib import Path
base=Path(__file__).resolve().parent;root=base.parents[2]
record_path=base/'focused-runner-build.json';out=Path('/tmp/rustred-source-portfolio-controls-20261010')
assert not record_path.exists() and not out.exists()
lib=Path('/tmp/librustred-source-portfolio-20261010.rlib');source=base/'run-focused-controls.rs'
def meta(p):return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
binding_path=base/'native-build/build-binding.json';binding=json.loads(binding_path.read_text());assert binding['exit_code']==0 and binding['output_sha256']==meta(lib)['sha256']
deps={'rustred':lib,'serde_json':root/'target/release/deps/libserde_json-84255895f7dd0f90.rlib'}
command=['rustc','--edition=2024','--crate-name','native_source_portfolio_controls','-C','opt-level=0','-C','debuginfo=0']
for name,path in deps.items():command+=['--extern',name+'='+str(path)]
command+=['-L','dependency='+str(root/'target/release/deps')]
for name in ['gmp-mpfr-sys-2e6d668657be102d','gmp-mpfr-sys-c789de89ed67d3db']:command+=['-L','native='+str(root/'target/release/build'/name/'out/lib')]
command+=[str(source),'-o',str(out)]
record={'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':meta(source),'dependencies':{n:meta(p) for n,p in deps.items()},'native_build_binding_sha256':meta(binding_path)['sha256'],'command':command,'cargo_invoked':False}
record_path.write_text(json.dumps(record,indent=2)+'\n');start=time.monotonic()
with (base/'focused-runner-build.log').open('w') as log:r=subprocess.run(command,cwd=root,env={**os.environ,'CARGO_CRATE_NAME':'native_source_portfolio_controls'},stdout=log,stderr=subprocess.STDOUT)
record['exit_code']=r.returncode;record['wall_seconds']=time.monotonic()-start
assert meta(source)==record['source']
for n,p in deps.items():assert meta(p)==record['dependencies'][n]
if r.returncode==0:record['executable']=meta(out)
record_path.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'exit_code':r.returncode,'wall_seconds':record['wall_seconds']}));raise SystemExit(r.returncode)
