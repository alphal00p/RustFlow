#!/usr/bin/env python3
"""Build one standalone pure-native point probe with an exact library binding; no Cargo."""
import argparse,datetime,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('tag',choices=['baseline','proposal']);a=p.parse_args()
base=Path(__file__).resolve().parent;root=base.parents[3]
out=Path('/tmp/rustred-role-subset-probe-'+a.tag)
binding_path=base/(a.tag+'-probe-build.json')
assert not out.exists() and not binding_path.exists()
lib=root/'target/release/deps/librustred-392075d0c6995fd1.rlib' if a.tag=='baseline' else Path('/tmp/librustred-direct-hit-lookahead-20261010.rlib')
def meta(path):return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
libmeta=meta(lib)
if a.tag=='baseline':assert libmeta['sha256']=='ac85472324303a84e909023a9aac6fe77f4993934985519fd88275cae0e1aad0'
else:
 d=json.loads((base/'native-build/build-binding.json').read_text());assert d['exit_code']==0 and d['output_sha256']==libmeta['sha256']
deps={'rustred':lib,'symbolica':root/'target/release/deps/libsymbolica-02ed301b2d4bd7f2.rlib','serde_json':root/'target/release/deps/libserde_json-84255895f7dd0f90.rlib','bincode':root/'target/release/deps/libbincode-30fae4fd0bc9d1b5.rlib'}
source=base/'probe.rs';command=['rustc','--edition=2024','--crate-name','direct_hit_lookahead_probe_'+a.tag,'-C','opt-level=0','-C','debuginfo=0']
if a.tag=='proposal':command+=['--cfg','direct_hit_lookahead']
for name,path in deps.items():command+=['--extern',name+'='+str(path)]
command+=['-L','dependency='+str(root/'target/release/deps')]
for directory in ['gmp-mpfr-sys-2e6d668657be102d','gmp-mpfr-sys-c789de89ed67d3db']:command+=['-L','native='+str(root/'target/release/build'/directory/'out/lib')]
command+=[str(source),'-o',str(out)]
record={'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'tag':a.tag,'source':meta(source),'dependencies':{name:meta(path) for name,path in deps.items()},'command':command,'cargo_invoked':False,'scope':'Standalone search/replay diagnostic only, using native arithmetic and elimination; no period or boundary evaluator.'}
binding_path.write_text(json.dumps(record,indent=2)+'\n');start=time.monotonic()
with (base/(a.tag+'-probe-build.log')).open('w') as log:result=subprocess.run(command,cwd=root,env={**os.environ,'CARGO_CRATE_NAME':'direct_hit_lookahead_probe_'+a.tag},stdout=log,stderr=subprocess.STDOUT)
record['exit_code']=result.returncode;record['wall_seconds']=time.monotonic()-start
assert meta(source)==record['source']
for name,path in deps.items():assert meta(path)==record['dependencies'][name]
if result.returncode==0:record['executable']=meta(out)
binding_path.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'tag':a.tag,'exit_code':result.returncode,'wall_seconds':record['wall_seconds']}));raise SystemExit(result.returncode)
