#!/usr/bin/env python3
"""Compile a standalone resume wrapper against a frozen native library; no Cargo."""
import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('engine',choices=['baseline','proposal'])
a=p.parse_args()
base=Path(__file__).resolve().parent
root=base.parents[3]
source=base/'pilot.rs'
out=base/(a.engine+'-build.json')
exe=Path('/tmp/rustred-source-portfolio-resume-'+a.engine+'-20261010')
assert not out.exists() and not exe.exists()
meta=lambda path:{'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
point_binding=json.loads((base.parent/(a.engine+'-generic-probe-build.json')).read_text())
assert point_binding['exit_code']==0
dependencies={name:Path(value['path']) for name,value in point_binding['dependencies'].items()}
assert all(meta(path)['sha256']==point_binding['dependencies'][name]['sha256'] for name,path in dependencies.items())
command=['rustc','--edition=2024','--crate-name','portfolio_resume_'+a.engine,'-C','opt-level=0','-C','debuginfo=0']
if a.engine=='proposal':command+=['--cfg','source_portfolio']
for name,path in dependencies.items():command+=['--extern',name+'='+str(path)]
command+=['-L','dependency='+str(root/'target/release/deps')]
for name in ['gmp-mpfr-sys-2e6d668657be102d','gmp-mpfr-sys-c789de89ed67d3db']:
    command+=['-L','native='+str(root/'target/release/build'/name/'out/lib')]
command+=[str(source),'-o',str(exe)]
record={'engine':a.engine,'source':meta(source),'dependencies':{k:meta(v)for k,v in dependencies.items()},
        'command':command,'cargo_invoked':False,'launcher':meta(Path(__file__)),
        'scope':'Standalone validation-only resume wrapper; unchanged native search and original global limits.'}
out.write_text(json.dumps(record,indent=2)+'\n')
start=time.monotonic()
with (base/(a.engine+'-build.log')).open('w')as log:
    result=subprocess.run(command,cwd=root,env={**os.environ,'CARGO_CRATE_NAME':'portfolio_resume_'+a.engine},stdout=log,stderr=subprocess.STDOUT)
record['exit_code']=result.returncode
record['wall_seconds']=time.monotonic()-start
assert meta(source)==record['source']
assert all(meta(path)==record['dependencies'][name]for name,path in dependencies.items())
if result.returncode==0:record['executable']=meta(exe)
out.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'engine':a.engine,'exit_code':result.returncode,'wall_seconds':record['wall_seconds']}))
raise SystemExit(result.returncode)
