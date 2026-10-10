from pathlib import Path
import os,json,hashlib,subprocess,time,glob
R=Path(__file__).resolve().parent; ROOT=R.parents[2]
def binding(p):
 p=Path(p);return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
deps={k:ROOT/'target/release/deps'/f for k,f in {'symbolica':'libsymbolica-02ed301b2d4bd7f2.rlib','serde_json':'libserde_json-84255895f7dd0f90.rlib','symbolica_amflow':'libsymbolica_amflow-bc386786e4706ace.rlib','rustred':'librustred-392075d0c6995fd1.rlib'}.items()}
cmd=['/nix/store/paaihvr59q1101200mabmk7y4yhw42g7-rustc-wrapper-1.97.1/bin/rustc','-C','linker=/nix/store/w88q44gqd1qg5wmkk7v0h97rpiqvam0l-gcc-wrapper-15.3.0/bin/cc','--edition=2024','--crate-name','integrated_moment_census','-C','opt-level=1','-C','debuginfo=0']
for k,v in deps.items():cmd+=['--extern',f'{k}={v}']
cmd+=['-L',f'dependency={ROOT}/target/release/deps']
for path in sorted(glob.glob(str(ROOT/'target/release/build/gmp-mpfr-sys-*/out/lib'))):cmd+=['-L',f'native={path}']
exe='/tmp/rustflow-integrated-moment-census-20261010';cmd+=[str(R/'source/census.rs'),'-o',exe]
record={'scope':'Isolated rustc compile; no Cargo/shared target mutation','command':cmd,'dependencies':{k:binding(v) for k,v in deps.items()},'sources':{str(p.relative_to(R)):binding(p) for p in sorted((R/'source').glob('*.rs'))},'runtime_limit_seconds':180}
env=os.environ.copy();env['CARGO_CRATE_NAME']='integrated_moment_census';t=time.monotonic()
with (R/'census-build.log').open('wb') as log:
 try: p=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180);record['exit_code']=p.returncode
 except subprocess.TimeoutExpired:record['exit_code']=124
record['wall_seconds']=time.monotonic()-t
if record['exit_code']==0:record['executable']=binding(exe)
(R/'census-build.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k in ['exit_code','wall_seconds','executable']}));raise SystemExit(record['exit_code'])
