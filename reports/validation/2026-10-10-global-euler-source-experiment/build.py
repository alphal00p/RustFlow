import json,hashlib,subprocess,time,os
from pathlib import Path
ROOT=Path.cwd(); BASE=ROOT/'reports/validation/2026-10-10-global-euler-source-experiment'
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=BASE/'probe.rs';out=Path('/tmp/rustflow-global-euler-source-probe')
deps={'symbolica':'02ed301b2d4bd7f2','serde_json':'84255895f7dd0f90','symbolica_amflow':'bc386786e4706ace','rustred':'392075d0c6995fd1','bincode':'30fae4fd0bc9d1b5'}
libs={k:ROOT/f'target/release/deps/lib{k}-{v}.rlib' for k,v in deps.items()}
bindings={k:{'path':str(p),'sha256':h(p),'bytes':p.stat().st_size} for k,p in libs.items()}
cmd=['rustc','--edition=2024','--crate-name','global_euler_source_probe','-C','opt-level=0','-C','debuginfo=0']
for k,p in libs.items():cmd+=['--extern',f'{k}={p}']
cmd+=['-L',f'dependency={ROOT}/target/release/deps']
for p in sorted((ROOT/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L',f'native={p}']
cmd += [str(source),'-o',str(out)]
start=time.monotonic();r=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'CARGO_CRATE_NAME':'global_euler_source_probe'});(BASE/'build.log').write_text(r.stdout+r.stderr)
result={'scope':'Standalone diagnostic linked to frozen production RustFlow/native dependencies; no Cargo or production modification.','command':cmd,'dependencies':bindings,'source':{'path':str(source),'sha256':h(source)},'exit_code':r.returncode,'wall_seconds':time.monotonic()-start,'cargo_invoked':False,'compile_environment_override':{'CARGO_CRATE_NAME':'global_euler_source_probe'}}
assert all(h(libs[k])==v['sha256'] for k,v in bindings.items())
if not r.returncode:result['executable']={'path':str(out),'sha256':h(out)}
(BASE/'build-binding.json').write_text(json.dumps(result,indent=2)+'\n');print(r.stdout+r.stderr);print('exit',r.returncode);raise SystemExit(r.returncode)
