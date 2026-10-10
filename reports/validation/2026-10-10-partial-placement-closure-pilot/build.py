import json,hashlib,subprocess,time,os
from pathlib import Path
ROOT=Path.cwd(); BASE=ROOT/'reports/validation/2026-10-10-partial-placement-closure-pilot'
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=BASE/'pilot.rs';out=Path('/tmp/rustflow-partial-placement-closure-pilot')
deps={'symbolica':'02ed301b2d4bd7f2','serde_json':'84255895f7dd0f90','symbolica_amflow':'bc386786e4706ace','rustred':'392075d0c6995fd1','bincode':'30fae4fd0bc9d1b5'}
libs={k:ROOT/f'target/release/deps/lib{k}-{v}.rlib' for k,v in deps.items()}
libs['rustred']=Path('/tmp/librustred-verified-program-reuse-20261010.rlib')
verified=ROOT/'reports/validation/2026-10-10-verified-program-reuse/library-build.json'
v=json.loads(verified.read_text());assert v['exit_code']==0 and h(libs['rustred'])==v['output']['sha256']
helper=BASE/'prepare.rs';helper_sha=h(helper)
bindings={k:{'path':str(p),'sha256':h(p),'bytes':p.stat().st_size} for k,p in libs.items()}
cmd=['rustc','--edition=2024','--crate-name','partial_placement_closure_pilot','-C','opt-level=0','-C','debuginfo=0']
for k,p in libs.items():cmd+=['--extern',f'{k}={p}']
cmd+=['-L',f'dependency={ROOT}/target/release/deps']
for p in sorted((ROOT/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L',f'native={p}']
cmd += [str(source),'-o',str(out)]
start=time.monotonic();r=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'CARGO_CRATE_NAME':'partial_placement_closure_pilot'});(BASE/'build.log').write_text(r.stdout+r.stderr)
result={'scope':'Standalone diagnostic linked to frozen production RustFlow/native dependencies; no Cargo or production modification.','command':cmd,'dependencies':bindings,'source':{'path':str(source),'sha256':h(source)},'exit_code':r.returncode,'wall_seconds':time.monotonic()-start,'cargo_invoked':False,'prepare_source':{'path':str(helper),'sha256':helper_sha},'isolated_native_build_binding':{'path':str(verified),'sha256':h(verified)},'final_audit':'independent with_terminals_replayed plus encode/decode target and derivative replay','compile_environment_override':{'CARGO_CRATE_NAME':'partial_placement_closure_pilot'}}
assert h(helper)==helper_sha
assert all(h(libs[k])==v['sha256'] for k,v in bindings.items())
if not r.returncode:result['executable']={'path':str(out),'sha256':h(out)}
(BASE/'build-binding.json').write_text(json.dumps(result,indent=2)+'\n');print(r.stdout+r.stderr);print('exit',r.returncode);raise SystemExit(r.returncode)
