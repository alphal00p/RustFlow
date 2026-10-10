from pathlib import Path
import json,hashlib,subprocess,os,time
ROOT=Path('/common/dev/rustflow_fermi');BASE=Path('/tmp/rustflow-partial-endpoint-20261010');REPORT=ROOT/'reports/validation/2026-10-10-partial-endpoint-capability-prototype'
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def newest(name):return max((ROOT/'target/release/deps').glob('lib'+name+'-*.rlib'),key=lambda p:p.stat().st_mtime_ns)
libs={n:newest(n) for n in ['symbolica','serde_json','serde','blake3','ahash','feynkit_graph','feynkit_kinematics']}
libs.update({'symbolica_amflow':BASE/'deps/libsymbolica_amflow-bc386786e4706ace.rlib','rustred':BASE/'deps/librustred-392075d0c6995fd1.rlib'})
bind={k:{'path':str(p),'sha256':h(p),'bytes':p.stat().st_size} for k,p in libs.items()}
cmd=['/nix/store/paaihvr59q1101200mabmk7y4yhw42g7-rustc-wrapper-1.97.1/bin/rustc','-C','linker=/nix/store/w88q44gqd1qg5wmkk7v0h97rpiqvam0l-gcc-wrapper-15.3.0/bin/cc','--edition=2024','--crate-name','partial_endpoint_prototype','--test','-C','opt-level=0','-C','debuginfo=0']
for k,p in libs.items():cmd+=['--extern',f'{k}={p}']
cmd+=['-L',f'dependency={BASE}/deps','-L',f'dependency={ROOT}/target/release/deps']
for p in sorted((ROOT/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L',f'native={p}']
cmd +=['-L','native=/nix/store/0bq2kpdgy8kimy82pzpjjw787h1sx6a9-python3-3.14.6-env/lib','-C','link-arg=-Wl,-rpath,/nix/store/0bq2kpdgy8kimy82pzpjjw787h1sx6a9-python3-3.14.6-env/lib']
cmd +=[str(BASE/'probe.rs'),'-o',str(BASE/'prototype-tests')]
sources_before={str(p.relative_to(BASE)):h(p) for p in sorted(BASE.glob("**/*.rs"))}
start=time.monotonic();r=subprocess.run(cmd,capture_output=True,text=True,env={**os.environ,'CARGO_CRATE_NAME':'partial_endpoint_prototype'});report={'command':cmd,'exit_code':r.returncode,'wall_seconds':time.monotonic()-start,'libraries':bind,'source_files':sources_before,'cargo_invoked':False}
assert sources_before=={str(p.relative_to(BASE)):h(p) for p in sorted(BASE.glob('**/*.rs'))}
log=REPORT/f'build-{len(list(REPORT.glob("build-*.json"))):02}';log.with_suffix('.log').write_text(r.stdout+r.stderr)
if not r.returncode:report['executable']={'path':str(BASE/'prototype-tests'),'sha256':h(BASE/'prototype-tests')}
log.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(r.stdout+r.stderr);print('EXIT',r.returncode);raise SystemExit(r.returncode)
