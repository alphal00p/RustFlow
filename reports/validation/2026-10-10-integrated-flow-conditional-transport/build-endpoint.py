import json,pathlib,hashlib,subprocess,os,time,shutil
R=pathlib.Path(__file__).resolve().parent; W=R.parents[2]; T=pathlib.Path('/tmp/rustflow-integrated-candidate-transport-20261010'); T.mkdir(exist_ok=True)
B=json.loads((R.parent/'2026-10-10-integrated-moment-expansion/gaussian-wide-build.json').read_text())
D={};args=['rustc','--edition=2024','--crate-name','integrated_candidate_transport','-C','opt-level=1','-C','debuginfo=0','-L',f'dependency={W}/target/release/deps']
for name in ['symbolica_amflow','serde_json']:
 x=B['dependencies'][name]; src=pathlib.Path(x['path']); dst=T/src.name
 assert hashlib.sha256(src.read_bytes()).hexdigest()==x['sha256']; shutil.copy2(src,dst)
 D[name]={'path':str(dst),'sha256':x['sha256']};args+=['--extern',f'{name}={dst}']
for native in [x/'out/lib' for x in (W/'target/release/build').glob('gmp-mpfr-sys-*')]+[pathlib.Path('/nix/store/rgnappqqc5vbq60gza5fflyk84sylwl6-python3-3.14.6/lib')]:
 if native.exists():args+=['-L',f'native={native}']
args+=['-C','link-arg=-Wl,-rpath,/nix/store/rgnappqqc5vbq60gza5fflyk84sylwl6-python3-3.14.6/lib',str(R/'endpoint.rs'),'-o',str(T/'endpoint')]
env=os.environ.copy();env['CARGO_CRATE_NAME']='integrated_candidate_transport'; start=time.monotonic()
with open(R/'endpoint-build.log','w') as f:r=subprocess.run(args,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
record={'command':args,'dependencies':D,'source_sha256':hashlib.sha256((R/'endpoint.rs').read_bytes()).hexdigest(),'exit_code':r.returncode,'seconds':time.monotonic()-start}
if not r.returncode: record['executable']={'path':str(T/'endpoint'),'sha256':hashlib.sha256((T/'endpoint').read_bytes()).hexdigest()}
(R/'endpoint-build.json').write_text(json.dumps(record,indent=2)+'\n');print(record);raise SystemExit(r.returncode)
