import json,hashlib,subprocess,time,os,sys
from pathlib import Path
B=Path(__file__).resolve().parent;R=B.parents[2]
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
lib=json.loads((B/'library-binding.json').read_text());attempt=sys.argv[1];out=Path('/tmp/rustflow-partial-boundary-certificate-validation-'+attempt)
cmd=['rustc','--edition=2024','--test','--crate-name','partial_boundary_certificate_validation','-C','opt-level=0','-C','debuginfo=0','-A','dead_code','-A','unused_imports']
for name,binding in lib['dependencies'].items():
 p=Path(binding['preserved_path']);assert h(p)==binding['sha256'];cmd+=['--extern',f'{name}={p}']
cmd+=['-L','dependency='+lib['dependency_search']]
for p in sorted((R/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L','native='+str(p)]
cmd += [str(B/'harness/main.rs'),'-o',str(out)]
inputs={str(p.relative_to(B)):h(p) for p in (B/'harness').glob('*.rs')}
start=time.monotonic();res=subprocess.run(cmd,capture_output=True,text=True,env={**os.environ,'CARGO_CRATE_NAME':'partial_boundary_certificate_validation','CARGO_MANIFEST_DIR':str(R)})
(B/(attempt+'-build.log')).write_text(res.stdout+res.stderr)
report={'scope':'Isolated source-copy unit executable against pinned K libraries; no Cargo or shared target/release writes','command':cmd,'environment':{'CARGO_CRATE_NAME':'partial_boundary_certificate_validation','CARGO_MANIFEST_DIR':str(R)},'library_binding_sha256':h(B/'library-binding.json'),'inputs_sha256':inputs,'exit_code':res.returncode,'wall_seconds':time.monotonic()-start}
assert all(h(B/p)==v for p,v in inputs.items())
if not res.returncode:report['executable']={'path':str(out),'sha256':h(out)}
(B/(attempt+'-build-binding.json')).write_text(json.dumps(report,indent=2)+'\n');print(res.stdout+res.stderr);raise SystemExit(res.returncode)
