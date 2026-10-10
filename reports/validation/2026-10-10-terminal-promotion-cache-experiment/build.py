import json,hashlib,subprocess,time,os,sys
from pathlib import Path
R=Path.cwd();B=Path(__file__).resolve().parent;h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();mode=sys.argv[1];attempt=sys.argv[2];previous=json.loads((R/'reports/validation/2026-10-10-final-reduction-phase-diagnostic/build-binding.json').read_text());deps=previous['dependencies'];out=Path('/tmp/rustflow-terminal-promotion-'+mode+'-'+attempt)
cmd=['rustc','--edition=2024','--crate-name','terminal_promotion_experiment','-C','opt-level=0','-C','debuginfo=0','-A','unused_imports']
if mode=='tests':cmd+=['--test']
for name,item in deps.items():assert h(Path(item['path']))==item['sha256'];cmd+=['--extern',name+'='+item['path']]
cmd+=['-L','dependency='+str(R/'target/isolated-partial-boundary-validation-libs-20261010/deps')]
for p in sorted((R/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L','native='+str(p)]
source=B/'probe.rs';sha=h(source);cmd += [str(source),'-o',str(out)];start=time.monotonic();run=subprocess.run(cmd,capture_output=True,text=True,env={**os.environ,'CARGO_CRATE_NAME':'terminal_promotion_experiment'});(B/(mode+'-'+attempt+'-build.log')).write_text(run.stdout+run.stderr)
j={'command':cmd,'dependencies':deps,'source_sha256':sha,'mode':mode,'exit_code':run.returncode,'wall_seconds':time.monotonic()-start,'cargo_invoked':False};assert h(source)==sha
if not run.returncode:j['executable']={'path':str(out),'sha256':h(out)}
(B/(mode+'-'+attempt+'-build-binding.json')).write_text(json.dumps(j,indent=2)+'\n');print(run.stdout+run.stderr);raise SystemExit(run.returncode)
