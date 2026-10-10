#!/usr/bin/env python3
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
PROOF=BASE.parent/'2026-10-10-partial-placement-support-proof'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
mode=sys.argv[1];assert mode in ['single-slot0','double-slots12','double-slots24']
build_path=BASE/'build-binding.json';build=json.loads(build_path.read_text());assert build['exit_code']==0
exe=Path(build['executable']['path']);assert sha(exe)==build['executable']['sha256']
assert sha(BASE/'pilot.rs')==build['source']['sha256']
assert sha(BASE/'prepare.rs')==build['prepare_source']['sha256']
for dep in build['dependencies'].values():assert sha(Path(dep['path']))==dep['sha256']
proof={str(p.relative_to(ROOT)):sha(p) for p in [PROOF/'proof.md',PROOF/'check.py',PROOF/'support-checks.json']}
proof_sha=hashlib.sha256(json.dumps(proof,sort_keys=True,separators=(',',':')).encode()).hexdigest()
physical=ROOT/'examples/finite_density/massless_three_loop_chain.json'
out=BASE/mode;assert not out.exists()
binding=BASE/f'{mode}-binding.json';assert not binding.exists()
env={'PILOT_REQUESTS':'4096','PILOT_ROUNDS':'16','PILOT_FRONTIER':'1024','PILOT_RULES':'65536','PILOT_ZERO_ATTEMPTS':'8192','PILOT_RAY_DOMAINS':'1','PILOT_POINT_DOMAINS':'1','PILOT_DEPTH':'3'}
command=['timeout','300',str(exe),str(physical),str(out),mode,proof_sha]
paths=[physical,exe,BASE/'pilot.rs',BASE/'prepare.rs',PROOF/'proof.md',PROOF/'check.py',PROOF/'support-checks.json']
inputs={str(p):sha(p) for p in paths}
report={'scope':'Fresh isolated native derivative-closure pilot with reviewed partial-support origin/scaleless sources. No old rule import, no numerical flow or production admission. Verified immutable scheduling unions; independent final source replay and persisted decode audit.',
        'command':command,'environment':env,'inputs_sha256':inputs,'proof_inputs_sha256':proof,'proof_binding_sha256':proof_sha,
        'build_binding_sha256':sha(build_path),'launcher_sha256':sha(Path(__file__))}
binding.write_text(json.dumps(report,indent=2)+'\n')
wrapper=ROOT/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'
r=subprocess.run([sys.executable,str(wrapper),str(BASE/f'{mode}-resources.json'),*command],cwd=ROOT,env={**os.environ,**env})
report.update(exit_code=r.returncode,inputs_unchanged=all(sha(Path(p))==h for p,h in inputs.items()))
assert report['inputs_unchanged'];binding.write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(r.returncode)
