#!/usr/bin/env python3
import difflib, hashlib, json, os, subprocess, sys, time
from pathlib import Path
ROOT=Path.cwd(); BASE=ROOT/'reports/validation/2026-10-10-verified-program-reuse'
NATIVE=ROOT/'vendor/rustred/crates/rustred-core'; ISOLATED=Path('/tmp/rustred-verified-program-reuse-20261010')
PRIOR=ROOT/'reports/validation/2026-10-10-occupation-order-experiment/build-binding.json'
old=json.loads(PRIOR.read_text()); kind=sys.argv[1]
assert kind in ('library','tests')
OUT=Path('/tmp/librustred-verified-program-reuse-20261010.rlib' if kind=='library' else '/tmp/rustred-verified-program-reuse-tests-20261010')
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert not (BASE/f'{kind}-build.json').exists()
base={str(p.relative_to(NATIVE)):h(p) for p in sorted(NATIVE.rglob('*')) if p.is_file()}
proposal={str(p.relative_to(ISOLATED)):h(p) for p in sorted(ISOLATED.rglob('*')) if p.is_file()}
patch=''
for name, digest in proposal.items():
    if base.get(name)==digest:continue
    p=NATIVE/name; q=ISOLATED/name
    patch+=''.join(difflib.unified_diff(p.read_text().splitlines(True) if p.exists() else [],q.read_text().splitlines(True),fromfile='a/'+name if p.exists() else '/dev/null',tofile='b/'+name))
patchfile=BASE/'isolated-native.patch'
if patchfile.exists():assert patchfile.read_text()==patch
else:patchfile.write_text(patch)
for v in old['externs'].values():assert h(Path(v['path']))==v['sha256']
assert h(Path(old['generated_arities']['path']))==old['generated_arities']['sha256']
cmd=[x.replace('/tmp/rustflow-occupation-order-core-20261010',str(ISOLATED)).replace('/tmp/librustred-occupation-degree-first-20261010.rlib',str(OUT)).replace('metadata=rustflow_occupation_degree_first_20261010','metadata=rustflow_verified_program_reuse_20261010') for x in old['command']]
cmd.insert(cmd.index('rustc'),'CARGO_CRATE_NAME=rustred')
if kind=='tests':
    i=cmd.index('--crate-type');cmd[i:i+2]=['--test']
    for p in sorted((ROOT/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L',f'native={p}']
record={'scope':'Isolated immutable verified-program composition; no production source/build modification.','base_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'base_sources':base,'isolated_sources':proposal,'externs':old['externs'],'generated_arities':old['generated_arities'],'patch_sha256':h(patchfile),'command':cmd}
start=time.monotonic()
with (BASE/f'{kind}-build.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
record.update(exit_code=r.returncode,wall_seconds=time.monotonic()-start)
record['production_sources_unchanged']=all(h(NATIVE/n)==v for n,v in base.items())
record['isolated_sources_unchanged']=all(h(ISOLATED/n)==v for n,v in proposal.items())
if OUT.exists():record['output']={'path':str(OUT),'sha256':h(OUT),'bytes':OUT.stat().st_size}
(BASE/f'{kind}-build.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('exit_code','wall_seconds','production_sources_unchanged','isolated_sources_unchanged')}),flush=True)
raise SystemExit(r.returncode)
