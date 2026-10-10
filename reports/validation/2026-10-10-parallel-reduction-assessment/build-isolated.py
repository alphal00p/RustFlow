#!/usr/bin/env python3
"""Compile the actual native crate with the isolated batch changes; no Cargo writes."""
import difflib, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path
ROOT=Path.cwd()
REPORT=ROOT/'reports/validation/2026-10-10-parallel-reduction-assessment'
NATIVE=ROOT/'vendor/rustred/crates/rustred-core'
DRAFT=REPORT/'implementation-draft/source/vendor/rustred/crates/rustred-core'
ISOLATED=Path('/tmp/rustred-ordered-unit-batch-20261010')
DEPS=Path('/tmp/rustred-ordered-unit-batch-deps-20261010')
kind, attempt=sys.argv[1:]
assert kind in ('tests','library')
prefix=REPORT/f'{kind}-build-{attempt}'
assert not prefix.with_suffix('.json').exists()
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def files(root):return {str(p.relative_to(root)):h(p) for p in sorted(root.rglob('*')) if p.is_file()}
old=json.loads((ROOT/'reports/validation/2026-10-10-occupation-order-experiment/build-binding.json').read_text())
base=files(NATIVE)
# Every attempt captures a fresh coherent copy; never writes production.
if ISOLATED.exists():shutil.rmtree(ISOLATED)
shutil.copytree(NATIVE,ISOLATED)
for p in DRAFT.rglob('*'):
    if p.is_file():
        q=ISOLATED/p.relative_to(DRAFT);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
proposal=files(ISOLATED)
DEPS.mkdir(exist_ok=True)
externs={}
for name,v in old['externs'].items():
    p=Path(v['path']);assert h(p)==v['sha256'],name
    q=DEPS/p.name
    if not q.exists():shutil.copy2(p,q)
    assert h(q)==v['sha256']
    externs[name]={'path':str(q),'original_path':str(p),'sha256':h(q)}
arities=Path(old['generated_arities']['path']);assert h(arities)==old['generated_arities']['sha256']
gen=DEPS/'generated';gen.mkdir(exist_ok=True);shutil.copy2(arities,gen/arities.name)
out=Path(f'/tmp/rustred-ordered-unit-batch-tests-20261010-{attempt}' if kind=='tests' else f'/tmp/librustred-ordered-unit-batch-20261010-{attempt}.rlib')
cmd=['timeout','600','env',f'OUT_DIR={gen}','CARGO_CRATE_NAME=rustred','rustc','--edition=2024','--crate-name','rustred']
cmd+=['--test'] if kind=='tests' else ['--crate-type','rlib']
cmd+=['--cfg','feature="native"','-C','opt-level=3','-C','debuginfo=0','-C','metadata=rustflow_ordered_unit_batch_20261010','-L',f'dependency={ROOT}/target/release/deps']
for name,v in externs.items():cmd+=['--extern',f'{name}={v["path"]}']
for p in sorted((ROOT/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L',f'native={p}']
cmd += [str(ISOLATED/'src/lib.rs'),'-o',str(out)]
record={'scope':'Isolated complete native crate and actual private memo/batch modules; no production or shared target writes. No old GuardedProgram alias.', 'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'base_sources':base,'isolated_sources':proposal,'draft_sources':files(DRAFT),'externs':externs,'generated_arities':{'path':str(gen/arities.name),'sha256':h(gen/arities.name)},'native_feature_fingerprint':json.loads((ROOT/'target/release/.fingerprint/rustred-392075d0c6995fd1/lib-rustred.json').read_text()),'command':cmd,'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
prefix.with_suffix('.pending.json').write_text(json.dumps(record,indent=2)+'\n')
start=time.monotonic()
with prefix.with_suffix('.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
record.update(exit_code=r.returncode,wall_seconds=time.monotonic()-start,production_sources_unchanged=files(NATIVE)==base,isolated_sources_unchanged=files(ISOLATED)==proposal)
if out.exists():record['output']={'path':str(out),'sha256':h(out),'bytes':out.stat().st_size}
prefix.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
prefix.with_suffix('.pending.json').unlink()
print(json.dumps({k:record[k] for k in ('exit_code','wall_seconds','production_sources_unchanged','isolated_sources_unchanged')}),flush=True)
raise SystemExit(r.returncode)
