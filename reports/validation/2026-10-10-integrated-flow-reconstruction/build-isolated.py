#!/usr/bin/env python3
"""Compile only this report's generic exact fitter; no Cargo/shared output writes."""
import datetime, hashlib, json, os, pathlib, resource, shutil, subprocess, time
ROOT=pathlib.Path('/common/dev/rustflow_fermi')
REPORT=ROOT/'reports/validation/2026-10-10-integrated-flow-reconstruction'
CACHE=pathlib.Path('/tmp/rustflow-integrated-fitter-20261010')
CACHE.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
subprocess.run(['rustfmt','--edition','2024',str(REPORT/'source/fit.rs')],check=True)
deps={}
for name,basename in {'symbolica':'libsymbolica-02ed301b2d4bd7f2.rlib','serde_json':'libserde_json-84255895f7dd0f90.rlib','blake3':'libblake3-d00e8c4219217a2f.rlib'}.items():
    src=ROOT/'target/release/deps'/basename;dst=CACHE/basename
    if not dst.exists():shutil.copy2(src,dst)
    assert sha(src)==sha(dst)
    deps[name]={'path':str(dst),'sha256':sha(dst)}
common=['rustc','--edition=2024',str(REPORT/'source/fit.rs'),'-C','opt-level=2','-L','dependency='+str(ROOT/'target/release/deps')]
for name,b in deps.items():common+=['--extern',name+'='+b['path']]
for path in sorted((ROOT/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):common+=['-L','native='+str(path)]
for mode in ('test','binary'):
    exe=CACHE/('fitter-tests' if mode=='test' else 'fitter')
    cmd=common+(['--test'] if mode=='test' else [])+['-o',str(exe)]
    start=time.monotonic();now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    proc=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=180)
    log=REPORT/f'build-{mode}.log';log.write_text(proc.stdout)
    out={'started_utc':now,'command':cmd,'exit_code':proc.returncode,'wall_seconds':time.monotonic()-start,'source':{'path':str(REPORT/'source/fit.rs'),'sha256':sha(REPORT/'source/fit.rs')},'direct_dependencies':deps,'shared_target_writes':False,'log':str(log),'compiler':subprocess.check_output(['rustc','--version'],text=True).strip()}
    if proc.returncode==0:out['executable']={'path':str(exe),'sha256':sha(exe)}
    (REPORT/f'build-{mode}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'mode':mode,'exit':proc.returncode,'wall':out['wall_seconds'],'log':str(log)}),flush=True)
    if proc.returncode:raise SystemExit(proc.returncode)
