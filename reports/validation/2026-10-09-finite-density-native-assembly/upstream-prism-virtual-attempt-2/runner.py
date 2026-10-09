from pathlib import Path
import datetime, hashlib, json, os, shutil, signal, subprocess, sys, tempfile, time

UPSTREAM=Path('/common/dev/amflow/target/upstream-performance/reproducible-build/amflow-26005517a288086c4cb4d1b26d829691bc088485')
KERNEL=Path('/nix/store/m4xs334vv5qlr31q6ayqr5r2x8lx43bn-mathematica-15.0.1/bin/WolframKernel')
KIRA=Path('/nix/store/rlxy7f90ifbzr77a6hvnsc7a13vmzzay-kira-3.1-git-dad964c/bin/kira')
FERMAT=Path('/nix/store/0582y8gjxv1j01phy7i7l3z07ky72bry-fermat-7.9b/bin/fer64')
SCRIPTS=Path('/common/dev/amflow/scripts')
HASHES={'AMFlow.m':'76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184','diffeq_solver/DESolver.m':'0dfa7f91bd8b8d160412dbc3a8c076413e488745bcefe58bc434999d03f76a46','ibp_interface/Kira/interface.m':'1c4476b32d8757ee18d6689d5b406bc00a2cfb00ad178926233f1768def7faaa'}
def sha(p):
    with Path(p).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
for name,expected in HASHES.items(): assert sha(UPSTREAM/name)==expected
out=Path(tempfile.mkdtemp(prefix='rustflow-prism-reduced-'))
print(out,flush=True)
shutil.copytree(UPSTREAM,out/'upstream')
work=out/'work';work.mkdir()
for name in ['wolfram_one_thread.py','wolfram_one_thread.wl']:
    shutil.copyfile(SCRIPTS/name,out/name)
(out/'wolfram_one_thread.py').chmod(0o755)
shutil.copyfile('/tmp/finite_density_prism_virtual_reduced_reference.wl',out/'oracle.wl')
shutil.copyfile(__file__,out/'runner.py')
(out/'upstream/ibp_interface/Kira/install.m').write_text(f'$KiraExecutable = "{KIRA}";\n$FermatExecutable = "{FERMAT}";\n$RatracerExecutable = "";\n')
fd,license_path=tempfile.mkstemp(prefix='rustflow-prism-mathpass-')
os.write(fd,b'!itplic.itp.unibe.ch\n');os.close(fd)
env=dict(os.environ)
env.update(WOLFRAMINIT='-pwfile '+license_path,RUSTFLOW_WOLFRAM_KERNEL=str(KERNEL),AMFLOW_ROOT=str(out/'upstream'),AMFLOW_CONTROLLED_KERNEL=str(out/'wolfram_one_thread.py'),AMFLOW_ORACLE_WORKDIR=str(work),FERMATPATH=str(FERMAT),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PATH=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH',''))
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
command=[str(out/'wolfram_one_thread.py'),'-noinit','-noprompt','-script',str(out/'oracle.wl')]
metadata={'scope':'independent validation-only prism virtual coefficient candidate; no production changes or answer data','source_sha256':HASHES,'driver_sha256':sha(out/'oracle.wl'),'runner_sha256':sha(__file__),'launcher_sha256':sha(out/'wolfram_one_thread.py'),'loader_sha256':sha(out/'wolfram_one_thread.wl'),'kernel':str(KERNEL),'kira':str(KIRA),'fermat':str(FERMAT),'command':command,'timeout_seconds':600,'cpu_affinity':sorted(os.sched_getaffinity(0)),'runtime_threads':1,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
start=time.monotonic()
try:
    with (out/'kernel.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            code=proc.wait(timeout=600);metadata['timed_out']=False
        except subprocess.TimeoutExpired:
            metadata['timed_out']=True
            try: os.killpg(proc.pid,signal.SIGTERM)
            except ProcessLookupError: pass
            try: proc.wait(timeout=20)
            except subprocess.TimeoutExpired: pass
            try: os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError: pass
            proc.wait();code=124
    metadata.update(exit_code=code,wall_seconds=time.monotonic()-start,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    if code==0:
        comparison=json.loads((work/'comparison.json').read_text())
        assert comparison['status']=='passed'
        metadata['comparison_sha256']=sha(work/'comparison.json')
    (out/'runner-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps({'output':str(out),'exit_code':code,'wall_seconds':metadata['wall_seconds']}),flush=True)
finally:
    Path(license_path).unlink(missing_ok=True)
sys.exit(code)
