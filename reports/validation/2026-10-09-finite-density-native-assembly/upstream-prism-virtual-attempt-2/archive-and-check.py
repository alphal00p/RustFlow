from pathlib import Path
import json,hashlib,shutil,os,sys,time,tempfile,subprocess,signal,datetime
source=Path('/tmp/rustflow-prism-reduced-wiswxb20')
out=Path('/common/dev/rustflow_fermi/reports/validation/2026-10-09-finite-density-native-assembly/upstream-prism-virtual-attempt-2');out.mkdir(exist_ok=True)
for name in ['oracle.wl','runner.py','wolfram_one_thread.py','wolfram_one_thread.wl','kernel.log','runner-metadata.json']:
 shutil.copyfile(source/name,out/name)
for name in ['result-18.wl','reduced-combinations.wl','exact-targets.wl','definition.json']:
 p=source/'work'/name
 if p.exists():shutil.copyfile(p,out/name)
shutil.copyfile('/tmp/check_prism_virtual_gamma.wl',out/'gamma-check.wl')
shutil.copyfile(__file__,out/'archive-and-check.py')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
status={'status':'partial_profile_completed_second_profile_timed_out','scope':'validation-only upstream virtual subgraph coefficients, not complete prism reference','completed_requested_digits':[18],'incomplete_requested_digits':[26],'complete_finite_density_references':0,'source_output':str(source),'definition_correction':'Both pair A [1,5] and pair B [1,7] jets plus pair C [5,7] raised neutral target are defined by the saved exact polynomial targets. The physical-input prose in the initial definition mentions pair A only; it is not the full five-target description.','artifacts_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
(out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
fd,license_path=tempfile.mkstemp(prefix='rustflow-prism-gamma-mathpass-');os.write(fd,b'!itplic.itp.unibe.ch\n');os.close(fd)
env=dict(os.environ);env.update(WOLFRAMINIT='-pwfile '+license_path,RUSTFLOW_WOLFRAM_KERNEL='/nix/store/m4xs334vv5qlr31q6ayqr5r2x8lx43bn-mathematica-15.0.1/bin/WolframKernel',RUSTFLOW_PRISM_REFERENCE_DIR=str(out),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PATH=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH',''))
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});(out/'wolfram_one_thread.py').chmod(0o755)
command=[str(out/'wolfram_one_thread.py'),'-noinit','-noprompt','-script',str(out/'gamma-check.wl')];start=time.monotonic()
try:
 with (out/'gamma-check.log').open('w') as log:
  proc=subprocess.Popen(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  try:code=proc.wait(timeout=60)
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGKILL);proc.wait();code=124
 metadata={'command':command,'exit_code':code,'wall_seconds':time.monotonic()-start,'timeout_seconds':60,'runtime_threads':1,'scope':'independent Gamma evaluation and comparison only','script_sha256':sha(out/'gamma-check.wl')}
 (out/'gamma-check-resources.json').write_text(json.dumps(metadata,indent=2)+'\n')
 print(json.dumps(metadata));print((out/'gamma-check.log').read_text()[-1000:])
finally:Path(license_path).unlink(missing_ok=True)
sys.exit(code)
