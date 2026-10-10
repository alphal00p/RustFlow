import datetime,hashlib,json,pathlib,subprocess,time,os
R=pathlib.Path('/common/dev/rustflow_fermi');P=R/'reports/validation/2026-10-10-integrated-angular-review';S=P/'angular-final.rs'
D=pathlib.Path('/tmp/rustflow-integrated-fitter-20261010/libsymbolica-02ed301b2d4bd7f2.rlib');O=pathlib.Path('/tmp/rustflow-angular-independent-review-tests')
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
cmd=['rustc','--edition=2024','--test',str(P/'review.rs'),'-C','opt-level=2','-L','dependency='+str(R/'target/release/deps'),'--extern','symbolica='+str(D)]
for d in sorted((R/'target/release/build').glob('gmp-mpfr-sys-*/out/lib')):cmd+=['-L','native='+str(d)]
cmd+=['-o',str(O)];s=time.monotonic();p=subprocess.run(cmd,capture_output=True,text=True,timeout=180,env=dict(os.environ,CARGO_CRATE_NAME="angular_independent_review"));(P/'build.log').write_text(p.stdout+p.stderr)
v={'environment':{'CARGO_CRATE_NAME':'angular_independent_review'},'command':cmd,'exit_code':p.returncode,'wall_seconds':time.monotonic()-s,'source':{'path':str(S),'sha256':sha(S)},'review':{'path':str(P/'review.rs'),'sha256':sha(P/'review.rs')},'dependency':{'path':str(D),'sha256':sha(D)},'shared_target_writes':False}
if p.returncode==0:v['executable']={'path':str(O),'sha256':sha(O)}
(P/'build.json').write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(v));raise SystemExit(p.returncode)
