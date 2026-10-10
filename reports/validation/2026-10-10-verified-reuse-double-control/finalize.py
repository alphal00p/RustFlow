"""Freeze a completed failed control, preserving byte-exact logical payloads."""
import gzip,hashlib,json,shutil,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def h(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
binding=json.loads((HERE/'run-binding.json').read_text());outcome=json.loads((HERE/'outcome.json').read_text())
assert binding['exit_code']==124 and binding['post_run_static_capsule_unchanged']
assert not outcome['closed_checkpoint_paths']and not outcome['saved_predictions']
proofroot=HERE/'double/native-closure'
last=max(proofroot.glob('round-*-provisional.json'));protected={last,last.with_suffix('.bin')}
buildpath=ROOT/'reports/validation/2026-10-10-verified-reuse-integration/native-digest-build.json'
build=json.loads(buildpath.read_text());exe=Path(build['executable']['path']);assert h(exe)==build['executable']['sha256']
shutil.copy2(buildpath,HERE/'evidence/native-digest-build.json')
verified=[]
for path in sorted(proofroot.glob('round-*-provisional.json')):
 record=json.loads(path.read_text());program=path.with_suffix('.bin');result=subprocess.run([str(exe),str(program)],capture_output=True,text=True,check=True);blake,size=result.stdout.split();assert blake==record['program_blake3']and int(size)==record['program_bytes']
 tx=record['active_state']['requested_discovery'];tp=Path(tx['path']);tp=tp if tp.is_absolute() else path.parent/tp
 result=subprocess.run([str(exe),str(tp)],capture_output=True,text=True,check=True);tb,ts=result.stdout.split();assert tb==tx['blake3']
 verified.append({'metadata':str(path.relative_to(HERE)),'metadata_sha256':h(path),'program':str(program.relative_to(HERE)),'program_sha256':h(program),'program_blake3':blake,'program_bytes':int(size),'requested_discovery':str(tp.relative_to(HERE)),'requested_discovery_sha256':h(tp),'requested_discovery_blake3':tb,'requested_discovery_bytes':int(ts)})
save(HERE/'checkpoint-integrity.json',{'scope':'Byte/hash and transaction-link integrity only; no new symbolic replay or closed-system claim.','native_digest_build_sha256':h(buildpath),'native_digest_executable_sha256':h(exe),'rounds':verified})
paths=sorted(p for p in proofroot.iterdir()if p.is_file()and p not in protected)+[HERE/'resources.log']
rows=[]
for p in paths:
 assert p.suffix!='.gz';q=p.with_name(p.name+'.gz');assert not q.exists();rawhash=h(p);size=p.stat().st_size
 with p.open('rb')as src,q.open('wb')as dst,gzip.GzipFile(filename='',fileobj=dst,mode='wb',mtime=0)as gz:shutil.copyfileobj(src,gz)
 with gzip.open(q,'rb')as f:restored=hashlib.file_digest(f,'sha256').hexdigest()
 assert restored==rawhash
 rows.append({'original_path':str(p.relative_to(HERE)),'original_sha256':rawhash,'original_bytes':size,'gzip_path':str(q.relative_to(HERE)),'gzip_sha256':h(q),'gzip_bytes':q.stat().st_size,'restored_sha256':restored});p.unlink()
save(HERE/'archive-map.json',{'scope':'Completed verbose native payloads and raw log only; deterministic gzip mtime=0; original paths remain logical historical references. Last published metadata/program remain raw. No active reader remained in this report.','files':rows,'protected_raw':[{'path':str(p.relative_to(HERE)),'sha256':h(p),'bytes':p.stat().st_size}for p in sorted(protected)],'raw_bytes':sum(x['original_bytes']for x in rows),'gzip_bytes':sum(x['gzip_bytes']for x in rows)})
print(json.dumps({'archived_files':len(rows),'raw_bytes':sum(x['original_bytes']for x in rows),'gzip_bytes':sum(x['gzip_bytes']for x in rows),'protected_raw':[str(p.relative_to(HERE))for p in sorted(protected)]}))
