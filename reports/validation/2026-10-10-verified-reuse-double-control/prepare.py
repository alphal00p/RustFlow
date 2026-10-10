"""Capture one coherent prebuilt runtime and immutable source/evidence capsule."""
import datetime,gzip,hashlib,importlib.util,json,os,shutil,subprocess,sys,tarfile
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path.cwd();HERE=Path(__file__).resolve().parent
VI=ROOT/'reports/validation/2026-10-10-verified-reuse-integration'
CAP=Path('/tmp/rustflow-verified-reuse-double-20261010')
def h(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
assert not CAP.exists() and not (HERE/'capsule-binding.json').exists()
spec=importlib.util.spec_from_file_location('frozen_sources',VI/'frozen_sources.py');owner=importlib.util.module_from_spec(spec);spec.loader.exec_module(owner)
sp=VI/'verified-reuse-source-hashes.json';bp=VI/'verified-reuse-build-provenance.json';gp=VI/'verified-reuse-root-gates.json'
snapshot=owner.verify_snapshot(sp);build=json.loads(bp.read_text());gates=json.loads(gp.read_text())
assert build['source_snapshot_sha256']==h(sp)
assert build['checks']['combined_release_build']=='passed'
assert gates['status']=='passed' and gates['build_provenance']['sha256']==h(bp) and gates['source_snapshot']['sha256']==h(sp)
assert gates['unique_test_count']==208
assert all(x['status']=='passed' and x['failed']==0 for x in gates['gates'])
assert {x['gate']for x in gates['gates']}=={'native','lib','flow-boundary','massless-sources','source-fingerprint','python','capacity','polynomial-sources'}
artifacts=[sp,bp,gp,VI/'verified-reuse-actual62-source-equivalence.json',VI/'run-three-loop-gate.py',VI/'run-resource-command.py']
for gate in gates['gates']:
 for item in gate['artifacts'].values():
  p=ROOT/item['path'];assert h(p)==item['sha256'];artifacts.append(p)
exes=[p for p in build['files']if Path(p).name.startswith('finite_density_runtime_flow-')];assert len(exes)==1
exe=ROOT/exes[0];assert h(exe)==build['files'][exes[0]]['sha256']
input_path='examples/finite_density/massless_three_loop_chain.json';assert h(ROOT/input_path)==snapshot['files'][input_path]
CAP.mkdir();(HERE/'evidence').mkdir()
shutil.copy2(exe,CAP/'runtime');os.chmod(CAP/'runtime',0o555)
shutil.copy2(ROOT/input_path,CAP/'input.json');shutil.copy2(ROOT/input_path,HERE/'input.json')
shutil.copy2(VI/'run-resource-command.py',CAP/'run-resource-command.py')
origin=[]
for p in dict.fromkeys(artifacts):
 q=HERE/'evidence'/p.name;shutil.copy2(p,q);assert h(q)==h(p)
 origin.append({'source_path':str(p.relative_to(ROOT)),'copy_path':str(q.relative_to(HERE)),'sha256':h(q),'bytes':q.stat().st_size})
# Source archive excludes binaries and output reports except explicit compile-time fixtures.
archive=HERE/'frozen-source-assets.tar.gz'
with archive.open('wb')as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0)as gz,tarfile.open(fileobj=gz,mode='w|')as tar:
 for name,digest in sorted(snapshot['files'].items()):
  p=ROOT/name;assert h(p)==digest
  info=tar.gettarinfo(str(p),arcname=name);info.mtime=0;info.uid=info.gid=0;info.uname=info.gname='';info.mode=0o644
  with p.open('rb')as f:tar.addfile(info,f)
with tarfile.open(archive,'r:gz')as tar:
 restored={}
 for member in tar:
  assert member.isfile();f=tar.extractfile(member);restored[member.name]=hashlib.file_digest(f,'sha256').hexdigest()
 assert restored==snapshot['files']
owner.verify_snapshot(sp)
assert h(CAP/'runtime')==build['files'][exes[0]]['sha256']
ldd=subprocess.run(['ldd',str(CAP/'runtime')],text=True,capture_output=True)
(HERE/'runtime-dynamic-libraries.txt').write_text(ldd.stdout+ldd.stderr)
assert ldd.returncode==0
record={'scope':'Static prebuilt VI runtime capsule for a fresh all-uncut-deformed selected double-cut control; no imported native rules, resume, analytic shortcut or reference values. Later ROOT source/target changes do not change this run.',
 'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'capsule_path':str(CAP),'original_runtime':exes[0],
 'source_snapshot_sha256':h(sp),'build_provenance_sha256':h(bp),'root_gates_sha256':h(gp),'compile_identities':gates['compile_identities'],'features':build['features'],
 'source_asset_count':len(snapshot['files']),'source_archive':{'path':archive.name,'sha256':h(archive),'bytes':archive.stat().st_size,'all_restored_member_hashes_match_snapshot':True},
 'evidence_copies':origin,'static_files':[{'path':str(p),'sha256':h(p),'bytes':p.stat().st_size}for p in sorted(CAP.iterdir())],
 'root_sources_verified_only_at_capture':True,'post_run_policy':'Check only immutable captured files, source archive and evidence copies; no live ROOT/source/target reads.',
 'prepare_script_sha256':h(Path(__file__)),'large_executable_storage':'/tmp only; do not commit executable'}
save(HERE/'capsule-binding.json',record)
print(json.dumps({'capsule':str(CAP),'runtime_sha256':h(CAP/'runtime'),'sources':len(snapshot['files']),'source_archive_bytes':archive.stat().st_size,'binding_sha256':h(HERE/'capsule-binding.json')}),flush=True)
