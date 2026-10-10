"""Check only static producer/comparator identities; never read prediction values."""
from pathlib import Path
import hashlib,json,tarfile,ast,datetime
ROOT=Path(__file__).resolve().parents[3];BASE=Path(__file__).resolve().parent
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
def main():
 cap=read(BASE/'capsule-binding.json');build=read(BASE/'build.json');harness=read(BASE/'harness-binding.json');pre=read(BASE/'numerical-prerequisites.json')
 for row in cap['static_files']:assert h(Path(row['path']))==row['sha256']
 for row in cap['evidence_copies']+pre['captured_files']:assert h(BASE/row['copy_path'])==row['sha256']
 assert h(BASE/'build-capsule.py')==cap['prepare_script_sha256']
 assert h(BASE/'build.json')==cap['standalone_build_sha256']
 assert h(BASE/'capture-prerequisites.py')==pre['capture_script_sha256']
 snap=read(BASE/'evidence/singleton-physical-zero-v2-source-hashes.json')
 archive=BASE/cap['source_archive']['path'];assert h(archive)==cap['source_archive']['sha256']
 with tarfile.open(archive,'r:gz')as tar:
  entries={m.name:hashlib.sha256(tar.extractfile(m).read()).hexdigest()for m in tar if m.isfile()};assert entries==snap['files']
  orig=tar.extractfile('tests/finite_density_runtime_flow.rs').read().decode()
 generated=(BASE/'joint-runtime.rs').read_text();assert generated.startswith(orig);joined=generated[len(orig):]
 for name,key in [('runtime_graph_full_amplitude','fixed_phase_body_sha256'),('runtime_graph_full_laurent','laurent_phase_body_sha256')]:
  body=orig.split('fn '+name+'() {',1)[1].split('\n#[test]',1)[0].rstrip()[:-1];prefix='\n    let run = Run::new();\n    let flow = run.full();';assert body.startswith(prefix);phase=body[len(prefix):]
  assert hashlib.sha256(phase.encode()).hexdigest()==harness[key];assert joined.count(phase)==1
 assert joined.count('let run = Run::new();')==joined.count('let flow = run.full();')==1
 assert h(BASE/'joint-runtime.rs')==harness['generated_sha256']==build['source_sha256']
 rootfp=read(BASE/'evidence/lib-symbolica_amflow.json')
 for name,row in build['direct_dependency_selection'].items():
  vals=[r[3]for r in rootfp['deps']if r[1]==name];assert vals==[row['root_dependency_fingerprint']]
  marker=BASE/'evidence'/Path(row['fingerprint_marker']).name
  assert h(marker)==row['fingerprint_sha256'];assert int.from_bytes(bytes.fromhex(marker.read_text().strip()),'little')==vals[0]
  assert h(marker.with_suffix('.json'))==row['metadata_sha256']
 plan=read(BASE/'comparison-plan.json')
 for key in ['frozen_comparator_and_reference_artifacts_sha256','report_local_metadata_validators_sha256']:
  for name,expected in plan[key].items():assert h(ROOT/name)==expected
 for path in BASE.glob('*.py'):ast.parse(path.read_text())
 old=ROOT/'tools/finite_density/compare_three_loop_reference.py'
 loops=lambda tree:[ast.dump(n,include_attributes=False)for n in ast.walk(tree)if isinstance(n,ast.For) and ast.unparse(n.target) in ['config','(i, expansion)','(j, power)','(sector, entries)','(i, target)','target','(key, actual)']]
 mains=[]
 for path in [old,BASE/'compare_three_loop_reference.py']:
  mains.append(next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef)and n.name=='main'))
 assert len(loops(mains[0]))==7 and loops(mains[0])==loops(mains[1])
 dig=read(BASE/'comparison-digest-binding.json');assert h(ROOT/dig['source_record'])==dig['source_record_sha256'];assert h(Path(dig['executable']['path']))==dig['executable']['sha256']
 files=['generate-harness.py','build-capsule.py','capture-prerequisites.py','run.py','capsule-binding.json','build.json','harness-binding.json','numerical-prerequisites.json','comparison-plan.json','comparison-readiness.json','joint_provenance.py','compare-joint.py','compare_three_loop_reference.py','certificate_validation.py','history_coverage.py','comparison-digest-binding.json','comparison-evidence/native-digest-build.json','check-static-comparison.py']
 report={'status':'static producer/capsule/comparison audit passed; complete numerical output still pending','reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files_sha256':{str((BASE/n).relative_to(ROOT)):h(BASE/n)for n in files},'restored_source_assets':len(entries),'exact_joined_phase_body_checks':2,'preparation_count':1,'cargo_dependency_fingerprints_exact':True,'development_library_symlinks_scope':'Used only during the frozen standalone link; static producer bytes and archived source/evidence are authenticated independently of later shared libraries.','packaging_failures_retained':['build-attempt00','build-attempt01'],'numeric_arithmetic_refinement_assembly_subloops_ast_equal':True,'numeric_subloops_checked':7,'tolerances_or_references_changed':False,'digest_executable_sha256':dig['executable']['sha256'],'comparison_scope':'Both complete fixed and Laurent phase markers, successful static post-run checks, exactly four fixed and five Laurent profiles, original input/normalization, explicit physical singleton zero certificates and actual multicut AMF proof/transactions, before reference comparison. No live production-source substitute.','remaining_limitations':['No result from the active joint computation is accepted by this static review.','The metadata validator binds the typed singleton theorem report; it is not an independent second theorem prover.','Individual Laurent cut coefficients are not inferred from the joined full-amplitude fit.'],'production_modified':False}
 p=BASE/'final-helper-comparison-review.json';assert not p.exists();p.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'sha256':h(p),'assets':len(entries)}))
if __name__=='__main__':main()
