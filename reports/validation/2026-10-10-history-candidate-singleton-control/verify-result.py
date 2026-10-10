"""Validate saved selected-sector output, then compare the independent reference.

No source import, numerical integration, native discovery or live-root binding.
The frozen native BLAKE3 helper is used only as a file hash utility.
"""
import hashlib,json,subprocess,sys,tarfile
from pathlib import Path
from fractions import Fraction
BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'tools/finite_density'))
from compare_three_loop_reference import reference_adapter,value
from compare_massive_reference import difference
sys.path.insert(0,str(ROOT/'reports/validation/2026-10-10-history-candidate-limit-integration'))
from history_coverage import verify_history_coverage
read=lambda p:json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 cap=read(BASE/'capsule-binding.json');binding=read(BASE/'run-binding.json');resource=read(BASE/'resources.json')
 assert resource['exit_code']==binding['exit_code']==0 and binding['post_run_static_capsule_unchanged']is True
 assert sha(BASE/'capsule-binding.json')==binding['capsule_binding_sha256']
 for item in cap['static_files']:assert sha(Path(item['path']))==item['sha256']
 for item in cap['evidence_copies']:assert sha(BASE/item['copy_path'])==item['sha256']
 assert sha(BASE/cap['source_archive']['path'])==cap['source_archive']['sha256']
 snap=read(BASE/'evidence/history-candidate-limit-source-hashes.json')
 with tarfile.open(BASE/cap['source_archive']['path'],'r:gz')as archive:
  members=[m for m in archive.getmembers()if m.isfile()]
  assert len(members)==len(snap['files'])
  for member in members:
   assert hashlib.sha256(archive.extractfile(member).read()).hexdigest()==snap['files'][member.name]
 directory=BASE/'single3';prediction_path=directory/'prediction-18-60-8.json';prediction=read(prediction_path)
 assert prediction['cut_slots']==[3]and prediction['full_amplitude']is False
 assert [prediction[k]for k in ['digits','series_order','occupied_start_scale','guard_digits']]==[18,60,8,40]
 assert prediction['epsilon']=='-5/4'and prediction['independent_reference_comparisons']==0
 assert prediction['normalization']=='unscaled Euclidean cut contribution'
 report=prediction['report'];assert report['construction']=='weighted_amf'and report['basis_size']==7
 assert report['source_options']=={'policy':'polynomial-closure-v1','positive_compact_energy_powers':False,'free_virtual_zero_sectors':True}
 assert report['shifted_slots']==[0,1,2,4]and report['massless_endpoint']is not None
 assert read(directory/'input.json')==read(BASE/'input.json')
 proofs=directory/'native-closure';closed_path=proofs/'round-018-closed.json';closed=read(closed_path)
 provisional_path=proofs/'round-018-provisional.json';provisional=read(provisional_path)
 assert closed['status']=='closed' and closed['round']==18
 coverage=verify_history_coverage(closed,provisional,0,read(proofs/'retired-unresolved.json'))
 digest_record=ROOT/'reports/validation/2026-10-10-history-candidate-limit-integration/native-digest-build.json'
 digest_data=read(digest_record);tool=Path(digest_data['executable']['path']);assert sha(tool)==digest_data['executable']['sha256']
 def native_hash(path):
  fields=subprocess.check_output([str(tool),str(path)],text=True).split();assert int(fields[1])==path.stat().st_size;return fields[0]
 proof_checks=[]
 for path in sorted(proofs.glob('round-*.json')):
  row=read(path);program=path.with_suffix('.bin')
  assert program.stat().st_size==row['program_bytes']and native_hash(program)==row['program_blake3']
  link=row['active_state']['requested_discovery'];transaction=path.parent/link['path'];assert native_hash(transaction)==link['blake3']
  proof_checks.append({'metadata':str(path.relative_to(BASE)),'metadata_sha256':sha(path),'program':str(program.relative_to(BASE)),'program_sha256':sha(program),'program_blake3':row['program_blake3'],'transaction':str(transaction.relative_to(BASE)),'transaction_sha256':sha(transaction),'transaction_blake3':link['blake3']})
 assert len(proof_checks)==20
 assert closed['measure_id'].split(':',3)[1:3]==[cap['compile_identities']['RUSTRED_SOURCE_DIGEST'],cap['compile_identities']['DEPENDENCY_SOURCE_DIGEST']]
 # Only now read independent reference data, after the complete selected
 # prediction and its source/proof bindings have passed.
 refpath=ROOT/'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-independent-reference/reference.json'
 ref,refinputs=reference_adapter(refpath)
 assert read(ROOT/ref['source']['input'])==read(BASE/'input.json')
 samples=[r for r in ref['samples']if Fraction(r['dimension'])==Fraction(13,2)and Fraction(r['chemical_potential'])==1]
 assert len(samples)==1
 expected=samples[0]['single_cut_3'];assert len(report['values'])==len(expected)==2
 comparisons=[{'target':name,'prediction':actual,'reference':wanted,**difference(value(actual),value(wanted))}for name,actual,wanted in zip(['scalar','raised_original_numerator'],report['values'],expected)]
 assert all(r['passed']for r in comparisons)
 output={'status':'passed selected-sector native pipeline and two independent reference checks','scope':'Original three-loop cut [3] only; actual native closure, recursive hard boundary, transport and physical endpoint. One profile; no full three-loop/four-loop acceptance, refinement or speed claim. Exact zero agreement alone does not validate nonzero hard boundary coefficients.',
 'full_amplitude':False,'cut_slots':[3],'dimension':'13/2','profile':[18,60,8],'source_snapshot_sha256':cap['source_snapshot_sha256'],'build_provenance_sha256':cap['build_provenance_sha256'],'capsule_binding_sha256':sha(BASE/'capsule-binding.json'),'runtime_sha256':next(x['sha256']for x in cap['static_files']if Path(x['path']).name=='runtime'),'source_archive_restored_hashes_verified':len(snap['files']),
 'native_rounds':19,'final_basis_size':7,'history_coverage':coverage,'native_proof_hash_checks':proof_checks,'native_hash_utility':{'record':str(digest_record.relative_to(ROOT)),'record_sha256':sha(digest_record),'executable_sha256':sha(tool)},
 'prediction':{'path':str(prediction_path.relative_to(BASE)),'sha256':sha(prediction_path)},'comparisons':comparisons,'independent_reference_comparison_count':2,'independent_refinement_comparison_count':0,'supplied_oracle_records_compared':0,'reference_inputs_sha256':{str(p.relative_to(ROOT)):sha(p)for p in refinputs},'comparison_arithmetic':'unchanged historical Decimal difference function; relative 1e-12, small threshold 1e-20, absolute 1e-25',
 'resources':resource,'run_binding_sha256':sha(BASE/'run-binding.json'),'resource_sha256':sha(BASE/'resources.json'),'verifier_sha256':sha(Path(__file__))}
 (BASE/'selected-reference-comparison.json').write_text(json.dumps(output,indent=2)+'\n')
 print(json.dumps({'status':output['status'],'rounds':19,'basis':7,'references':2,'wall_seconds':resource['wall_seconds']}))
if __name__=='__main__':main()
