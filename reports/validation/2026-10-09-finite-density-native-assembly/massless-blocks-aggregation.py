from pathlib import Path
import sys,json,hashlib,gzip,shutil,datetime
from decimal import Decimal
ROOT=Path('/common/dev/rustflow_fermi');R=ROOT/'reports/validation/2026-10-09-finite-density-native-assembly'
sys.path.insert(0,str(ROOT/'tools/finite_density'))
from compare_massive_reference import complex_decimal,difference,ZERO
from compare_complete_massive_reference import component_keys,SECTORS,TARGETS

def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def rel(p):return str(p.relative_to(R))
newpath=R/'massless-blocks-massive-regression/prediction-28-80-12.json';oldpath=R/'full-sunset-frontier-sectors/prediction-28-80-12.json';refpath=R/'independent-reference/quadrature-64-60.json';refmeta=R/'independent-reference/independent-reference.json';new=read(newpath);old=read(oldpath);ref=read(refpath)['contributions'];meta=read(refmeta)
assert read(newpath.parent/'input.json')==meta['definition']==read(oldpath.parent/'input.json')
for obj in [new,old]:
 assert (obj['digits'],obj['series_order'],obj['occupied_start_scale'],obj['epsilon'])==(28,80,12,'4/5')
 assert obj['normalization']==meta['normalization']=='unscaled Euclidean amplitude'
 assert obj['guard_digits']==40 and obj['search_frontier_sectors'] is True and obj['independent_reference_comparisons']==0
 assert len(obj['contributions'])==4
comparison=[];regression=[];assembly=[]
for sector,nvals,ovals in zip(SECTORS,new['contributions']+[new['assembled']],old['contributions']+[old['assembled']]):
 for target,n,o in zip(TARGETS,nvals,ovals):
  actual=complex_decimal(n);previous=complex_decimal(o);keys=component_keys(target,sector);expected=(sum((Decimal(ref[k])for k in keys),ZERO),ZERO)
  comparison.append({'sector':sector,'target':target,'prediction':n,'reference_component_keys':keys,'reference_real':str(expected[0]),**difference(actual,expected),'imaginary_zero_passed':abs(actual[1])<=Decimal('1e-25')})
  regression.append({'sector':sector,'target':target,'identical_decimal_string':n==o,**difference(actual,previous)})
for i,target in enumerate(TARGETS):
 summed=tuple(sum((complex_decimal(c[i])[j] for c in new['contributions']),ZERO)for j in range(2))
 assembly.append({'target':target,**difference(summed,complex_decimal(new['assembled'][i]))})
assert all(c['passed']and c['imaginary_zero_passed']for c in comparison)and all(c['passed']for c in regression+assembly)
paths=[newpath,oldpath,refpath,refmeta,newpath.parent/'input.json',oldpath.parent/'input.json']
massive={'schema':1,'status':'passed','scope':'one current-source D12/5 profile, both original targets, vacuum/all cuts/full amplitude; independent reference check plus same-profile historical regression, not a new four-profile refinement study','configuration':[28,80,12],'guard_digits':40,'reference_comparison_count':10,'historical_same_profile_comparison_count':10,'independent_refinement_comparison_count':0,'supplied_oracle_numerical_records_compared':0,'reference_error_estimate':meta['error_estimate'],'criterion':'relative <=1e-12 for reference magnitude >=1e-20; absolute <=1e-25 otherwise; imaginary absolute <=1e-25','inputs_sha256':{rel(p):sha(p)for p in paths},'reference_comparisons':comparison,'historical_regressions':regression,'assembly_checks':assembly}
save(R/'massless-blocks-massive-regression-comparison.json',massive)
# Archive only completed roots owned by this task. No E7/source pilot files.
archivepath=R/'massless-blocks-numerical-archives.json'
assert not archivepath.exists(),'do not overwrite an existing archive manifest'
roots=['massless-blocks-fixed','massless-blocks-laurent','massless-blocks-d7','massless-blocks-massive-regression']
archive=[];selected=[]
for name in roots:
 root=R/name
 for p in root.rglob('*'):
  if p.is_file()and (p.suffix=='.bin' or (p.suffix=='.json' and any(part in ['native-closure','cut-0','cut-1','cut-01']for part in p.relative_to(root).parts[:-1]))):selected.append(p)
for name in roots+['massless-blocks-python']:
 p=R/(name+'-resources.log')
 if p.is_file():selected.append(p)
for p in sorted(set(selected)):
 rawsha=sha(p);rawsize=p.stat().st_size;dest=p.with_name(p.name+'.gz');assert not dest.exists()
 with p.open('rb')as a,gzip.open(dest,'wb',compresslevel=9)as b:shutil.copyfileobj(a,b)
 with gzip.open(dest,'rb')as b:assert hashlib.file_digest(b,'sha256').hexdigest()==rawsha
 archive.append({'original_path':rel(p),'original_sha256':rawsha,'original_bytes':rawsize,'archive_path':rel(dest),'archive_sha256':sha(dest),'archive_bytes':dest.stat().st_size,'decompression_verified':True});p.unlink()
save(archivepath,{'schema':1,'scope':'completed current-source massless and massive regression logs/native corpora only; numerical predictions retained as readable JSON','files':archive,'file_count':len(archive)})
# Current checkpoint provenance, preserving the original snapshot and explicit test delta.
snapshot=R/'massless-blocks-source-hashes.json';snap=read(snapshot);changes=[]
for name,oldhash in snap['files'].items():
 p=ROOT/name;newhash=sha(p)
 if newhash!=oldhash:changes.append({'path':name,'snapshot_sha256':oldhash,'current_sha256':newhash})
assert all(c['path']=='tests/finite_density_massless_sources.rs'for c in changes),changes
runs=[]
for name in roots+['massless-blocks-python']:
 resource=R/(name+'-resources.json');provenance=R/(name+'-resources.provenance.json')
 runs.append({'run':name,'resources':rel(resource),'resources_sha256':sha(resource),'result':read(resource),'executable_provenance':rel(provenance),'executable_provenance_sha256':sha(provenance),'captured_before_run':read(provenance)})
provenance={'schema':1,'scope':'current singleton-germ/rank-one-block/vacuum-zero implementation numerical checkpoint','source_snapshot':rel(snapshot),'source_snapshot_sha256':sha(snapshot),'snapshot_file_count':len(snap['files']),'base_commit':snap['base_commit'],'post_snapshot_changes':changes,'post_snapshot_change_scope':'test-only decoding of the nested stored measure identity before asserting equality; production sources unchanged','native_runs':runs,'archive_index':rel(archivepath),'archive_index_sha256':sha(archivepath),'supplied_oracle_numerical_records_compared':0,'four_loop_native_predictions_compared':0}
save(R/'massless-blocks-numerical-provenance.json',provenance)

def summary(name):
 path=R/name;v=read(path)
 def maxdiff(key,criterion):
  vals=[Decimal(c['relative_difference'] if criterion=='relative' else c['absolute_difference'])for c in v[key]if c['criterion']==criterion and (criterion!='relative' or c['relative_difference']is not None)]
  return str(max(vals))if vals else None
 return {'report':name,'sha256':sha(path),'status':v['status'],'reference_comparison_count':v['reference_comparison_count'],'independent_refinement_comparison_count':v['independent_refinement_comparison_count'],'largest_reference_relative_difference':maxdiff('comparisons','relative'),'largest_reference_small_absolute_difference':maxdiff('comparisons','absolute'),'largest_refinement_relative_difference':maxdiff('refinements','relative'),'largest_refinement_small_absolute_difference':maxdiff('refinements','absolute')}
fixed=summary('massless-blocks-fixed-reference-comparison.json');laurent=summary('massless-blocks-laurent-reference-comparison.json')
aggregate={'schema':1,'status':'two_loop_numerical_regressions_passed_four_loop_closure_unresolved','source_snapshot':rel(snapshot),'source_snapshot_sha256':sha(snapshot),'scope':'current sealed multi-loop singleton germ and rank-one-block origin evidence, native zero-vacuum certificates, with complete massless two-loop numerical regression and one complete massive regression profile','massless_fixed_dimension':{'dimension':'13/2',**fixed},'massless_laurent':{'orders':[-2,-1,0],**laurent},'massive_fixed_dimension_regression':{'dimension':'12/5','report':'massless-blocks-massive-regression-comparison.json','sha256':sha(R/'massless-blocks-massive-regression-comparison.json'),'reference_comparison_count':10,'historical_same_profile_comparison_count':10,'independent_refinement_comparison_count':0},'total_current_independent_reference_comparisons':80,'total_current_independent_refinement_comparisons':54,'same_profile_historical_regressions':10,'native_massless_occupied_basis_sizes':[2,2,3],'D7_specialized_sample':{'status':'unsupported_additional_indicial_resonance','stage':'double-cut large-eta Frobenius recurrence after both singleton occupied cuts succeeded and native zero-vacuum pruning skipped ordinary vacuum flow','accepted_predictions':0,'historical_difference':'previous checkpoint failed during ordinary recursive vacuum-boundary preparation; historical evidence remains unchanged'},'python':{'status':'passed','test_count':1,'resources':'massless-blocks-python-resources.json'},'four_loop_status':'E7 endpoint/origin class is admitted, but bounded weighted closure remains unresolved; no four-loop native numerical predictions or acceptance. E8 common contour and prism coupled two-virtual class remain unadmitted.','supplied_oracle_numerical_records_compared':0,'four_loop_native_predictions_compared':0,'provenance':'massless-blocks-numerical-provenance.json','archives':rel(archivepath),'reference_uncertainty':'Observed agreement and parameter refinements are empirical numerical evidence, not rigorous interval error bounds.'}
save(R/'massless-blocks-numerical-validation.json',aggregate)
print(json.dumps({'massless_sample':fixed,'massless_laurent':laurent,'massive_reference_comparisons':10,'same_profile_regressions':10,'archived':len(archive),'source_delta':changes},indent=2))
