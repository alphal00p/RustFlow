"""Final integrity, deterministic archive and manifest for completed SZ2 gates."""
from pathlib import Path
import gzip,hashlib,json,sys
sys.dont_write_bytecode=True
from frozen_sources import BASE,ROOT,digest,verify_snapshot
PREFIX='singleton-physical-zero-v2'
read=lambda p:json.loads(p.read_text())
def main():
 assert not(BASE/'artifact-manifest.json').exists()
 sp=BASE/(PREFIX+'-source-hashes.json');snapshot=verify_snapshot(sp)
 bp=BASE/(PREFIX+'-build-provenance.json');build=read(bp)
 for name,row in build['files'].items():assert digest(ROOT/name)==row['sha256']
 gate=read(BASE/(PREFIX+'-root-gates.json'));assert gate['status']=='passed'and gate['build_provenance']['sha256']==digest(bp)
 comparisons={}
 for name in ['two-loop-fixed','two-loop-laurent','massive-regression']:
  stem=PREFIX+'-'+name;cp=BASE/(stem+'-comparison.json');cb=BASE/(stem+'-comparison-binding.json');run=BASE/(stem+'-binding.json');rr=BASE/(stem+'-resources.json')
  c,b,r,resource=map(read,[cp,cb,run,rr]);assert c['status']==b['status']=='passed'and b['comparison_sha256']==digest(cp)
  assert r['exit_code']==resource['exit_code']==0 and r['post_run_source_and_executable_unchanged']is True
  for file,h in b['inputs_sha256'].items():assert digest(ROOT/file)==h
  comparisons[name]={'comparison':cp.name,'comparison_sha256':digest(cp),'binding_sha256':digest(cb),'reference_checks':c['reference_comparison_count'],'refinement_checks':c['independent_refinement_comparison_count'],'historical_checks':c.get('historical_same_profile_comparison_count',0),'assembly_checks':len(c['assembly_checks']),'wall_seconds':resource['wall_seconds'],'peak_child_rss_kib':resource['peak_child_rss_kib']}
 # Keep every file required by the saved successful comparison raw. Archive
 # only obsolete proof/checkpoint payloads plus completed numerical stdout.
 protected=set()
 for name in comparisons:
  b=read(BASE/(PREFIX+'-'+name+'-comparison-binding.json'))
  for n in b['inputs_sha256']:
   p=ROOT/n
   if p.is_relative_to(BASE):protected.add(p)
 selected=[]
 for name in comparisons:
  directory=BASE/(PREFIX+'-'+name)
  selected.extend(p for p in (directory/'native-closure').rglob('*')if p.is_file()and p not in protected)
  selected.append(BASE/(PREFIX+'-'+name+'-resources.log'))
 rows=[]
 for p in sorted(set(selected)):
  q=p.with_name(p.name+'.gz');assert not q.exists();original=digest(p);size=p.stat().st_size
  with p.open('rb')as source,q.open('wb')as dest,gzip.GzipFile(filename='',mode='wb',fileobj=dest,mtime=0,compresslevel=6)as zipped:
   while block:=source.read(1024*1024):zipped.write(block)
  with gzip.open(q,'rb')as restored:assert hashlib.file_digest(restored,'sha256').hexdigest()==original
  rows.append({'logical_raw_path':str(p.relative_to(BASE)),'raw_sha256':original,'raw_bytes':size,'gzip_path':str(q.relative_to(BASE)),'gzip_sha256':digest(q),'gzip_bytes':q.stat().st_size,'restored_sha256':original});p.unlink()
 archive={'status':'all archived payloads restored byte-exactly','gzip_mtime':0,'files':rows,'raw_bytes':sum(r['raw_bytes']for r in rows),'gzip_bytes':sum(r['gzip_bytes']for r in rows),'preserved_comparison_inputs':{str(p.relative_to(BASE)):digest(p)for p in sorted(protected)}}
 (BASE/'archive-map.json').write_text(json.dumps(archive,indent=2)+'\n')
 # The numerical comparison bindings remain directly verifiable after storage
 # cleanup; no archived input is hidden behind a changed success assertion.
 for name in comparisons:
  for file,h in read(BASE/(PREFIX+'-'+name+'-comparison-binding.json'))['inputs_sha256'].items():assert digest(ROOT/file)==h
 verify_snapshot(sp)
 (BASE/'final-integrity.json').write_text(json.dumps({'status':'passed','source_assets_verified':len(snapshot['files']),'coherent_artifacts_verified':len(build['files']),'all_saved_comparison_inputs_unchanged':True,'source_snapshot_sha256':digest(sp),'build_provenance_sha256':digest(bp),'archive_files':len(rows),'source_change_after_last_run':False},indent=2)+'\n')
 summary={'status':'passed corrected coherent integration and all lower-loop comparisons','scope':'Distinct physical singleton endpoint-zero certificates with unchanged actual multicut AMF and all massive flows. No full three-loop/four-loop acceptance. Direct occupied-flow owner remains available; independently run HC cut-3 evidence is separate.',
 'root_gates':{'path':PREFIX+'-root-gates.json','sha256':digest(BASE/(PREFIX+'-root-gates.json')),'passed_tests':218,'native_passed':91,'root_passed':127,'ignored_native':1},'comparisons':comparisons,'source_snapshot_sha256':digest(sp),'build_provenance_sha256':digest(bp),'prior_failure':{'report':'../2026-10-10-singleton-physical-zero-integration','manifest_sha256':digest(BASE.parent/'2026-10-10-singleton-physical-zero-integration/artifact-manifest.json'),'reason':'one test expected explicit upper-support labels in original raised Cn/H0 targets; corrected test/comment only, algorithm unchanged'},'shared_host_timing_is_not_a_speed_benchmark':True,'archive':{'files':len(rows),'raw_bytes':archive['raw_bytes'],'gzip_bytes':archive['gzip_bytes']}}
 (BASE/'checkpoint-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 (BASE/'README.md').write_text('''# Singleton physical-zero integration: validated corrected attempt\n\n**Passed 218 coherent tests and all lower-loop numerical comparisons.** The gates include 91 native tests and 127 RustFlow/interface/source/capacity tests; one native test remains ignored. The build captures 13 artifacts against 2,485 source assets, including the direct massless-flow test target.\n\nThe massless two-loop fixed-D gate passed 40 independent value checks and 30 profile refinements. Its Laurent gate passed 30 independent coefficient checks and 24 refinements. The single massive profile passed 10 independent and 10 historical value checks plus two assembly checks; it is not a new precision-refinement gate. Exact source/build/executable identities and saved predictions precede every reference comparison.\n\nThe massless assembly retains all four physical cuts. Two singleton contributions use separately reported physical endpoint-zero certificates, with original coefficient conditions and finite-jet witnesses. They are not finite-eta source zeros, native closures, or transported AMF solutions. The remaining double cut retains actual native source closure, recursive boundary construction, transport and endpoint validation. Massive occupied cuts retain their actual flows. The direct PreparedOccupiedFlow path is unchanged.\n\nThe original test failure is preserved in the sibling initial-attempt report. The corrected fixture asserts the actual Cn/H0 representation and its n-1 mass jets, rather than expecting explicit upper-support indices in original input terms. The adjacent comment was corrected; the production algorithm was unchanged between attempts. Report-local metadata validators preserve historical reference definitions, numerical arithmetic and tolerances; old strict-AMF comparator files are unchanged.\n\n`checkpoint-summary.json` records counts/resources; `final-integrity.json` verifies the final source/build state and every successful comparison input. `archive-map.json` provides deterministic gzip and restored hashes for obsolete completed proof payloads and numerical logs. Current closed/provisional proofs and all directly bound comparison inputs remain raw. Timings are shared-host measurements, not speed comparisons.\n\nNo full three-loop or four-loop numerical acceptance is claimed. The independent HC selected cut-3 AMF control and the new standalone joint three-loop attempt are separate reports.\n''')
 files={str(p.relative_to(BASE)):{'sha256':digest(p),'bytes':p.stat().st_size}for p in BASE.rglob('*')if p.is_file()and '__pycache__'not in p.parts and p.name!='artifact-manifest.json'}
 (BASE/'artifact-manifest.json').write_text(json.dumps({'status':'frozen validated integration checkpoint','files':files,'logical_archived_files':rows,'production_modified_by_validation':False},indent=2)+'\n')
 print(json.dumps({'files':len(files),'manifest_sha256':digest(BASE/'artifact-manifest.json'),'archived':len(rows)}))
if __name__=='__main__':main()
