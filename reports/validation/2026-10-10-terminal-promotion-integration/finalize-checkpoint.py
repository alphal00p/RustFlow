"""Final integrity, deterministic archive and manifest for completed terminal-promotion gates."""
from pathlib import Path
import gzip,hashlib,json,sys
sys.dont_write_bytecode=True
from frozen_sources import BASE,ROOT,digest,verify_snapshot
PREFIX='terminal-promotion'
read=lambda p:json.loads(p.read_text())
def main():
 assert not(BASE/'artifact-manifest.json').exists()
 sp=BASE/(PREFIX+'-source-hashes.json');snapshot=verify_snapshot(sp)
 bp=BASE/(PREFIX+'-build-provenance.json');build=read(bp)
 for name,row in build['files'].items():assert digest(ROOT/name)==row['sha256']
 gate=read(BASE/(PREFIX+'-root-gates.json'));assert gate['status']=='passed'and gate['build_provenance']['sha256']==digest(bp)
 assert gate['unique_test_count']==231 and sum(r['passed']for r in gate['gates'])==231
 comparisons={}
 for name in ['two-loop-fixed','two-loop-laurent','massive-regression']:
  stem=PREFIX+'-'+name;cp=BASE/(stem+'-comparison.json');cb=BASE/(stem+'-comparison-binding.json');run=BASE/(stem+'-binding.json');rr=BASE/(stem+'-resources.json')
  c,b,r,resource=map(read,[cp,cb,run,rr]);assert c['status']==b['status']=='passed'and b['comparison_sha256']==digest(cp)
  assert r['exit_code']==resource['exit_code']==0 and r['post_run_source_and_executable_unchanged']is True
  for file,h in b['inputs_sha256'].items():assert digest(ROOT/file)==h
  comparisons[name]={'comparison':cp.name,'comparison_sha256':digest(cp),'binding_sha256':digest(cb),'reference_checks':c['reference_comparison_count'],'refinement_checks':c['independent_refinement_comparison_count'],'historical_checks':c.get('historical_same_profile_comparison_count',0),'assembly_checks':len(c['assembly_checks']),'wall_seconds':resource['wall_seconds'],'peak_child_rss_kib':resource['peak_child_rss_kib'],'memo':b['native_provenance']['unit_reduction_memo']}
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
 summary={'status':'passed coherent memo integration and all lower-loop comparisons','scope':'Optional bounded complete-unit reduction memo, default disabled; final complete source replay, original weighted-target audit and basis-derivative audit remain mandatory. No three-loop/four-loop numerical acceptance or speed claim.',
 'root_gates':{'path':PREFIX+'-root-gates.json','sha256':digest(BASE/(PREFIX+'-root-gates.json')),'passed_tests':231,'native_passed':101,'root_passed':130,'ignored_native':1},'comparisons':comparisons,'source_snapshot_sha256':digest(sp),'build_provenance_sha256':digest(bp),'predecessor':{'report':'../2026-10-10-singleton-physical-zero-integration-v2','manifest_sha256':digest(BASE.parent/'2026-10-10-singleton-physical-zero-integration-v2/artifact-manifest.json')},'reviewed_draft':{'path':'../2026-10-10-terminal-promotion-integration-draft/artifact-manifest.json','sha256':digest(BASE.parent/'2026-10-10-terminal-promotion-integration-draft/artifact-manifest.json')},'shared_host_timing_is_not_a_speed_benchmark':True,'archive':{'files':len(rows),'raw_bytes':archive['raw_bytes'],'gzip_bytes':archive['gzip_bytes']}}
 (BASE/'checkpoint-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 (BASE/'README.md').write_text("""# Optional terminal-promotion memo: validated integration

**Passed 231 coherent tests and all lower-loop numerical comparisons.** The gates include 101 native tests and 130 RustFlow/interface/source/capacity tests; one native test remains ignored. The coherent build binds 13 artifacts to 2,487 source assets. All ten new native memo tests and three root integration tests pass.

The memo-enabled massless two-loop fixed-D gate passed 40 independent reference checks and 30 profile refinements. Its Laurent gate passed 30 coefficient checks and 24 refinements. The massive one-profile gate kept the option disabled and passed 10 independent and 10 historical checks plus two assembly checks; no new massive precision-refinement claim is made.

The option reuses only complete unit calls in one immutable program, with a narrow replayed terminal-promotion handoff. It remains disabled by default. Encoded entries obey aggregate retained-data quotas; no descendant-tail substitution, rule discovery, source change, ordering change or epsilon specialization is performed. Full final source replay and all original weighted-target and basis-derivative audits remain. The exact root fixture compares enabled/disabled basis, matrix, target rows, conditions, program bytes and logical application counts. Native cancellation, guard failures, work limits and decoding errors remain explicit.

Saved structured counters distinguish logical, performed and memoized native reduction applications. They do not measure proof replay, terminal-admission checks, encoding/decoding or total CAS work. Peak owned-byte quotas exclude allocator overhead, the program and scratch; they are not process-memory limits. The fixed and Laurent controls exercise real memo hits. `checkpoint-summary.json` records exact counters and resources without a speed claim.

Physical singleton zero certificates remain separate from the actual nonzero multicut AMF. All original coefficient conditions, native multicut proofs, proof digests, transaction/history checks and independent reference definitions are preserved. The report-only digest helper initially lacked its copied Rust source; this precompile packaging failure and correction are recorded in `digest-preflight-correction.json`. The coherent production build and regression suites had no failure.

`final-integrity.json` verifies source/build identities and every successful comparison input. `archive-map.json` provides deterministic gzip and restored hashes for obsolete completed proof payloads and numerical logs. All directly comparison-bound proofs and metadata remain raw. Older frozen reports are unchanged. No new full-three-loop or four-loop numerical acceptance is claimed; the independent joint three-loop attempt belongs to a separate captured executable/report.
""")
 files={str(p.relative_to(BASE)):{'sha256':digest(p),'bytes':p.stat().st_size}for p in BASE.rglob('*')if p.is_file()and '__pycache__'not in p.parts and p.name!='artifact-manifest.json'}
 (BASE/'artifact-manifest.json').write_text(json.dumps({'status':'frozen validated integration checkpoint','files':files,'logical_archived_files':rows,'production_modified_by_validation':False},indent=2)+'\n')
 print(json.dumps({'files':len(files),'manifest_sha256':digest(BASE/'artifact-manifest.json'),'archived':len(rows)}))
if __name__=='__main__':main()
