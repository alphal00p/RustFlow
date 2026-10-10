import json,hashlib,pathlib
r=pathlib.Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=[];resources=[]
for mode in ['single3','double']:
 vs={d:json.loads((r/mode/f'depth-{d}/result.json').read_text())for d in [3,4,5]}
 for i,p in enumerate(vs[3]['probes']):
  row={'sector':mode,'selection':p['selection'],'historical_status':vs[3]['selection_observations'][i]['status'],'by_depth':{}}
  for d,v in vs.items():
   q=v['probes'][i];a=q['application']
   assert q['roundtrip'] and v['all_rhs_labels_admitted']
   row['by_depth'][str(d)]={'ray_status':q['ray']['status'],'point_fallback':q['point_fallback'],'final_status':a['status'],'rhs_terms':len(a['rhs']),'exact_rhs_status_and_ordered_conditions_equal_depth3':a==p['application'],'gaps':len(q['gaps'])}
  rows.append(row)
 for d in [3,4,5]:
  p=r/mode/f'depth-{d}-resources.json';v=json.loads(p.read_text());assert v['exit_code']==0
  resources.append({'sector':mode,'depth':d,'path':str(p.relative_to(r)),'sha256':sha(p),'wall_seconds':v['wall_seconds'],'peak_child_rss_kib':v['peak_child_rss_kib']})
out={'scope':'Fresh bounded native source-only discovery. Historical programs are decoded only to verify selected statuses, then discarded. No historical rules are used in fresh discovery. No period/closure claim.','result':'All six actual NoApplicableRule leaves and one historical Applied scalar control are covered by depth3 exact-point fallback; depth4/5 change none of their exact RHS/status/ordered conditions.','selection':'Lexicographic representatives of source-role/cut-power/surface strata from completed round003 frontiers; no reference values. Seven labels total.','budgets':{'ray_domains_per_label':1,'fallback_point_domains_per_label':1,'depths':[3,4,5],'timeout_seconds_per_case':300,'cases':6},'points':rows,'resources':resources,'higher_depth_benefit_observed':False,'all_roundtrips_passed':True,'program_count':21,'performance_scope':'Cases shared host resources. Timing includes historical replay selection verification; no speed superiority inferred.'}
(r/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
(r/'README.md').write_text('All seven selected labels have identical native outputs at depths 3, 4 and 5. Six are actual NoApplicableRule leaves in completed larger-frontier checkpoints; the seventh is an already-covered scalar control. Every fresh anchored-ray search misses at its one-domain cap, and exact-point fallback discovers an applicable rule at depth 3. Raising depth changes neither the exact RHS nor the ordered nonzero conditions.\n\nSingle-cut RHS counts are 4, 6, 8 and 1; double-cut counts are 6, 19 and 1. All 21 encoded native programs decode and replay, and all RHS labels stay admitted. Unresolved search boxes remain recorded. These results support point coverage for this sample; they do not establish a finite closed basis or recommend a larger search depth.\n\nThe selection is deterministic by cut-power/occupation strata, without numerical references. The saved historical program is independently replayed only to verify selection status, then discarded before fresh discovery. Each case uses its complete unchanged source corpus (54 or 232 rows), ray cap 1, point cap 1, sample seed 0 and a 300-second timeout. Source, library, executable, copied input and per-run bindings are retained. No production code or Cargo build changed.\n\nSee summary.json for exact selections and resources, input-bindings.json for original checkpoint copies, build-binding.json for dependency hashes, and each sector/depth result for full coefficients, conditions, proof sources and unresolved domains. All proof binaries remain raw. Runtimes are shared-host observations, including historical replay, and are not a performance comparison.\n')
files=[p for p in sorted(r.rglob('*'))if p.is_file()and p.name!='manifest.json']
manifest={'schema':'sha256-file-manifest-v1','files':[{'path':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':sha(p)}for p in files]}
(r/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('files',len(files),'manifest',sha(r/'manifest.json'))
