from pathlib import Path
import hashlib,json
root=Path('/common/dev/rustflow_fermi');base=root/'reports/validation/2026-10-10-native-source-portfolio-experiment';data=base/'role-subset-controls'
rel=json.loads((base/'control-artifact-relocation.json').read_text())
def resolve(p):
 s=str(p)
 if s.startswith(rel['old_absolute_prefix']+'/'):return Path(rel['new_absolute_prefix']+s[len(rel['old_absolute_prefix']):])
 return Path(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name,m in rel['files'].items():assert sha(data/name)==m['sha256']
paths={mode:data/f'baseline-{mode}-19-points.json' for mode in ['original','strict','cut-nonincreasing','index-dependent']}
runs={m:json.loads(p.read_text()) for m,p in paths.items()}
fixture=json.loads((base/'points.json').read_text());points=[tuple(p) for p in fixture];assert len(set(points))==19
full=runs['original'];source={r['original_ordinal']:r for r in full['original_source_rows']};assert len(source)==52
reports={};proofs={}
for mode,run in runs.items():
 assert not run['historical_rules_imported'] and not run['experiment_enabled']
 assert run['max_depth']==3 and run['max_domains_per_point']==1
 assert [tuple(p['target']) for p in run['probes']]==points
 assert run['roles']==full['roles'] and run['index_variable_positions']==full['index_variable_positions']
 assert len(run['complete_original_ordinal_map'])==run['source_count']
 for row in run['original_source_rows']:assert row==source[row['original_ordinal']]
 for p in run['probes']:
  assert p['persistence_roundtrip_replay']
  path=resolve(p['proof_path']);assert path.is_file();proofs[str(path.relative_to(base))]=sha(path)
  assert all(s['original_ordinal'] in run['complete_original_ordinal_map'] for s in p['rule_sources'])
 rows=[]
 for a,b in zip(full['probes'],run['probes']):
  applied=b['status'].startswith('Applied')
  improve=applied and len(b['rhs'])<len(a['rhs'])
  rows.append({'target':a['target'],'baseline_status':a['status'],'trial_status':b['status'],'baseline_rhs_terms':len(a['rhs']),'trial_rhs_terms':len(b['rhs']),'trial_applied':applied,'strictly_shorter':improve,'rhs_identical':a['rhs']==b['rhs'],'conditions_identical':a['conditions']==b['conditions'],'baseline_conditions':a['conditions'],'trial_conditions':b['conditions'],'trial_gaps':b['discovery_gaps'],'baseline_sources':a['rule_sources'],'trial_sources':b['rule_sources'],'baseline_stats':a['search_stats'],'trial_stats':b['search_stats'],'selected_arm':'trial' if improve else 'baseline','selected_rhs_terms':len(b['rhs']) if improve else len(a['rhs']),'selected_conditions':b['conditions'] if improve else a['conditions']})
 resources=json.loads((data/f'baseline-{mode}-19-points-resources.json').read_text());assert resources['exit_code']==0
 reports[mode]={'source_count':run['source_count'],'applied':sum(r['trial_applied'] for r in rows),'unresolved':sum(not r['trial_applied'] for r in rows),'source_original_ordinals':run['complete_original_ordinal_map'],'returned_rhs_terms':sum(r['trial_rhs_terms'] for r in rows),'strictly_shorter_applied_points':sum(r['strictly_shorter'] for r in rows),'retained_baseline_plus_shorter_trial_terms':sum(r['selected_rhs_terms'] for r in rows),'selected_changed_condition_lists':sum(r['strictly_shorter'] and not r['conditions_identical'] for r in rows),'rows':rows,'resources':resources}
assert reports['strict']['source_count']==18 and reports['cut-nonincreasing']['source_count']==34 and reports['index-dependent']['source_original_ordinals']==list(range(20,52))
combined=[]
for i,point in enumerate(points):
 candidates=[('original',full['probes'][i])]+[(m,runs[m]['probes'][i]) for m in ['strict','cut-nonincreasing','index-dependent']]
 candidates=[(m,p) for m,p in candidates if p['status'].startswith('Applied')]
 mode,p=min(candidates,key=lambda x:len(x[1]['rhs']))
 combined.append({'target':list(point),'selected_arm':mode,'rhs_terms':len(p['rhs']),'conditions':p['conditions'],'baseline_rhs_terms':len(full['probes'][i]['rhs'])})
report={'schema':1,'status':'fresh_native_subset_controls_completed','scope':'Four independent source-only searches/replay on exact same19points, baseline native library. Subset contexts preserve a recorded original-ordinal bijection but are separately indexed public-API contexts. This is not a full-corpus native portfolio implementation. Post-hoc shortest selection is diagnostic; conditions remain binding and no uniform parameter/domain coverage is asserted.','native_core_changed':False,'integral_order_changed':False,'source_equations_added':False,'known_values_used':False,'closure_or_physical_acceptance':False,'original_points':19,'original_rhs_terms':192,'selectors':reports,'combined_shortest_diagnostic':{'tie_break_order':['original','strict','cut-nonincreasing','index-dependent'],'rows':combined,'rhs_terms':sum(r['rhs_terms'] for r in combined),'scope':'Unimplemented four-arm selection among native-applied point rules; no claim of improved full closure or preserved exceptional loci.'},'proof_sha256':proofs,'input_sha256':{str(p.relative_to(base)):sha(p) for p in [*paths.values(),base/'points.json',base/'control-artifact-relocation.json',data/'baseline-probe-build.json',data/'probe.rs']}}
out=base/'selector-controls-comparison.json';assert not out.exists();out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({m:{k:v for k,v in x.items() if k in ['source_count','applied','unresolved','returned_rhs_terms','strictly_shorter_applied_points','retained_baseline_plus_shorter_trial_terms','selected_changed_condition_lists']} for m,x in reports.items()},indent=2));print('combined',report['combined_shortest_diagnostic']['rhs_terms'])
