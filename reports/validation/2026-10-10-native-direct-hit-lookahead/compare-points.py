#!/usr/bin/env python3
"""Compare saved native search presentations; never evaluate a physical period."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('baseline',type=Path);p.add_argument('proposal',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
if a.output.exists():raise ValueError('preserve prior comparison')
b=json.loads(a.baseline.read_text());q=json.loads(a.proposal.read_text())
assert b['mode']==q['mode'];assert b['source_count']==q['source_count']
assert b['original_source_rows']==q['original_source_rows'],'source corpus/order changed'
assert b['experiment_enabled'] is False and q['experiment_enabled'] is True
assert b['historical_rules_imported'] is False and q['historical_rules_imported'] is False
assert b['max_depth']==q['max_depth']==3 and b['max_domains_per_point']==q['max_domains_per_point']==1
bp={tuple(r['target']):r for r in b['probes']};qp={tuple(r['target']):r for r in q['probes']}
assert len(bp)==b['points']==len(b['probes']) and len(qp)==q['points']==len(q['probes'])
assert set(bp)==set(qp)
fixtures=[Path(b['point_fixture']),Path(q['point_fixture'])]
expected=[tuple(r) for r in json.loads(fixtures[0].read_text())]
assert len(set(expected))==len(expected) and set(bp)==set(expected)
assert json.loads(fixtures[0].read_text())==json.loads(fixtures[1].read_text())
rows=[];proofs={}
for point,old in bp.items():
 new=qp[point]
 assert old['persistence_roundtrip_replay'] and new['persistence_roundtrip_replay']
 for r in [old,new]:
  f=Path(r['proof_path']);assert f.is_file();proofs[str(f)]={'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size}
 for stat in new['search_stats']:
  extra=stat['lookahead']
  if extra is None:continue
  assert 0<=extra['extra_attempted_rows']<=512
  if extra['completion']=='strictly-smaller-indirect':
   assert extra['indirect_rhs_terms']<extra['first_direct_rhs_terms']
   assert extra['indirect_trace_rows']<=64 and extra['indirect_trace_terms']<=16384
   assert not stat['direct_hit']
  if 'fallback' in extra['completion']:assert stat['direct_hit']
 oa=old['status'].startswith('Applied');na=new['status'].startswith('Applied')
 rows.append({'target':list(point),'baseline_status':old['status'],'proposal_status':new['status'],'baseline_rhs_terms':len(old['rhs']),'proposal_rhs_terms':len(new['rhs']),'both_applied':oa and na,'strictly_shorter_applied_rhs':oa and na and len(new['rhs'])<len(old['rhs']),'rhs_identical':old['rhs']==new['rhs'],'conditions_identical':old['conditions']==new['conditions'],'baseline_conditions':old['conditions'],'proposal_conditions':new['conditions'],'baseline_gaps':old['discovery_gaps'],'proposal_gaps':new['discovery_gaps'],'baseline_sources':old['rule_sources'],'proposal_sources':new['rule_sources'],'baseline_search_stats':old['search_stats'],'proposal_search_stats':new['search_stats']})
report={'schema':1,'status':'saved_search_results_compared','scope':'Exact native point-rule presentation comparison, including exceptional conditions and certification/discovery gaps. Shorter RHS does not imply equal conditional coverage, differential closure, or a physical amplitude.','mode':b['mode'],'points':len(rows),'baseline_applied':sum(r['baseline_status'].startswith('Applied') for r in rows),'proposal_applied':sum(r['proposal_status'].startswith('Applied') for r in rows),'strictly_shorter_applied_rhs_points':sum(r['strictly_shorter_applied_rhs'] for r in rows),'identical_rhs_points':sum(r['rhs_identical'] for r in rows),'identical_condition_lists':sum(r['conditions_identical'] for r in rows),'baseline_total_rhs_terms':sum(r['baseline_rhs_terms'] for r in rows),'proposal_total_rhs_terms':sum(r['proposal_rhs_terms'] for r in rows),'rows':rows,'proof_files':proofs,'inputs_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [a.baseline,a.proposal,Path(__file__),*fixtures]},'production_recommendation':False,'physical_numerical_acceptance':False}
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','mode','points','baseline_applied','proposal_applied','strictly_shorter_applied_rhs_points','baseline_total_rhs_terms','proposal_total_rhs_terms']}))
