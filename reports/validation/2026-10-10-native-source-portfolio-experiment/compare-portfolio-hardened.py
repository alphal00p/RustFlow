#!/usr/bin/env python3
"""Compare saved exact native point rules, including conservative condition coverage."""
import argparse,hashlib,importlib.util,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('baseline',type=Path);p.add_argument('proposal',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
assert not a.output.exists(),'preserve previous comparisons'
base=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('condition_polynomials',base/'condition-polynomials-v2.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
b=json.loads(a.baseline.read_text());q=json.loads(a.proposal.read_text())
assert b['mode']==q['mode']=='original' and b['full_corpus_preserved'] and q['full_corpus_preserved']
assert not b['experiment_enabled'] and q['experiment_enabled']
assert not b['historical_rules_imported'] and not q['historical_rules_imported']
assert b['source_count']==q['source_count']==b['original_source_count']==q['original_source_count']
assert b['original_source_rows']==q['original_source_rows']
assert b['max_depth']==q['max_depth']==3 and b['max_domains_per_point']==q['max_domains_per_point']==1
fixture=Path(b['point_fixture']);expected=json.loads(fixture.read_text());assert expected==json.loads(Path(q['point_fixture']).read_text())
assert len(set(map(tuple,expected)))==len(expected)
assert [r['target'] for r in b['probes']]==[r['target'] for r in q['probes']]==expected
rows=[];proofs={};sources=set(r['original_ordinal'] for r in b['original_source_rows'])
for old,new in zip(b['probes'],q['probes']):
 for record in [old,new]:
  assert record['persistence_roundtrip_replay']
  path=Path(record['proof_path']);assert path.is_file();proofs[str(path)]={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
  assert all(s['original_ordinal'] in sources for s in record['rule_sources'])
 oa=old['status'].startswith('Applied');na=new['status'].startswith('Applied')
 covered,matches=helper.condition_gate(old['conditions'],new['conditions'])
 stats=[s['portfolio'] for s in new['search_stats'] if s['portfolio'] is not None]
 for stat in stats:
  assert stat['baseline_rhs_terms']==len(old['rhs'])
  assert stat['selected_rhs_terms']==len(new['rhs'])
  assert abs(stat['total_elapsed']-stat['baseline_elapsed']-stat['policy_elapsed'])<1e-8
  for t in stat['trials']:
   assert t['attempted_rows']<=2048 and t['accepted_rows']<=512
   assert t['accepted_rows']+t['guard_rejected_rows']==t['attempted_rows']
   if t['completion']=='selected-strictly-shorter':assert t['exact_trace_rows']<=64 and t['exact_trace_terms']<=16384
  if stat['selected_arm']==0:assert old['rhs']==new['rhs'] and old['conditions']==new['conditions']
  else:assert na and len(new['rhs'])<len(old['rhs']) and covered
 rows.append({'target':old['target'],'baseline_status':old['status'],'proposal_status':new['status'],'baseline_rhs_terms':len(old['rhs']),'proposal_rhs_terms':len(new['rhs']),'strict_improvement':oa and na and len(new['rhs'])<len(old['rhs']),'rhs_identical':old['rhs']==new['rhs'],'conditions_identical':old['conditions']==new['conditions'],'baseline_conditions':old['conditions'],'proposal_conditions':new['conditions'],'condition_gate_pass':covered,'condition_matches':matches,'baseline_gaps':old['discovery_gaps'],'proposal_gaps':new['discovery_gaps'],'portfolio_stats':stats,'baseline_sources':old['rule_sources'],'proposal_sources':new['rule_sources']})
report={'schema':1,'status':'saved_native_point_rules_compared','scope':'Bounded isolated point-only original-source portfolio, independent native replay/roundtrip and exact rational-associate condition check. No active-closure or physical amplitude claim.','points':len(rows),'source_count':b['source_count'],'baseline_applied':sum(r['baseline_status'].startswith('Applied') for r in rows),'proposal_applied':sum(r['proposal_status'].startswith('Applied') for r in rows),'baseline_rhs_terms':sum(r['baseline_rhs_terms'] for r in rows),'proposal_rhs_terms':sum(r['proposal_rhs_terms'] for r in rows),'strict_improvements':sum(r['strict_improvement'] for r in rows),'identical_rhs':sum(r['rhs_identical'] for r in rows),'condition_gate_pass_points':sum(r['condition_gate_pass'] for r in rows),'rows':rows,'proofs':proofs,'inputs_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [a.baseline,a.proposal,fixture,Path(__file__),base/'condition-polynomials-v2.py']},'production_modified':False,'closure_or_numerical_acceptance':False}
a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['points','source_count','baseline_applied','proposal_applied','baseline_rhs_terms','proposal_rhs_terms','strict_improvements','condition_gate_pass_points']}))
