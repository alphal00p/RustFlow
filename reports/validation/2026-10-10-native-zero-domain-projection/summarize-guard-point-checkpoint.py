#!/usr/bin/env python3
"""Summarize completed saved gates without generating any predictions or references."""
import gzip,json,re
from decimal import Decimal
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
P='guard-point'
snapshot_path=BASE/(P+'-source-hashes.json');snapshot=verify_snapshot(snapshot_path)
build_path=BASE/(P+'-build-provenance.json');build=json.loads(build_path.read_text())
assert build['source_snapshot_sha256']==digest(snapshot_path)

def binding(path):
 return {'path':str(path.relative_to(ROOT)),'sha256':digest(path)}
def load(name):return json.loads((BASE/name).read_text())
def logtext(name):
 path=BASE/name
 if path.exists():return path.read_text()
 return gzip.open(str(path)+'.gz','rt').read()
def check_comparison(name):
 d=load(name);assert d['status']=='passed'
 for path,sha in d.get('inputs_sha256',{}).items():assert digest(ROOT/path)==sha,path
 for key in ['comparisons','refinements','reference_comparisons','historical_regressions','assembly_checks']:
  for row in d.get(key,[]):
   if 'passed' in row:assert row['passed'],(name,row)
 return d
def maximum(rows,key,criterion=None):
 values=[Decimal(row[key]) for row in rows if row.get(key) is not None and (criterion is None or row.get('criterion')==criterion)]
 return str(max(values)) if values else None
def run_info(name):
 r=load(name+'-resources.json');b=load(name+'-binding.json')
 assert b['source_snapshot_sha256']==digest(snapshot_path)
 assert b['build_provenance_sha256']==digest(build_path)
 if name.endswith('massive-regression'):
  post=load(name+'-post-run-binding.json');assert post['unchanged_source_and_executable_verified'];assert post['resources_sha256']==digest(BASE/(name+'-resources.json'))
  exesha=post['executable_sha256']
 else:
  assert b['post_run_source_and_executable_unchanged'];assert b['exit_code']==r['exit_code'];exesha=b['executable_sha256']
 assert exesha==build['files']['target/release/deps/finite_density_runtime_flow-44939937ffe8534f']['sha256']
 matches=re.findall(r'test result: .*?finished in ([0-9.]+)s',logtext(name+'-resources.log'));assert len(matches)==1
 return {'resources':binding(BASE/(name+'-resources.json')),'launch':binding(BASE/(name+'-resources.provenance.json')),'binding':binding(BASE/(name+'-binding.json')),'exit_code':r['exit_code'],'harness_seconds':float(matches[0]),'wall_seconds':r['wall_seconds'],'peak_child_rss_kib':r['peak_child_rss_kib'],'executable_sha256':exesha}
fixed=check_comparison(P+'-two-loop-fixed-comparison.json');laurent=check_comparison(P+'-two-loop-laurent-comparison.json');massive=check_comparison(P+'-massive-regression-comparison.json');default=check_comparison(P+'-massive-default-path-comparison.json')
assert (fixed['reference_comparison_count'],fixed['independent_refinement_comparison_count'])==(40,30)
assert (laurent['reference_comparison_count'],laurent['independent_refinement_comparison_count'])==(30,24)
assert (massive['reference_comparison_count'],massive['historical_same_profile_comparison_count'],len(massive['assembly_checks']))==(10,10,2)
assert default['identical_decimal_strings']==default['numerical_comparisons']==10
massive_summary={'schema':1,'status':'passed','scope':'One D=12/5 profile (28 digits, order80, start12, guard40), both original targets, vacuum/all occupied cuts/full amplitude. Unchanged default closure configuration; no new precision-refinement claim.','normalization':'unscaled Euclidean amplitude','execution':run_info(P+'-massive-regression'),'comparison':binding(BASE/(P+'-massive-regression-comparison.json')),'immediate_predecessor_comparison':binding(BASE/(P+'-massive-default-path-comparison.json')),'independent_reference_values_compared':10,'historical_values_compared':10,'assembly_identities_checked':2,'immediate_predecessor_identical_decimal_strings':10,'largest_relative_reference_difference':maximum(massive['reference_comparisons'],'relative_difference','relative'),'reference_uncertainty':'Empirical independent quadrature/precision refinement, not a rigorous error interval.','independent_refinement_comparison_count':0,'native_occupied_basis_sizes':[6,6,8],'source_snapshot':binding(snapshot_path),'build_provenance':binding(build_path),'supplied_oracle_numerical_records_compared':0}
(BASE/(P+'-massive-regression-summary.json')).write_text(json.dumps(massive_summary,indent=2)+'\n')
name=P+'-three-loop-fixed';failure=load(name+'/failure.json');resources=run_info(name)
assert resources['exit_code']==101 and failure['numerical_acceptance'] is False
assert not list((BASE/name).glob('prediction*.json'))
text=logtext(name+'-resources.log');rows=re.findall(r'DifferentialClosure \{ round: (\d+), requested: (\d+), basis_size: (\d+), new_derivatives: (\d+) \}',text)
match=re.search(r'rounds=(\d+), provisional=(\d+), unresolved_terms=(\d+)',failure['error']);assert match
rounds,frontier,unresolved=map(int,match.groups());assert (rounds,frontier,unresolved)==(8,1025,9850)
trajectory=[int(row[2]) for row in rows]+[frontier];assert trajectory==[7,29,99,256,420,641,914,1025]
three={'schema':1,'status':'failed_bounded_native_closure','numerical_acceptance':False,'saved_prediction_count':0,'reference_comparison_count':0,'scope':'Full three-loop attempt stopped during preparation of first occupied singleton [0]; later cuts and all four configured numerical profiles were not evaluated.','cut_slots_at_failure':[0],'rounds':rounds,'frontier_limit':1024,'provisional_frontier':frontier,'frontier_trajectory':trajectory,'unresolved_terms':unresolved,'failure':binding(BASE/name/'failure.json'),'configuration':binding(BASE/name/'configuration.json'),'input':binding(BASE/name/'input.json'),'execution':resources,'source_snapshot':binding(snapshot_path),'build_provenance':binding(build_path),'compile_identities':build['compile_identities'],'improvement_claim':False}
(BASE/(P+'-three-loop-failure-summary.json')).write_text(json.dumps(three,indent=2)+'\n')
gates=load(P+'-root-gates.json');assert gates['status']=='passed' and gates['unique_test_count']==179
runs={}
for label,d in [('two-loop-fixed',fixed),('two-loop-laurent',laurent)]:
 runs[label]={'status':'passed','execution':run_info(P+'-'+label),'comparison':binding(BASE/(P+'-'+label+'-comparison.json')),'reference_comparison_count':d['reference_comparison_count'],'independent_refinement_comparison_count':d['independent_refinement_comparison_count'],'largest_relative_reference_difference':maximum(d['comparisons'],'relative_difference','relative'),'largest_relative_refinement_difference':maximum(d['refinements'],'relative_difference','relative'),'largest_absolute_difference_under_absolute_criterion':maximum(d['comparisons']+d['refinements'],'absolute_difference','absolute'),'relative_tolerance':d['relative_tolerance'],'small_magnitude_threshold':d['small_magnitude_threshold'],'small_absolute_tolerance':d['small_absolute_tolerance'],'imaginary_zero_absolute_tolerance':d['imaginary_zero_absolute_tolerance'],'profile_count':len(d['configurations'])}
status={'schema':2,'status':'two_loop_regression_scope_validated_three_loop_closure_incomplete','checkpoint_regression_scope_validated':True,'full_feature_validated':False,'source_snapshot':binding(snapshot_path),'source_delta':binding(BASE/(P+'-source-delta.json')),'build_provenance':binding(build_path),'regression_gates':binding(BASE/(P+'-root-gates.json')),'source_files':len(snapshot['files']),'build_artifacts':len(build['files']),'unique_regression_tests_passed':179,'ignored':1,'compile_identities':build['compile_identities'],'numerical_gates':runs,'massive_regression':binding(BASE/(P+'-massive-regression-summary.json')),'three_loop_failure':binding(BASE/(P+'-three-loop-failure-summary.json')),'runtime_configuration_distinction':'Massless fixed/Laurent and three-loop controls explicitly use guard_refinement.max_passes=1 with active closure/rays; prior failed projection controls used0. Massive regression retains its unchanged default closure, including max_passes=3 and inactive active-closure/ray options.','laurent_provenance_limitation':'Laurent predictions serialize total expansions, not per-sector construction metadata. The full-amplitude runtime implementation, exact source/build/executable binding, matching fixed-dimension cut coverage, and saved configuration bind that run; no additional per-sector Laurent evidence is invented.','historical_failures':{'unchanged_manifest':binding(BASE/'pre-refinement-artifacts.json'),'archive_resolution':binding(BASE/'released-native-archives.json'),'scope':'Prior valid conditional rule was unavailable on a coupled exceptional index face. Both failed massless controls and earlier three-loop failures remain preserved; the fresh point refinement supplies new exact native proofs instead of suppressing ConditionVanished.'},'three_loop_full_amplitude_validated':False,'native_four_loop_acceptance':False,'supplied_oracle_numerical_records_compared':0,'reference_uncertainty':'Empirical quadrature, precision, and native refinement evidence; no rigorous numerical interval bound.','timing_scope':'Prebuilt execution on shared host, compilation/Nix startup excluded; no isolated speedup claim.','summary_script':binding(Path(__file__).resolve())}
(BASE/(P+'-validation-status.json')).write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps({'status':status['status'],'regression_tests':179,'massless_fixed':[40,30],'massless_laurent':[30,24],'massive':[10,10,2],'three_loop':'failed: frontier1025,8rounds,9850unresolved; no predictions'}))
