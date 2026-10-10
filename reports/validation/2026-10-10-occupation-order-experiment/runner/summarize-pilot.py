"""Summarize only completed checkpoints; never infer closure after a timeout."""
import argparse, collections, datetime, hashlib, json, pathlib
parser=argparse.ArgumentParser();parser.add_argument('--tags', nargs='+', required=True);parser.add_argument('--output', required=True);args=parser.parse_args()
p=pathlib.Path('reports/validation/2026-10-10-occupation-order-experiment')
def meta(f):
 f=pathlib.Path(f);return {'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size}
result={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Bounded standalone full-source native active closure; no old rule import, no oracle data, no handpicked masters, no numerical transport.','baseline_library':'native v2 with explicit zero-domain projection; not the earlier canonical v1 library.','settings':{'rounds':16,'frontier':1024,'cumulative_requested':4096,'rules':65536,'depth':3,'ray_domains_per_point':1,'point_fallback_domains':1,'direct_zero_attempts_per_round':8192,'timeout_seconds':300},'runtime_scope':'Prebuilt opt-level0 wrappers with release dependencies. Compilation and Nix startup excluded; concurrent processes/standalone compiler load can affect timings. Interpret proof/count outcomes first.','pilot_source':meta(p/'runner/active-pilot.rs'),'arms':{}}
for tag in args.tags:
 directory=p/(tag+'-pilot');resources=json.loads((p/(tag+'-pilot-resources.json')).read_text());final=directory/'result.json';status=json.loads(final.read_text())['status'] if final.exists() else ('timeout' if resources['exit_code']==124 else 'process-failure')
 rounds=[]
 for f in sorted(directory.glob('round-*.json')):
  x=json.loads(f.read_text());original={tuple(t['indices']) for row in x['original_target_rows']for t in row};j=collections.Counter((i[0],i[9])for i in original)
  rounds.append({'round':x['round'],'requested':len(x['requested']),'cumulative_requested':len(x['attempted_exact_points']),'frontier':len(x['frontier']),'rules':x['rule_count'],'zero_first_rules':x['zero_first_rule_count'],'processed_actual_derivatives':len(x['derivative_rows']),'original_output_row_terms':list(map(len,x['original_target_rows'])),'original_output_distinct_leaves':len(original),'joint_cut_upper_surface_distribution':[{'cut':c,'upper_surface':h,'labels':n}for(c,h),n in sorted(j.items())],'conditions':len(x['conditions_encountered']),'seconds':x['seconds'],'checkpoint':meta(f)})
 result['arms'][tag]={'status':status,'closure_proved':status=='closed','resources':resources,'completed_rounds':rounds,'final_result':meta(final) if final.exists()else None,'interrupted_round_incomplete':resources['exit_code']==124,'native_decode_and_replay_required':True,'builds':str(p/(tag+'-runner-builds.json'))}
pathlib.Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
