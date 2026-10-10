import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
assert not (HERE/'artifact-manifest.json').exists()
binding=read(HERE/'run-binding.json');resources=read(HERE/'resources.json');data=read(HERE/'profile/result.json');selection=read(HERE/'points.json');meta=read(HERE/'input.json')
assert binding['exit_code']==resources['exit_code']==0 and binding['post_run_inputs_unchanged']
assert all(h(Path(p))==sha for p,sha in binding['inputs_sha256'].items())
assert selection['terminals']==meta['frontier'] and [p['point']for p in selection['points'][:26]]==meta['active_state']['next_needed']
assert len(selection['points'])==len(data['points'])==42
assert all(a==b['selection']for a,b in zip(selection['points'],data['points']))
phases={p['phase']:p for p in data['timings']}
rows=[]
for i,p in enumerate(data['points']):
 r=p['reduction'];rows.append({'selection':p['selection'],'native_seconds':phases[f'point-{i:02}-reduce-native']['wall_seconds'],'conversion_seconds':phases[f'point-{i:02}-reduce-expression-json']['wall_seconds'],'rule_applications':r['rule_applications'],'terms':len(r['terms']),'unresolved':len(r['unresolved']),'conditions':len(r['conditions'])})
summary={'scope':data['scope'],'status':'bounded_profile_completed','source_replay_seconds':phases['decode_generated_and_replay_all_rules']['wall_seconds'],
 'final_frontier_labels':len(selection['terminals']),'next_needed_labels':26,'sampled_historical_labels':16,'complete_history_size':selection['history_count'],
 'rules':data['rules'],'limits':data['limits'],'slow_required_component':rows[25],
 'other_required_native_seconds':{'minimum':min(r['native_seconds']for r in rows[:25]),'maximum':max(r['native_seconds']for r in rows[:25])},
 'historical_sample_native_seconds':{'minimum':min(r['native_seconds']for r in rows[26:]),'maximum':max(r['native_seconds']for r in rows[26:])},
 'rows':rows,'resources':resources,'compiler_wall_seconds':read(HERE/'build-binding.json')['wall_seconds'],
 'inference_limits':['The required raised-target component itself costs substantial time; optional historical-map collection is not the sole demonstrated issue.','No full-history extrapolation, native closure, final original weighted-sum audit, boundary or numerical acceptance follows from these individual applications.','Per-rule scan, substitution, rational accumulation and condition retention are combined inside native reduction; finer attribution remains pending.','Timing samples overlap other shared-host work and do not establish an end-to-end speed ratio.'],
 'source_and_input_checks':'passed'}
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(HERE/'README.md').write_text('A copied final provisional three-loop singleton program was cold-decoded with complete source replay, then bound to its exact seven candidate terminals. The profile covers all 26 next-needed labels and 16 deterministic index-stratified historical labels. It does not claim the original weighted target sums or final closure.\n\nFull proof replay took 111.028 seconds. Required raised-target component `[2,1,1,1,1,0,-1,-1,0,0,0,0,0,0,0,0]` took 313.801 seconds inside native reduction: 15,021 applications, seven returned terms, no unresolved remainder, and 3,499 conditions. Its expression/JSON conversion took 0.0264 seconds. The other required labels took at most 0.0572 seconds each; the 16 historical samples ranged from 0.0018 to 3.34 seconds. These samples do not estimate the full historical pass.\n\nThe process finished in 431.694 seconds with 304,152 KiB peak child RSS, excluding compilation/Nix startup and prior native prototype compilation. This was shared-host work, not a controlled speed comparison. The tested immutable native prototype was linked statically; full cold decode and unchanged native reduction semantics were retained. No discovery, source-policy or production modification occurs in this report.\n\nThe original one-hour baseline remains a timeout without a closed checkpoint. This separate profile demonstrates expensive required-label reduction and motivates finer attribution; it does not replace the missing closure or numerical gate.\n')
files={str(p.relative_to(HERE)):{'sha256':h(p),'bytes':p.stat().st_size}for p in sorted(HERE.rglob('*'))if p.is_file()and'__pycache__'not in p.parts}
(HERE/'artifact-manifest.json').write_text(json.dumps({'scope':data['scope'],'files':files,'file_count':len(files)},indent=2)+'\n')
print(json.dumps({'files':len(files),'manifest_sha256':h(HERE/'artifact-manifest.json')}))
