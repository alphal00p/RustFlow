"""Summarize a completed immutable-runtime control without evaluating references."""
import gzip,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
binding=json.loads((HERE/'run-binding.json').read_text())
assert 'exit_code' in binding and binding['post_run_static_capsule_unchanged'] is True
resources=json.loads((HERE/'resources.json').read_text())
rounds=[]
metadata=set((HERE/'double/native-closure').glob('round-*-provisional.json'))
metadata.update(p.with_suffix('')for p in (HERE/'double/native-closure').glob('round-*-provisional.json.gz'))
for path in sorted(metadata):
 raw=path.read_bytes()if path.exists()else gzip.decompress(path.with_name(path.name+'.gz').read_bytes())
 value=json.loads(raw)
 rounds.append({'path':str(path.relative_to(HERE)),'sha256':hashlib.sha256(raw).hexdigest(),'round':value['round'],'status':value['status'],'rules':value['native_rule_count'],'frontier_count':len(value['frontier']),'requested_count':len(value['requested']),'program_bytes':value['program_bytes'],'program_blake3':value['program_blake3']})
predictions=sorted((HERE/'double').glob('prediction*.json'))
closed=sorted((HERE/'double/native-closure').glob('round-*-closed.json'))
record={'scope':binding['scope'],'run_binding_sha256':h(HERE/'run-binding.json'),'capsule_binding_sha256':h(HERE/'capsule-binding.json'),'resources_sha256':h(HERE/'resources.json'),'exit_code':binding['exit_code'],'resources':resources,'published_provisional_rounds':rounds,'closed_checkpoint_paths':[str(p.relative_to(HERE))for p in closed],'saved_predictions':[{'path':str(p.relative_to(HERE)),'sha256':h(p)}for p in predictions],'numerical_comparison_performed':False,'full_three_loop_acceptance':False,'four_loop_acceptance':False,'timing_comparison_claimed':False,'source_and_static_executable_unchanged':True,'incomplete_work_scope':'Only fully published checkpoint files are counted. A timeout may leave unfinished later work; no continuation or completed round is inferred.','summarizer_sha256':h(Path(__file__))}
(HERE/'outcome.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'exit_code':binding['exit_code'],'published_rounds':len(rounds),'last_round':rounds[-1]if rounds else None,'predictions':len(predictions),'closed_checkpoints':len(closed)},indent=2))
