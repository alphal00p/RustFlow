#!/usr/bin/env python3
"""Check resume admission and rejection using completed parents; no native search."""
import copy
import hashlib
import json
from pathlib import Path
from checkpoint import validate

base=Path(__file__).resolve().parent
output=base/'checkpoint-controls.json'
assert not output.exists()
passed=[]
fixtures={}
for engine in ['baseline','proposal']:
    parent=base.parent/'active-pilot'/('polynomial62-'+engine)
    read=lambda name:json.loads((parent/name).read_text())
    progress=read('progress.json')
    fixture={'engine':engine,'config':read('configuration.json'),'binding':read('run-binding.json'),
             'resources':read('resources.json'),'progress':progress,
             'state':read(f"round-{progress[-1]['round']:03}.json")}
    assert validate(**fixture)==len(progress)
    passed.append(engine+'_completed_timeout_metadata_admitted')
    fixtures[engine]=fixture

def reject(name,edit):
    fixture=copy.deepcopy(fixtures['baseline'])
    edit(fixture)
    try:validate(**fixture)
    except (AssertionError,KeyError,TypeError):passed.append(name)
    else:raise AssertionError('invalid checkpoint admitted: '+name)

reject('unfinished_parent_rejected',lambda f:f['binding'].pop('exit_code'))
reject('completed_non_timeout_rejected',lambda f:f['resources'].update(exit_code=0))
reject('failed_non_timeout_rejected',lambda f:f['resources'].update(exit_code=101))
reject('wrong_engine_rejected',lambda f:f.update(engine='proposal'))
reject('unverified_parent_inputs_rejected',lambda f:f['binding'].update(post_run_inputs_unchanged=False))
reject('result_file_rejected',lambda f:f.update(result_exists=True))
reject('closed_program_rejected',lambda f:f.update(closed_exists=True))
reject('enlarged_round_budget_rejected',lambda f:f['config'].update(max_rounds=17))
reject('changed_request_cap_rejected',lambda f:f['config'].update(max_requested=4097))
reject('changed_search_depth_rejected',lambda f:f['config'].update(depth=4))
reject('changed_environment_rejected',lambda f:f['binding']['settings'].update(PILOT_ROUNDS='17'))
reject('historical_rules_import_rejected',lambda f:f['config'].update(old_proofs_reused=True))
reject('different_target_count_rejected',lambda f:f['config']['original_targets'].pop())
reject('noncontiguous_history_rejected',lambda f:f['progress'][0].update(round=1))
reject('last_state_mismatch_rejected',lambda f:f['state'].update(round=0))
reject('frontier_exhaustion_rejected',lambda f:f['state'].update(frontier_budget_exhausted=True))
reject('no_pending_obligation_rejected',lambda f:f['state'].update(pending_new_derivative_labels=[]))
reject('rule_exhaustion_rejected',lambda f:f['state'].update(rule_count=65537))
reject('duplicate_frontier_rejected',lambda f:f['state']['frontier'].append(f['state']['frontier'][0]))
reject('restored_frontier_not_needed_rejected',lambda f:f['state'].update(next_needed=[]))
reject('request_history_missing_rejected',lambda f:f['state'].update(attempted_exact_points=[]))
reject('request_exhaustion_rejected',lambda f:f['state'].update(next_needed=f['state']['next_needed']+[[i+10000]+[0]*11 for i in range(4097)]))
reject('global_round_exhaustion_rejected',lambda f:(f['progress'].extend(
    [{'round':i}for i in range(len(f['progress']),16)]),f['state'].update(round=15)))
report={'status':'passed','passed':len(passed),'controls':passed,
        'scope':'Metadata-only timeout-resume admission controls; no native replay or search executed. Native decode and exact saved-row replay occur in the compiled wrapper before continuing.',
        'source_sha256':{name:hashlib.sha256((base/name).read_bytes()).hexdigest()for name in ['checkpoint.py','test-checkpoint.py','pilot.rs']}}
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':len(passed),'status':'passed'}))
