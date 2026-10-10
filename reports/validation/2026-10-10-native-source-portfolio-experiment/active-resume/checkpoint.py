"""Pure metadata checks shared by the timeout-resume launcher and its controls."""
LIMITS={'max_requested':4096,'max_rounds':16,'max_frontier':1024,'max_rules':65536,
        'zero_attempts_per_round':8192,'depth':3,'ray_domains_per_requested_point':1,
        'point_domains_per_fallback':1}
SETTINGS={'PILOT_REQUESTS':'4096','PILOT_ROUNDS':'16','PILOT_FRONTIER':'1024',
          'PILOT_RULES':'65536','PILOT_ZERO_ATTEMPTS':'8192','PILOT_DEPTH':'3',
          'PILOT_RAY_DOMAINS':'1','PILOT_POINT_DOMAINS':'1','PILOT_REPLAY_OLD':'false'}

def validate(engine,config,binding,resources,progress,state,*,result_exists=False,closed_exists=False):
    assert engine in ['baseline','proposal'] and binding['engine']==engine
    assert binding['exit_code']==resources['exit_code']==124
    assert binding['post_run_inputs_unchanged'] is True
    assert not result_exists and not closed_exists
    assert binding['settings']==SETTINGS
    assert all(config[key]==value for key,value in LIMITS.items())
    assert config['old_proofs_reused'] is False and config['old_proof_optional_replay'] is None
    assert config['shifted_ordinary_slots']==[1,2,3,4]
    assert config['source_count']==62
    assert len(config['original_targets'])==2 and sum(map(len,config['original_targets']))==3
    assert progress and [row['round']for row in progress]==list(range(len(progress)))
    last=progress[-1]
    assert state['round']==last['round'] and state['round']+1<LIMITS['max_rounds']
    assert state['frontier_budget_exhausted'] is False
    sets={key:set(map(tuple,state[key]))for key in ['frontier','attempted_exact_points','next_needed','requested','pending_new_derivative_labels']}
    assert all(len(sets[key])==len(state[key])for key in sets)
    assert all(len(point)==12 for points in sets.values()for point in points)
    assert sets['requested']<=sets['attempted_exact_points']
    assert sets['frontier']<=sets['next_needed']
    assert sets['pending_new_derivative_labels']
    assert len(sets['frontier'])<=LIMITS['max_frontier']
    assert len(sets['attempted_exact_points']|sets['next_needed'])<=LIMITS['max_requested']
    assert state['rule_count']<=LIMITS['max_rules']
    assert last['frontier']==len(sets['frontier']) and last['rules']==state['rule_count']
    assert last['requested']==len(sets['requested'])
    return state['round']+1
