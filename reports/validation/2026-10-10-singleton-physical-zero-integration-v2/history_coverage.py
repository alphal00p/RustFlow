"""Exact coverage checks; unassessed history is never promoted to a proof."""

def verify_history_coverage(closed, provisional, expected_limit, retired):
    assert expected_limit is None or (type(expected_limit) is int and expected_limit >= 0)
    assert closed['status']=='closed' and provisional['status']=='provisional'
    assert closed['active_target_closure'] is True and provisional['active_target_closure'] is True
    assert closed['max_history_candidate_maps']==provisional['max_history_candidate_maps']==expected_limit
    # Only stopping-set proof bytes, status and final coverage reporting differ.
    def shared(row):
        return {k:v for k,v in row.items()if k not in ('status','program_bytes','program_blake3','active_state')}
    assert shared(closed)==shared(provisional),'preclosure metadata changed'
    state=dict(closed['active_state'])
    coverage=state.pop('historical_candidate_collection')
    assert state==provisional['active_state'],'preclosure active state changed'
    assert state['history_is_coverage_certificate'] is False
    history=state['historical_requests'];capacity=closed['storage_capacity']
    assert all(isinstance(label,list)and len(label)==capacity and all(type(n)is int for n in label)for label in history)
    assert history==sorted(history) and len({tuple(x)for x in history})==len(history)
    count=len(history)if expected_limit is None else min(expected_limit,len(history))
    assert coverage['max_history_candidate_maps']==expected_limit
    assert coverage['selection']=='ascending stored-index lexicographic prefix'
    assert coverage['selected']==history[:count]
    assert coverage['unassessed']==history[count:]
    for key in ('assessed_count','certified_map_count','assessed_unresolved_count'):
        assert type(coverage[key])is int and coverage[key]>=0
    assert coverage['assessed_count']==count
    assert coverage['certified_map_count']+coverage['assessed_unresolved_count']==count
    assert coverage['unassessed_are_unresolved'] is False
    assert coverage['history_is_coverage_certificate'] is False
    assert len(retired)==coverage['assessed_unresolved_count']
    retired_labels=[entry['requested']for entry in retired]
    assert len({tuple(x)for x in retired_labels})==len(retired_labels)
    assert all(label in coverage['selected']for label in retired_labels)
    assert all(isinstance(entry['unresolved'],list)and entry['unresolved']for entry in retired)
    return {'max_history_candidate_maps':expected_limit,'history_count':len(history),
        'selected_count':count,'assessed_count':coverage['assessed_count'],
        'certified_map_count':coverage['certified_map_count'],
        'assessed_unresolved_count':coverage['assessed_unresolved_count'],
        'unassessed_count':len(coverage['unassessed']),
        'shared_preclosure_metadata_and_active_state_identical':True,
        'unassessed_implies_no_value_or_failure_classification':True}
