#!/usr/bin/env python3
"""Read-only query reconstruction: no native search, replay or reduction."""
from pathlib import Path
import collections, gzip, hashlib, json
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
FH=ROOT/'reports/validation/2026-10-10-final-history-work-profile'
FP=ROOT/'reports/validation/2026-10-10-final-reduction-phase-diagnostic'
LL=ROOT/'reports/validation/2026-10-10-longer-singleton-control'
bindings={}
def read(path):
    raw=path.read_bytes(); bindings[str(path.relative_to(ROOT))]={'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)
def derivative(label):
    return [tuple(list(label[:a])+[label[a]+1]+list(label[a+1:])) for a in SHIFTED if label[a]!=0]
SHIFTED=[0,1,2,4]
final=read(FH/'input.json'); selected=read(FH/'points.json'); observations=read(FP/'profile/terminal-observations.json')
profile=read(FH/'profile/result.json'); phase=read(FP/'summary.json'); history_scope=read(LL/'summary.json')
assert 'native Di-eta on [0, 1, 2, 4]' in final['measure_id']
frontier=[tuple(x)for x in final['frontier']]
assert [x['point']for x in observations]==final['frontier']
assert all(x['status']=='Unresolved(NoApplicableRule)' and x['terms']==0 for x in observations)
needed=set(map(tuple,final['active_state']['next_needed']))
roots=sorted(needed-set(frontier)-{d for f in frontier for d in derivative(f)})
assert len(roots)==3
# Independently recover the three roots from the initial request + first-derivative union.
physical_input=read(LL/'single3/input.json')
assert [(t['powers'],t['numerator'])for t in physical_input['targets']]==[([1,1,1,1,1],'1'),([2,1,1,1,1],'g1_2+u1*u2')]
initial=read(LL/'single3/native-closure/round-000-provisional.json')
initial_requests=set(map(tuple,initial['active_state']['submitted_requests']))
minimal=[]
for x in initial_requests:
    has_parent=False
    for a in SHIFTED:
        y=list(x); y[a]-=1
        if y[a]!=0 and tuple(y) in initial_requests: has_parent=True
    if not has_parent:minimal.append(x)
assert sorted(minimal)==roots
assert initial_requests==set(roots)|{d for r in roots for d in derivative(r)}
rounds=[]
for path in sorted((LL/'single3/native-closure').glob('round-*-provisional.json*')):
    m=read(path); h=set(map(tuple,m['active_state']['historical_requests'])); calls=[{'kind':'original-target-component','point':list(x)}for x in roots]; skipped=[]
    for index,f in enumerate(m['frontier']):
        f=tuple(f); ds=derivative(f)
        if f not in h or any(d not in h for d in ds): skipped.append(list(f));continue
        calls.extend({'kind':'basis-derivative-component','basis_index':index,'point':list(d)}for d in ds)
    counts=collections.Counter(tuple(c['point'])for c in calls)
    rounds.append({'round':m['round'],'frontier':len(m['frontier']),'program_blake3':m['program_blake3'],'rules':m['native_rule_count'],'reconstructed_calls':len(calls),'distinct_calls':len(counts),'duplicate_calls':len(calls)-len(counts),'deferred_basis_labels':len(skipped),'calls':calls})
assert len(rounds)==19 and rounds[-1]['round']==18
assert not rounds[-1]['deferred_basis_labels'] and not final['active_state']['deferred_requests']
assert len(set(r['program_blake3']for r in rounds))==19
final_calls=[{'kind':'original-target-component','point':list(x)}for x in roots]
final_calls.extend({'kind':'basis-derivative-component','basis_index':i,'point':list(d)}for i,f in enumerate(frontier)for d in derivative(f))
assert final_calls==rounds[-1]['calls']
points=[tuple(x['point'])for x in selected['points']]
assert len(points)==len(set(points))==42
assert len(final_calls)==len({tuple(c['point'])for c in final_calls})==19
by_point={tuple(p['selection']['point']):p for p in profile['points']}
assert all(tuple(c['point'])in by_point for c in final_calls)
for filename in ['src/finite_density/reduction/active.rs','src/finite_density/reduction.rs','src/finite_density/guarded.rs','vendor/rustred/crates/rustred-core/src/solver/guarded/lifecycle.rs',str((FH/'probe.rs').relative_to(ROOT)),str((FP/'probe.rs').relative_to(ROOT))]:
    p=ROOT/filename; raw=p.read_bytes(); bindings[filename]={'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
result={'scope':'Read-only code/metadata reconstruction, not a newly executed native trace or full workflow timing. Root multiplicity is checked against initial request construction and the actual two-target fixture. No production edits or native run.',
 'roots':[list(x)for x in roots],'shifted_slots':SHIFTED,'selected_profile_calls':42,'selected_profile_distinct_calls':42,
 'rounds':rounds,'round_total_calls':sum(r['reconstructed_calls']for r in rounds),'round_total_duplicate_calls':sum(r['duplicate_calls']for r in rounds),'distinct_round_programs':19,
 'final_provisional_calls':rounds[-1]['calls'],'scheduled_final_calls':final_calls,'scheduled_final_exact_repeats':19,
 'final_call_execution_caveat':'All nineteen provisional calls precede the published round18 checkpoint. Source control flow schedules the same nineteen calls after final terminal rebinding. The original run timed out without a closed checkpoint, so this does NOT establish that all final repeated calls executed or completed.',
 'final_call_binding_change':'Same ordered native rules/source/order and application limits; empty old terminals become the seven explicitly verified NoApplicableRule labels. Cross-terminal cache reuse is not automatic.',
 'all_added_terminals_previously_no_rule':True,'historical_original_run_status':history_scope['status'],
 'profile_slow_point_seconds':313.80065236,'profile_native_apply_seconds':phase['operation_seconds']['native-apply'],
 'profile_native_apply_fraction_instrumented_wall':phase['operation_seconds']['native-apply']/phase['instrumented_reduction_wall_seconds'],
 'descendant_overlap':'Unavailable: FH saves final maps/counts; FP saves operation totals and overwrites a last-operation progress file, not the full per-label trace. Common output terminals are not evidence of overlapping expensive reduced tails.',
 'recommendation':'Do not substitute generic fully reduced unit tails. Investigate only bounded whole-call reuse under an identical immutable binding, optionally with the separately checked NoRule-only terminal-promotion equivalence. Keep final independent rule replay and weighted target/derivative audit.'}
(OUT/'observations.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'input-bindings.json').write_text(json.dumps({'files':bindings},indent=2)+'\n')
print(json.dumps({k:result[k]for k in ['round_total_calls','round_total_duplicate_calls','distinct_round_programs','scheduled_final_exact_repeats','profile_native_apply_fraction_instrumented_wall']}))
