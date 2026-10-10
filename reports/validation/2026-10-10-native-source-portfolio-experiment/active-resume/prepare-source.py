#!/usr/bin/env python3
"""Derive a timeout-resume wrapper without changing the native closure algorithm."""
import hashlib
import json
from pathlib import Path

base = Path(__file__).resolve().parent
original = base.parent / 'active-pilot/pilot.rs'
output = base / 'pilot.rs'
assert not output.exists()
source = original.read_text()

def replace(old, new):
    global source
    assert source.count(old) == 1, old[:100]
    source = source.replace(old, new)

replace('let mut portfolio_stats=Vec::new();', 'let mut portfolio_stats=Vec::<serde_json::Value>::new();')
replace('let mut history=Vec::new();', 'let mut history=Vec::<serde_json::Value>::new();')
replace('    for rule in exact.rules() { if let Some(s)=rule.candidate().stats.guarded_portfolio {',
        '    #[cfg(source_portfolio)]\n    for rule in exact.rules() { if let Some(s)=rule.candidate().stats.guarded_portfolio {')
replace('let old_proof_replay=if std::env::var("PILOT_REPLAY_OLD").as_deref()==Ok("true"){Some(GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap().rules().len())}else{None};',
        'assert_ne!(std::env::var("PILOT_REPLAY_OLD").as_deref(),Ok("true"));\n let old_proof_replay:Option<usize>=None;')
replace(' std::fs::write(out.join("configuration.json"),serde_json::to_vec_pretty(&serde_json::json!(',
        ' let configuration=serde_json::json!(')
replace('"measure_id":sources.measure_id()})).unwrap()).unwrap();',
        '"measure_id":sources.measure_id()});\n std::fs::write(out.join("configuration.json"),serde_json::to_vec_pretty(&configuration).unwrap()).unwrap();')
needle=' let mut attempted_points=BTreeSet::new();\n for round in 0..max_rounds {'
replacement=r''' let mut attempted_points=BTreeSet::new();
 let mut first_round=0;
 let mut resume_evidence=serde_json::Value::Null;
 if let Some(parent_arg)=args.get(3) {
  let parent=std::path::Path::new(parent_arg);
  let read=|name:&str| -> serde_json::Value {serde_json::from_slice(&std::fs::read(parent.join(name)).unwrap()).unwrap()};
  let resource=read("resources.json");let run_binding=read("run-binding.json");
  assert_eq!(resource["exit_code"],124);assert_eq!(run_binding["exit_code"],124);
  assert_eq!(run_binding["post_run_inputs_unchanged"],true);
  assert_eq!(run_binding["engine"],if cfg!(source_portfolio){"proposal"}else{"baseline"});
  assert!(!parent.join("result.json").exists()&&!parent.join("closed.bin").exists(),"only unfinished timeout checkpoints may resume");
  assert_eq!(read("configuration.json"),configuration,"source identity, targets and every limit must match");
  assert_eq!((max_rounds,max_requested,max_frontier,max_rules),(16,4096,1024,65536));
  assert_eq!((zero_attempts,depth,ray_budget,point_budget),(8192,3,1,1));
  history=serde_json::from_value(read("progress.json")).unwrap();
  assert!(!history.is_empty());
  for (i,item) in history.iter().enumerate(){assert_eq!(item["round"].as_u64(),Some(i as u64));}
  let last=history.last().unwrap();let last_round=last["round"].as_u64().unwrap()as usize;
  first_round=last_round+1;assert!(first_round<max_rounds,"global round budget already exhausted");
  let round_file=format!("round-{last_round:03}.json");
  let proof_file=format!("round-{last_round:03}.bin");
  let state=read(&round_file);assert_eq!(state["round"].as_u64(),Some(last_round as u64));
  assert_eq!(state["frontier_budget_exhausted"],false);
  let points=|field:&str| -> BTreeSet<Point> {
   let raw:Vec<Point>=serde_json::from_value(state[field].clone()).unwrap();
   let set=raw.iter().copied().collect::<BTreeSet<_>>();assert_eq!(set.len(),raw.len());
   assert!(set.iter().all(|i|admitted.contains(i)));set
  };
  frontier=points("frontier");attempted_points=points("attempted_exact_points");requested=points("next_needed");
  assert!(roots.is_subset(&requested)&&frontier.is_subset(&requested));
  assert!(points("requested").is_subset(&attempted_points));
  assert!(!points("pending_new_derivative_labels").is_empty(),"a completed closing round must not resume");
  assert!(frontier.len()<=max_frontier&&attempted_points.union(&requested).count()<=max_requested);
  all_conditions=serde_json::from_value(state["conditions_encountered"].clone()).unwrap();
  program=GuardedProgram::decode_generated(&std::fs::read(parent.join(&proof_file)).unwrap(),sources.clone(),limits).unwrap();
  assert!(program.terminals().is_empty(),"a resumable program cannot contain declared master terminals");
  assert!(program.rules().len()<=max_rules);
  assert_eq!(state["rule_count"].as_u64(),Some(program.rules().len()as u64));
  assert_eq!(last["rules"],state["rule_count"]);assert_eq!(last["frontier"].as_u64(),Some(frontier.len()as u64));
  assert_eq!(last["requested"].as_u64(),Some(points("requested").len()as u64));
  let mut replayed_targets=Vec::new();
  for row in&targets {let replay=reduce(&program,row,&admitted);assert!(replay.failures.is_empty());assert!(replay.conditions.is_subset(&all_conditions));replayed_targets.push(display(&replay.leaves));}
  assert_eq!(serde_json::json!(replayed_targets),state["original_target_rows"]);
  for row in state["derivative_rows"].as_array().unwrap(){
   let basis:Point=serde_json::from_value(row["basis"].clone()).unwrap();assert!(frontier.contains(&basis));
   let raw=derivative(basis,&shifted);assert!(raw.keys().all(|i|attempted_points.contains(i)));
   assert_eq!(serde_json::json!(display(&raw)),row["raw_derivative"]);
   let replay=reduce(&program,&raw,&admitted);assert!(replay.failures.is_empty());assert!(replay.conditions.is_subset(&all_conditions));
   assert_eq!(serde_json::json!(display(&replay.leaves)),row["reduced_actual_derivative"]);
  }
  resume_evidence=serde_json::json!({"parent":parent_arg,"last_completed_round":last_round,
   "first_continuation_round":first_round,"restored_rules":program.rules().len(),
   "restored_frontier":frontier.len(),"restored_attempted_exact_points":attempted_points.len(),
   "restored_requested":requested.len(),"restored_conditions":all_conditions,
   "native_decode_and_saved_target_derivative_replay":true,
   "original_corpus_rules_imported":false,"only_this_engine_parent_discovered_rules_restored":true});
  std::fs::write(out.join("resume-state.json"),serde_json::to_vec_pretty(&resume_evidence).unwrap()).unwrap();
  std::fs::write(out.join("progress.json"),serde_json::to_vec_pretty(&history).unwrap()).unwrap();
 }
 for round in first_round..max_rounds {'''
replace(needle,replacement)
replace('let result=serde_json::json!({"status":stop,"closed":closed,',
        'let result=serde_json::json!({"status":stop,"closed":closed,"resume":resume_evidence,')
output.write_text(source)
(base/'source-binding.json').write_text(json.dumps({
    'scope':'Validation-only checkpoint resume scaffolding and cfg-gated statistics; the native search/reduction/closure algorithm is retained.',
    'base_source':str(original), 'base_source_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),
    'generated_source_sha256_before_formatting':hashlib.sha256(output.read_bytes()).hexdigest(),
    'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
},indent=2)+'\n')
