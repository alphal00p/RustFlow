#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::Decode;
use rustred::persistence::{inspect_program, BinaryIoLimits, SectionTag, DecodedCoefficientTable, CoefficientId};
use rustred::solver::guarded::{GuardedSource, GuardedSourceSystem, GuardedProgram, IndexRole, IndexBounds, IndexDomain};
use rustred::solver::{Integral, Power, Term, SearchOptions};
use std::sync::Arc;
type Domain=Vec<(Option<i64>,Option<i64>)>;
type Label=Vec<(bool,i16)>;
#[derive(Decode)] struct Tr { integral: Label, coefficient:usize }
#[derive(Decode)] struct Sr { id:String, domain:Domain, conditions:Vec<usize>, terms:Vec<Tr> }
#[derive(Decode)] struct Sd { row:usize, integral:Label, shifts:Vec<i16> }
#[derive(Decode)] struct Rr { fixed:Vec<Option<i16>>,target:Label,rhs:Vec<Tr>,sources:Vec<Sd>,domain:Domain,discovery_domain:Domain,conditions:Vec<usize>,sector:Vec<bool>,permutation:Option<Vec<usize>> }
#[derive(Decode)] struct Rec { schema:String,measure:String,roles:Vec<u8>,indices:Vec<usize>,sources:Vec<Sr>,zero_domains:Vec<Domain>,rules:Vec<Rr>,terminals:Vec<Vec<i64>> }
fn domain<const N:usize>(v:&Domain)->IndexDomain<N>{IndexDomain::new(v.iter().map(|&(lo,hi)|IndexBounds::new(lo,hi).unwrap()).collect::<Vec<_>>().try_into().unwrap()).unwrap()}
fn label<const N:usize>(v:&Label)->Integral<N>{Integral::new(v.iter().map(|&(s,p)|Power::new(s,p).unwrap()).collect::<Vec<_>>().try_into().unwrap())}
fn bounds<const N:usize>(v:&IndexDomain<N>)->Domain{v.bounds().iter().map(|b|(b.lower(),b.upper())).collect()}

use symbolica::prelude::*;
use std::collections::BTreeMap;
fn main() {
 let args:Vec<_>=std::env::args().collect();
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let (record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;16]=record.roles.iter().map(|r|match r {0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let indices:[usize;16]=record.indices.clone().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,indices,rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let started=std::time::Instant::now();
 let historical=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap().with_terminals_replayed([],100_000).unwrap();
 let historical_replay_seconds=started.elapsed().as_secs_f64();
 let points:Vec<[i64;16]>=[2_i64,3].into_iter().map(|power|{let mut point=[0_i64;16];point[..6].copy_from_slice(&[1,power,0,1,0,-2]);point}).collect();
 let ordinary=sources.solve_domains(points.iter().map(|p|IndexDomain::new(p.map(IndexBounds::fixed)).unwrap()).collect(),SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},points.len()).unwrap();
 assert!(ordinary.unresolved.is_empty());
 let historical=GuardedProgram::new(sources.clone(),ordinary.rules,[]).unwrap().union_replayed(historical,[],100_000).unwrap();
 let alignments_per_point=record.sources.iter().map(|s|s.terms.len()).sum::<usize>();
 let maximum=alignments_per_point.checked_mul(points.len()).unwrap();
 let mut per_point=Vec::new();
 let mut minimum=0_usize;
 for point in &points {
  let started=std::time::Instant::now();
  let found=sources.direct_zero_rules_at_points(&[*point],alignments_per_point).unwrap();
  assert_eq!(found.solution.rules.len(),1,"actual even target needs one-row proof");
  assert!(found.solution.unresolved.is_empty());assert_eq!(found.completed_points,vec![*point]);
  minimum=minimum.checked_add(found.attempted_rows).unwrap();
  let rule=&found.solution.rules[0];
  assert!(rule.candidate().rhs.is_empty());sources.replay_rule(rule).unwrap();
  let traces=rule.candidate().sources.iter().map(|s|serde_json::json!({"source_ordinal":s.basis_row,"source_id":record.sources[s.basis_row].id,"seed":(0..16).map(|i|s.seed.integral[i].value()).collect::<Vec<_>>(),"shifts":s.seed.shifts.to_vec()})).collect::<Vec<_>>();
  per_point.push(serde_json::json!({"target":point.to_vec(),"global_attempt_budget":alignments_per_point,"attempted_rows":found.attempted_rows,"seconds":started.elapsed().as_secs_f64(),"traces":traces,"nonzero_conditions":rule.nonzero_conditions().iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"skipped_seeds":found.skipped_seeds.iter().map(|s|format!("{s:?}")).collect::<Vec<_>>()}));
 }
 let started=std::time::Instant::now();
 let found=sources.direct_zero_rules_at_points(&points,minimum).unwrap();
 let combined_seconds=started.elapsed().as_secs_f64();
 assert_eq!(found.attempted_rows,minimum);assert_eq!(found.completed_points,points);assert_eq!(found.solution.rules.len(),2);assert!(found.solution.unresolved.is_empty());
 let budget=sources.direct_zero_rules_at_points(&points,minimum-1).unwrap();
 assert_eq!(budget.solution.rules.len(),1);assert_eq!(budget.completed_points,vec![points[0]]);assert_eq!(budget.solution.unresolved.len(),1);
 assert!(matches!(budget.solution.unresolved[0].reason,rustred::solver::guarded::GuardedUnresolvedReason::DomainBudget));
 let zero=GuardedProgram::new(sources.clone(),found.solution.rules,[]).unwrap();
 let zero_bytes=zero.encode_native(limits).unwrap();
 let zero=GuardedProgram::decode_generated(&zero_bytes,sources.clone(),limits).unwrap();
 let mut before=Vec::new();
 for point in &points {
  let applied=historical.apply(point).unwrap();
  assert!(!applied.terms.is_empty(),"control must retain the old nonzero recurrence");
  before.push(serde_json::json!({"target":point.to_vec(),"status":format!("{:?}",applied.status),"terms":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"nonzero_conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()}));
 }
 let started=std::time::Instant::now();
 let union=zero.union_replayed(historical,[],100_000).unwrap();
 let union_replay_seconds=started.elapsed().as_secs_f64();
 let mut after=Vec::new();
 for point in &points {
  let applied=union.apply(point).unwrap();
  assert!(matches!(applied.status,rustred::solver::guarded::GuardedApplicationStatus::Applied{..}));
  assert!(applied.terms.is_empty());
  let reduced=union.reduce(*point,Default::default()).unwrap();assert!(reduced.unresolved.is_empty());assert!(reduced.terms.is_empty());
  after.push(serde_json::json!({"target":point.to_vec(),"status":format!("{:?}",applied.status),"terms":[],"rule_applications":reduced.rule_applications,"nonzero_conditions":reduced.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()}));
 }
 let result=serde_json::json!({"status":"pass","scope":"Exact native direct-one-original-row zero discovery on the unchanged historical three-loop source corpus. No analytic zero, no custom row elimination, no amplitude value.","source_program":args[1],"source_count":record.sources.len(),"historical_rules":record.rules.len(),"historical_replay_seconds":historical_replay_seconds,"maximum_all_alignments_for_two_points":maximum,"minimum_global_attempts_for_both_in_this_order":minimum,"combined_search_seconds":combined_seconds,"per_point":per_point,"one_less_budget":{"attempts":budget.attempted_rows,"proved_rules":budget.solution.rules.len(),"completed_points":budget.completed_points.iter().map(|p|p.to_vec()).collect::<Vec<_>>(),"gaps":budget.solution.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>()},"zero_program_persistence_replayed":true,"old_recurrences":before,"zero_first_union":after,"union_rules":union.rules().len(),"union_replay_seconds":union_replay_seconds,"analytic_values_read":0,"oracle_answers_read":0});
 std::fs::write(&args[2],serde_json::to_vec_pretty(&result).unwrap()).unwrap();
 println!("{}",serde_json::json!({"status":"pass","points":points.len(),"attempts":minimum,"search_seconds":combined_seconds,"replay_seconds":historical_replay_seconds+union_replay_seconds,"output":args[2]}));
}
