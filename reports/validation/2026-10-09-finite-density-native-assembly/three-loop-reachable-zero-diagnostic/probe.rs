#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::{Decode,Encode};
use rustred::persistence::{inspect_program, BinaryIoLimits, SectionTag, DecodedCoefficientTable, CoefficientId};
use rustred::solver::guarded::{GuardedSource, GuardedSourceSystem, GuardedProgram, IndexRole, IndexBounds, IndexDomain};
use rustred::solver::{Integral, Power, Term};
use std::sync::Arc;
type Domain=Vec<(Option<i64>,Option<i64>)>;
type Label=Vec<(bool,i16)>;
#[derive(Decode,Encode)] struct Tr { integral: Label, coefficient:usize }
#[derive(Decode,Encode)] struct Sr { id:String, domain:Domain, conditions:Vec<usize>, terms:Vec<Tr> }
#[derive(Decode,Encode)] struct Sd { row:usize, integral:Label, shifts:Vec<i16> }
#[derive(Decode,Encode)] struct Rr { fixed:Vec<Option<i16>>,target:Label,rhs:Vec<Tr>,sources:Vec<Sd>,domain:Domain,discovery_domain:Domain,conditions:Vec<usize>,sector:Vec<bool>,permutation:Option<Vec<usize>> }
#[derive(Decode,Encode)] struct Rec { schema:String,measure:String,roles:Vec<u8>,indices:Vec<usize>,sources:Vec<Sr>,zero_domains:Vec<Domain>,rules:Vec<Rr>,terminals:Vec<Vec<i64>> }
fn domain<const N:usize>(v:&Domain)->IndexDomain<N>{IndexDomain::new(v.iter().map(|&(lo,hi)|IndexBounds::new(lo,hi).unwrap()).collect::<Vec<_>>().try_into().unwrap()).unwrap()}
fn label<const N:usize>(v:&Label)->Integral<N>{Integral::new(v.iter().map(|&(s,p)|Power::new(s,p).unwrap()).collect::<Vec<_>>().try_into().unwrap())}
fn bounds<const N:usize>(v:&IndexDomain<N>)->Domain{v.bounds().iter().map(|b|(b.lower(),b.upper())).collect()}

fn main(){
 let a=std::env::args().collect::<Vec<_>>();let output=std::path::Path::new(&a[2]);std::fs::create_dir_all(output).unwrap();
 let bytes=std::fs::read(&a[1]).unwrap();let limits=BinaryIoLimits::default();let env=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(env.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,env.section(SectionTag::PROGRAM).unwrap().len());
 let table=DecodedCoefficientTable::import_generated_normalized(env.section(SectionTag::SYMBOLICA_STATE).unwrap(),env.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let original=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();assert!(original.terminals().is_empty());

 let targets:Vec<[i64;12]>=serde_json::from_slice(&std::fs::read(&a[3]).unwrap()).unwrap();
 let targets=targets.into_iter().take(15).collect::<Vec<_>>();
 let history:serde_json::Value=serde_json::from_slice(&std::fs::read(&a[4]).unwrap()).unwrap();
 let history=history["historical_requests"].as_array().unwrap().iter().map(|v|serde_json::from_value::<[i64;12]>(v.clone()).unwrap()).collect::<std::collections::BTreeSet<_>>();
 let mut pending=targets.iter().copied().collect::<std::collections::BTreeSet<_>>();
 let mut visited=std::collections::BTreeSet::new();
 let mut applications=Vec::new();
 const NODE_LIMIT:usize=10000;
 while let Some(point)=pending.pop_first(){
  if !visited.insert(point){continue}
  let applied=original.apply(&point).unwrap();
  for child in applied.terms.keys(){if !visited.contains(child){pending.insert(*child);}}
  applications.push(serde_json::json!({"point":point.to_vec(),"status":format!("{:?}",applied.status),"children":applied.terms.len(),"previously_requested":history.contains(&point)}));
  if visited.len()==NODE_LIMIT{break}
 }
 let hidden=visited.difference(&history).copied().collect::<Vec<_>>();
 let zero=sources.direct_zero_rules_at_points(&hidden,2_000_000).unwrap();
 let new_zero_points=zero.solution.rules.iter().map(|r|r.candidate().target.powers().iter().map(|p|{assert!(!p.is_symbolic());i64::from(p.value())}).collect::<Vec<_>>()).collect::<Vec<_>>();
 let rules_count=zero.solution.rules.len();
 let completed_points=zero.completed_points.len();
 let zero_program=GuardedProgram::new(sources.clone(),zero.solution.rules,[]).unwrap();
 let combined=zero_program.union_replayed(GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap(),[],65536).unwrap();
 std::fs::write(output.join("zero-first.bin"),combined.encode_native(limits).unwrap()).unwrap();
 let audit=|program:&GuardedProgram<12>|->serde_json::Value{serde_json::Value::Array(targets.iter().map(|target|{
  let result=program.reduce(*target,Default::default()).unwrap();
  serde_json::json!({"target":target.to_vec(),"rule_applications":result.rule_applications,"terms":result.terms.iter().map(|(i,c)|serde_json::json!({"point":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"unresolved":result.unresolved.iter().map(|r|serde_json::json!({"point":r.integral.to_vec(),"coefficient":r.coefficient.to_expression().to_string(),"reason":format!("{:?}",r.reason)})).collect::<Vec<_>>(),"conditions":result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})
 }).collect())};
 let report=serde_json::json!({"scope":"structural reachable-node audit only; traversal is an overestimate before weighted cancellation; all reduction and zero proofs owned by native RustRed with full original-source replay; no amplitude or boundary evaluated","node_limit":NODE_LIMIT,"traversal_complete":pending.is_empty(),"visited_nodes":visited.len(),"hidden_nodes":hidden.len(),"zero_attempt_budget":2_000_000,"attempted_rows":zero.attempted_rows,"completed_points":completed_points,"new_zero_rules":rules_count,"new_zero_points":new_zero_points,"discovery_gaps":zero.solution.unresolved.iter().map(|r|serde_json::json!({"domain":bounds(&r.domain),"reason":format!("{:?}",r.reason),"detail":r.detail})).collect::<Vec<_>>(),"skipped_alignment_count":zero.skipped_seeds.len(),"applications":applications,"original":audit(&original),"zero_first":audit(&combined)});
 std::fs::write(output.join("result.json"),serde_json::to_vec_pretty(&report).unwrap()).unwrap();
 println!("Reachable nodes {}, hidden {}, new exact zero rules {}",visited.len(),hidden.len(),rules_count);
}
