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

use std::collections::BTreeSet;
use std::time::Instant;
fn main() {
 let args=std::env::args().collect::<Vec<_>>();let output=std::path::Path::new(&args[3]);std::fs::create_dir_all(output).unwrap();
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let (record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;16]=record.roles.iter().map(|r|match r {0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let previous=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();let previous_rules=previous.rules().len();assert!(previous.terminals().is_empty(),"prototype requires the saved provisional program with no stopping labels");
 let metadata:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();assert_eq!(metadata["physical_arity"],13);assert_eq!(metadata["storage_capacity"],16);
 let frontier:Vec<[i64;16]>=serde_json::from_value(metadata["frontier"].clone()).unwrap();let requested:Vec<[i64;16]>=serde_json::from_value(metadata["requested"].clone()).unwrap();
 let setting=|key:&str,default:usize|std::env::var(key).ok().map(|s|s.parse::<usize>().unwrap()).unwrap_or(default);
 let per=setting("PROBE_PER_RAY",32);let total=setting("PROBE_TOTAL",8192);let depth=setting("PROBE_DEPTH",3) as u32;assert!(per>0&&total>0);
 let started=Instant::now();let mut remaining=total;let mut rules=Vec::new();let mut ray_reports=Vec::new();let mut omitted=Vec::new();let mut non_no_rule=Vec::new();
 for target in &frontier {
  assert!(target[13..].iter().all(|&n|n==0));
  let status=previous.apply(target).unwrap().status;
  if !matches!(status,rustred::solver::guarded::GuardedApplicationStatus::Unresolved(rustred::solver::guarded::GuardedApplicationFailure::NoApplicableRule)) {
   non_no_rule.push(serde_json::json!({"target":target.to_vec(),"status":format!("{status:?}")}));continue;
  }
  let ray=IndexDomain::new(std::array::from_fn(|i|if target[i]==0 {IndexBounds::fixed(0)}else if target[i]>0 {IndexBounds::new(Some(target[i]),None).unwrap()}else {IndexBounds::new(None,Some(target[i])).unwrap()})).unwrap();
  if remaining==0 {omitted.push(serde_json::json!({"target":target.to_vec(),"domain":bounds(&ray),"reason":"not submitted: global allocated-domain budget exhausted"}));continue;}
  let allocated=remaining.min(per);remaining-=allocated;let found=sources.solve_domains(vec![ray.clone()],SearchOptions{max_depth:Some(depth),sample_seed:0,..Default::default()},allocated).unwrap();
  ray_reports.push(serde_json::json!({"target":target.to_vec(),"domain":bounds(&ray),"allocated_domains":allocated,"rules":found.rules.len(),"unresolved":found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>()}));rules.extend(found.rules);
  if ray_reports.len()%32==0 {eprintln!("searched {} residual rays, {} allocated domains, {} new native rules",ray_reports.len(),total-remaining,rules.len());}
 }
 let discovery_seconds=started.elapsed().as_secs_f64();let new_rules=rules.len();
 std::fs::write(output.join("discovery.json"),serde_json::to_vec_pretty(&serde_json::json!({"scope":"bounded symbolic ray search; allocated visits are conservative upper bounds, not measured native work","previous_rules":previous_rules,"input_frontier":frontier.len(),"original_requested":requested.len(),"per_ray":per,"total_allocated_cap":total,"allocated":total-remaining,"discovery_seconds":discovery_seconds,"new_rules":new_rules,"rays":ray_reports,"unsubmitted":omitted,"non_no_rule_input_statuses":non_no_rule})).unwrap()).unwrap();
 let replay_started=Instant::now();let batch=GuardedProgram::new(sources.clone(),rules,[]).unwrap();
 let mut batch_applies=0;for item in &ray_reports {let target:[i64;16]=serde_json::from_value(item["target"].clone()).unwrap();if matches!(batch.apply(&target).unwrap().status,rustred::solver::guarded::GuardedApplicationStatus::Applied{..}) {batch_applies+=1;}}
 let combined=previous.union_replayed(batch,[],65536).unwrap();let replay_seconds=replay_started.elapsed().as_secs_f64();
 std::fs::write(output.join("combined.bin"),combined.encode_native(limits).unwrap()).unwrap();eprintln!("native batch replay+old-first union: {}+{} rules; {} searched targets now directly apply; auditing {} original requests",previous_rules,new_rules,batch_applies,requested.len());
 let mut unique=BTreeSet::new();let mut conditions=BTreeSet::new();let mut records=Vec::new();let mut unresolved_count=0;let mut applications=0;let mut other_failures=0;
 let audit_started=Instant::now();
 for target in &requested {
  let reduced=combined.reduce(*target,Default::default()).unwrap();applications+=reduced.rule_applications;unresolved_count+=reduced.unresolved.len();
  for c in reduced.nonzero_conditions {conditions.insert(c.to_expression().to_string());}
  for term in &reduced.unresolved {if term.reason==rustred::solver::guarded::GuardedApplicationFailure::NoApplicableRule {unique.insert(term.integral);}else{other_failures+=1;}}
  records.push(serde_json::json!({"requested":target.to_vec(),"terminal_terms":reduced.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"unresolved":reduced.unresolved.iter().map(|t|serde_json::json!({"indices":t.integral.to_vec(),"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>()}));
  if records.len()%128==0 {eprintln!("audited {}/{} requests, {} unique residual labels",records.len(),requested.len(),unique.len());std::fs::write(output.join("progress.json"),serde_json::to_vec_pretty(&serde_json::json!({"audited":records.len(),"requested":requested.len(),"unique_no_rule_labels":unique.len(),"other_failures":other_failures,"applications":applications,"audit_seconds":audit_started.elapsed().as_secs_f64(),"complete":false})).unwrap()).unwrap();}
 }
 let result=serde_json::json!({"scope":"source replay and old-first retained rule union; all original requests reduced independently; no master, closure, boundary or numerical acceptance claim","input_frontier":frontier.len(),"original_requested":requested.len(),"searched_rays":ray_reports.len(),"allocated_domains":total-remaining,"previous_rules":previous_rules,"new_rules":new_rules,"combined_rules":combined.rules().len(),"searched_targets_directly_applied":batch_applies,"discovery_seconds":discovery_seconds,"replay_union_seconds":replay_seconds,"audit_seconds":audit_started.elapsed().as_secs_f64(),"unique_no_rule_labels":unique.len(),"unresolved_terms":unresolved_count,"non_no_rule_failures":other_failures,"rule_applications":applications,"remaining_labels":unique.into_iter().map(|x|x.to_vec()).collect::<Vec<_>>(),"conditions":conditions,"reductions":records,"complete":true});
 std::fs::write(output.join("result.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();eprintln!("Result: {} -> {} unique residuals, {} unresolved terms, {} other failures",frontier.len(),result["unique_no_rule_labels"],unresolved_count,other_failures);
}
