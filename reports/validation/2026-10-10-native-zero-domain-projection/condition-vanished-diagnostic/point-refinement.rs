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
fn main(){
 let args:Vec<_>=std::env::args().collect();let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v2");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;9]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let p=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();let point=[1,1,2,-1,-1,0,0,0,0];let applied=p.apply(&point).unwrap();
 let rules=p.rules().iter().enumerate().filter(|(_,r)|r.domain().contains(&point)).map(|(n,r)|{let c=r.candidate();serde_json::json!({"ordinal":n,"domain":bounds(r.domain()),"target":c.target.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"conditions":r.nonzero_conditions().iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"rhs":c.rhs.iter().map(|t|serde_json::json!({"integral":t.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"coefficient":t.coefficient.to_expression().to_string()})).collect::<Vec<_>>(),"proof":c.sources.iter().map(|s|{let sr=&record.sources[s.basis_row];serde_json::json!({"ordinal":s.basis_row,"source_id":sr.id,"source_domain":sr.domain,"source_conditions":sr.conditions.iter().map(|&id|polynomial(id).to_expression().to_string()).collect::<Vec<_>>(),"source_terms":sr.terms.iter().map(|t|serde_json::json!({"integral":t.integral,"coefficient":polynomial(t.coefficient).to_expression().to_string()})).collect::<Vec<_>>(),"seed":s.seed.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"seed_shifts":s.seed.shifts.to_vec()})}).collect::<Vec<_>>()})}).collect::<Vec<_>>();
 let found=sources.solve_domains(vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},32).unwrap();
 let gaps=found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
 let fresh=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();let fresh_result=fresh.apply(&point).unwrap();
 let fresh_report=serde_json::json!({"rule_count":fresh.rules().len(),"status":format!("{:?}",fresh_result.status),"conditions":fresh_result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"rhs":fresh_result.terms.iter().map(|(i,c)|serde_json::json!({"integral":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"gaps":gaps});
 let union=p.union_replayed(fresh,[],65536).unwrap();let union_result=union.apply(&point).unwrap();let reduced=union.reduce(point,Default::default()).unwrap();
 let union_report=serde_json::json!({"precedence":"retained-first","rule_count":union.rules().len(),"status":format!("{:?}",union_result.status),"conditions":union_result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"rhs":union_result.terms.iter().map(|(i,c)|serde_json::json!({"integral":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"reduction_applications":reduced.rule_applications,"reduction_conditions":reduced.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"reduction_terms":reduced.terms.iter().map(|(i,c)|serde_json::json!({"integral":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"reduction_unresolved":reduced.unresolved.iter().map(|t|serde_json::json!({"integral":t.integral,"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>()});
 std::fs::write(std::path::Path::new(&args[2]).with_extension("bin"),union.encode_native(limits).unwrap()).unwrap();
 let out=serde_json::json!({"point_search":{"max_depth":3,"max_domains":32,"result":fresh_report},"retained_first_union":union_report,"scope":"exact fresh native v2 decode/replay; rule inspection only","program":args[1],"measure_id":record.measure,"point":point,"rule_count":record.rules.len(),"apply_status":format!("{:?}",applied.status),"applied_conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"domain_matching_rules":rules});
 std::fs::write(&args[2],serde_json::to_vec_pretty(&out).unwrap()).unwrap();println!("{}",serde_json::json!({"status":out["apply_status"],"matching_rules":out["domain_matching_rules"].as_array().unwrap().len()}));
}
