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
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let p=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();let point=[1,0,0,2,2,-1,-1,-1,0,0,0,0];let applied=p.apply(&point).unwrap();

 let requested_file=std::fs::read(&args[3]).unwrap();let metadata:serde_json::Value=serde_json::from_slice(&requested_file).unwrap();let requests:Vec<[i64;12]>=serde_json::from_value(metadata["requested"].clone()).unwrap();
 let mut failed=Vec::new();let mut no_rule=std::collections::BTreeSet::new();let mut all_failed=std::collections::BTreeSet::new();let mut applications=0;
 for target in &requests {let r=p.reduce(*target,Default::default()).unwrap();applications+=r.rule_applications;for u in r.unresolved {if u.reason==rustred::solver::guarded::GuardedApplicationFailure::NoApplicableRule {no_rule.insert(u.integral);}else {all_failed.insert(u.integral);failed.push(serde_json::json!({"target":target,"integral":u.integral,"reason":format!("{:?}",u.reason),"coefficient":u.coefficient.to_expression().to_string()}));}}}
 let out=serde_json::json!({"scope":"Full native replay reduction of actual discovery-request set; no search","program":args[1],"requested_metadata":args[3],"requests":requests.len(),"point_directly_requested":requests.contains(&point),"point_in_failed_descendants":all_failed.contains(&point),"rule_applications":applications,"no_rule_labels":no_rule,"failed_labels":all_failed,"failures":failed});
 std::fs::write(&args[2],serde_json::to_vec_pretty(&out).unwrap()).unwrap();println!("requests={} applications={} failures={} no_rule={}",requests.len(),applications,failed.len(),no_rule.len());
}
