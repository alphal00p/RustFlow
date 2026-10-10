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
 let point=[1,1,1,1,1,0,0,0,0,0,0,0];let mut reports=Vec::new();
 let modes=vec![("full62",(0..record.sources.len()).collect::<Vec<_>>()),("row4",vec![4]),("row60",vec![60]),("rows4_60",vec![4,60])];
 for (name,ordinals) in modes {
 let rows=ordinals.iter().map(|&n| {let s=&record.sources[n];GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())}).collect();
 let sources=Arc::new(GuardedSourceSystem::new(format!("{};ray-point-fixture-subset={name}",record.measure),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 for phase in ["ray","point"]{
 let request=if phase=="point" {IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()}else{IndexDomain::new(point.map(|p|if p==0{IndexBounds::fixed(0)}else{IndexBounds::new(Some(p),None).unwrap()})).unwrap()};
 let found=sources.solve_domains(vec![request],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},1).unwrap();let gaps=found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();let program=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();let result=program.apply(&point).unwrap();let bytes=program.encode_native(limits).unwrap();let decoded=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();let replay=decoded.apply(&point).unwrap();assert_eq!(format!("{:?}",result.status),format!("{:?}",replay.status));assert_eq!(result.terms,replay.terms);assert_eq!(result.nonzero_conditions,replay.nonzero_conditions);
 std::fs::write(std::path::Path::new(&args[2]).join(format!("{name}-{phase}.bin")),bytes).unwrap();reports.push(serde_json::json!({"mode":name,"original_source_ordinals":ordinals,"phase":phase,"target":point,"source_count":sources.sources().len(),"rules":program.rules().len(),"status":format!("{:?}",result.status),"gaps":gaps,"rhs":result.terms.iter().map(|(i,c)|serde_json::json!({"integral":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"roundtrip_exact":true}));
 }}
 std::fs::write(std::path::Path::new(&args[2]).join("result.json"),serde_json::to_vec_pretty(&reports).unwrap()).unwrap();for row in &reports {println!("{} {} {} rhs={}",row["mode"],row["phase"],row["status"],row["rhs"].as_array().unwrap().len());}
}
