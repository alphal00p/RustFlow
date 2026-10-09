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
fn run<const N:usize>(bytes:Vec<u8>, target:[i64;N], mode:&str, output:&str){
 let limits=BinaryIoLimits::default();let env=inspect_program(&bytes,limits).unwrap();
 let (r,used):(Rec,usize)=bincode::decode_from_slice(env.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,env.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(r.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(env.section(SectionTag::SYMBOLICA_STATE).unwrap(),env.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let p=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(p.denominator.is_one());p.numerator};
 let roles:[IndexRole;N]=r.roles.iter().map(|r|match r {0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=r.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(r.measure.clone(),roles,r.indices.try_into().unwrap(),rows).unwrap().with_zero_domains(r.zero_domains.iter().map(domain).collect()).unwrap());
 let replay=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();
 let application=replay.apply(&target).unwrap();
 eprintln!("Historical replay passed: {} sources, {} rules; target {:?}",sources.sources().len(),replay.rules().len(),application.status);
 let point=IndexDomain::new(target.map(IndexBounds::fixed)).unwrap();
 let d=match mode {
 "point"=>point,
 "one-ray"=>{let mut b=*point.bounds();b[6]=IndexBounds::new(Some(target[6]),None).unwrap();IndexDomain::new(b).unwrap()},
 "zero-face"=>IndexDomain::new(std::array::from_fn(|i|if i>=if N==16{14}else{14} {IndexBounds::fixed(target[i])}else if target[i]==0 {IndexBounds::fixed(0)}else if target[i]>0 {IndexBounds::new(Some(1),None).unwrap()}else {IndexBounds::new(None,Some(-1)).unwrap()})).unwrap(),
 "sector"=>IndexDomain::new(std::array::from_fn(|i|if roles[i]==IndexRole::Occupation||i>=if N==16{16}else{18}{IndexBounds::fixed(target[i])}else if target[i]>0 {IndexBounds::new(Some(1),None).unwrap()}else{IndexBounds::new(None,Some(0)).unwrap()})).unwrap(),_=>panic!("unknownmode")};
 let depth=std::env::var("PROBE_DEPTH").ok().map(|s|s.parse().unwrap()).unwrap_or(3);let budget=std::env::var("PROBE_DOMAINS").ok().map(|s|s.parse().unwrap()).unwrap_or(256);
 let started=std::time::Instant::now();let found=sources.solve_domains(vec![d.clone()],SearchOptions{max_depth:Some(depth),sample_seed:0,..Default::default()},budget).unwrap();let elapsed=started.elapsed().as_secs_f64();
 let gaps=found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
 let program=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();let a=program.apply(&target).unwrap();let reduction=program.reduce(target,Default::default()).unwrap();
 let rules=program.rules().iter().enumerate().map(|(i,r)|serde_json::json!({"ordinal":i,"domain":bounds(r.domain()),"target":r.candidate().target.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"sources":r.candidate().sources.iter().map(|s|serde_json::json!({"id":sources.sources()[s.basis_row].id,"integral":s.seed.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"shifts":s.seed.shifts.to_vec()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
 let report=serde_json::json!({"scope":"historical generated source context reconstructed losslessly, then native replay and bounded native discovery; no production changes or numerical predictions","capacity":N,"historical_rule_count":replay.rules().len(),"historical_application":format!("{:?}",application.status),"source_count":sources.sources().len(),"mode":mode,"requested_domain":bounds(&d),"target":target.to_vec(),"depth":depth,"max_domains":budget,"discovery_seconds":elapsed,"rules":rules,"unresolved":gaps,"application":format!("{:?}",a.status),"application_terms":a.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"reduction_rule_applications":reduction.rule_applications,"reduction_unresolved":reduction.unresolved.iter().map(|t|serde_json::json!({"indices":t.integral.to_vec(),"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>(),"conditions":reduction.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()});
 std::fs::write(output,serde_json::to_vec_pretty(&report).unwrap()).unwrap();std::fs::write(format!("{output}.bin"),program.encode_native(limits).unwrap()).unwrap();eprintln!("Result {:?}: {} rules, {} gaps, {} unresolved reductions in {:.3}s",a.status,program.rules().len(),report["unresolved"].as_array().unwrap().len(),reduction.unresolved.len(),elapsed);
}
fn main(){let a=std::env::args().collect::<Vec<_>>();let bytes=std::fs::read(&a[1]).unwrap();match a[2].as_str(){"single"=>run::<16>(bytes,[1,1,0,1,1,1,2,-2,0,0,0,0,0,0,0,0],&a[3],&a[4]),"double"=>run::<20>(bytes,[1,1,1,2,1,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0],&a[3],&a[4]),_=>panic!()}}
