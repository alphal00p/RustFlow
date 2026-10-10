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
 let selected:Vec<_>=record.sources.iter().enumerate().filter(|(_,s)|s.id.starts_with("lorentz/2/0/")).map(|(i,_)|i).collect();assert!(!selected.is_empty());
 let mut specializations=Vec::new();
 for &ordinal in &selected {
  let source=&record.sources[ordinal];
  for power in [1_i64,2] {
   let mut seed=[0_i64;16];seed[..5].copy_from_slice(&[1,power,0,1,0]);seed[5]=-1;
   let variables=polynomial(source.terms[0].coefficient).variables().clone();
   let replacements:BTreeMap<_,_>=indices.iter().enumerate().map(|(axis,&variable)|(variables[variable].to_atom(),Atom::num(seed[axis]))).collect();
   let specialize=|input:Atom| input.replace_map(|view,_,out|{if matches!(view,AtomView::Var(_)){if let Some(value)=replacements.get(&view.to_owned()){**out=value.clone();}}}).together().cancel();
   let mut terms=BTreeMap::<Vec<i64>,Atom>::new();
   for term in &source.terms {
    let coefficient=specialize(polynomial(term.coefficient).to_expression());
    if coefficient.is_zero(){continue;}
    let point=term.integral.iter().enumerate().map(|(axis,&(symbolic,power))|if symbolic {seed[axis]+i64::from(power)}else{i64::from(power)}).collect();
    *terms.entry(point).or_insert_with(Atom::new)+=coefficient;
   }
   terms.retain(|_,c|{*c=c.together().cancel();!c.is_zero()});
   specializations.push(serde_json::json!({"source_ordinal":ordinal,"source_id":source.id,"seed":seed.to_vec(),"source_guard_contains_seed":domain::<16>(&source.domain).contains(&seed),"source_domain":source.domain,
      "nonzero_conditions":source.conditions.iter().map(|&id|specialize(polynomial(id).to_expression()).to_string()).collect::<Vec<_>>(),
      "terms":terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_string()})).collect::<Vec<_>>(),
      "original_terms":source.terms.iter().map(|t|serde_json::json!({"index_offsets":t.integral,"coefficient":polynomial(t.coefficient).to_expression().to_string()})).collect::<Vec<_>>() }));
  }
 }
 let mode=args.get(3).map(String::as_str).unwrap_or("dump");
 let mut probes=Vec::new();
 if mode!="dump" {
  let mut order:Vec<usize>=(0..record.sources.len()).collect();
  match mode {"lorentz-first"=>order.sort_by_key(|i|!selected.contains(i)),"lorentz-only"=>order=selected.clone(),"original"=>{},_=>panic!("invalid mode")}
  let rows=order.iter().map(|&i|{let s=&record.sources[i];GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())}).collect();
  let sources=Arc::new(GuardedSourceSystem::new(format!("{}; existing-source-order-probe={mode}",record.measure),roles,indices,rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
  for power in [2_i64,3] {
   let mut target=[0_i64;16];target[..6].copy_from_slice(&[1,power,0,1,0,-2]);
   let point=IndexDomain::new(std::array::from_fn(|i|IndexBounds::fixed(target[i]))).unwrap();
   for depth in [1_u32,2,3] {
    let started=std::time::Instant::now();
    let found=sources.solve_domains(vec![point.clone()],SearchOptions{max_depth:Some(depth),sample_seed:0,..Default::default()},1).unwrap();
    let unresolved=found.unresolved.iter().map(|u|serde_json::json!({"domain":bounds(&u.domain),"reason":format!("{:?}",u.reason),"detail":u.detail})).collect::<Vec<_>>();
    let program=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();
    let applied=program.apply(&target).unwrap();
    probes.push(serde_json::json!({"mode":mode,"target":target.to_vec(),"max_depth":depth,"max_domains":1,"rules":program.rules().len(),"unresolved":unresolved,"status":format!("{:?}",applied.status),"terms":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"nonzero_conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"seconds":started.elapsed().as_secs_f64()}));
   }
  }
 }
 let result=serde_json::json!({"scope":"exact specialization of original native source rows; optional bounded native solve and original-source replay with only existing source rows/order changed; no installed analytic zero or custom elimination", "source_program":args[1],"source_count":record.sources.len(),"selected_source_ordinals":selected,"specializations":specializations,"probes":probes});
 std::fs::write(&args[2],serde_json::to_vec_pretty(&result).unwrap()).unwrap();
 for s in result["specializations"].as_array().unwrap(){println!("source {} seed{} guard{} ->{}",s["source_id"],s["seed"],s["source_guard_contains_seed"],s["terms"]);}
 for p in result["probes"].as_array().unwrap(){println!("probe mode{} target{} depth{} {} ->{}",p["mode"],p["target"],p["max_depth"],p["status"],p["terms"]);}
}
