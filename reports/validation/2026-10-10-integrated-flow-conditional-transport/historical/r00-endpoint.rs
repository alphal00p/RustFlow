//! Conditional candidate endpoint through existing Frobenius/matching/limit APIs.
use serde_json::{json,Value};
use symbolica_amflow::{BoundaryData,DifferentialSystem,FlowOptions,Precision,RunContext,symbolica};
use symbolica::prelude::*;
fn atom(s:&str)->Atom{Atom::parse(s,"integrated_candidate_endpoint",Default::default()).unwrap()}
fn main(){
 let a:Vec<_>=std::env::args().collect(); assert_eq!(a.len(),4);
 let meta:Value=serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
 let finite:Value=serde_json::from_str(&std::fs::read_to_string(&a[2]).unwrap()).unwrap();
 let system=DifferentialSystem{variable:symbol!("integrated_candidate_endpoint::eta"),matrix:meta["matrix_eta"].as_array().unwrap().iter().map(|r|r.as_array().unwrap().iter().map(|x|atom(x.as_str().unwrap())).collect()).collect()};
 let run=RunContext::default();let started=std::time::Instant::now();
 let prepared=match system.prepare_frobenius(&run){Ok(v)=>v,Err(e)=>{std::fs::write(&a[3],json!({"status":"candidate_endpoint_preparation_failed","error":format!("{e:?}"),"identity_certified":false}).to_string()).unwrap();return;}};
 let mut records=vec![];
 for r in finite["records"].as_array().unwrap().iter().filter(|r|r["end_eta"]==1 && r.get("state").is_some()) {
  let digits=r["digits"].as_u64().unwrap() as u32;let terms=r["boundary_terms"].as_u64().unwrap() as usize;let start=r["start_eta"].as_i64().unwrap();
  let p=Precision::decimal(digits+20).unwrap();let compiled=system.compile(p,&Default::default()).unwrap();
  let basis=match prepared.evaluate(p,&Default::default(),terms,&run){Ok(b)=>b,Err(e)=>{records.push(json!({"digits":digits,"terms":terms,"start_eta":start,"error":format!("{e:?}")}));continue;}};
  let boundary=BoundaryData{point:p.i(1),values:r["state"].as_array().unwrap().iter().map(|x|p.parse(x["re"].as_str().unwrap(),x["im"].as_str().unwrap()).unwrap()).collect()};
  let matches=if digits==50 && terms==64 && start==32 {vec![8,16]} else {vec![8]};
  for denominator in matches {
   let begin=std::time::Instant::now();let point=p.scale(&p.i(1),1,denominator);
   let result=(||->symbolica_amflow::Result<Value>{
    let out=compiled.transport(&boundary,&[point.clone()],&FlowOptions{digits,guard_digits:20,series_order:64,max_steps:128-r["profile_steps"].as_u64().unwrap() as usize,max_precision_attempts:1,max_boundary_attempts:1,..Default::default()},&run)?;
    let constants=basis.match_values(&point,&out.values,&Default::default())?;
    let endpoint=basis.physical_limit(&constants,symbol!("unused_epsilon"))?;
    Ok(json!({"endpoint":endpoint.iter().map(|x|json!({"re":x.re.to_string(),"im":x.im.to_string()})).collect::<Vec<_>>(),"matching_constants":constants.iter().map(|x|json!({"re":x.re.to_string(),"im":x.im.to_string()})).collect::<Vec<_>>(),"steps":out.diagnostics.steps}))
   })();
   records.push(json!({"digits":digits,"working_digits":digits+20,"boundary_terms":terms,"frobenius_terms":terms,"start_eta":start,"match_eta":format!("1/{denominator}"),"seconds":begin.elapsed().as_secs_f64(),"result":match result{Ok(v)=>v,Err(e)=>json!({"error":format!("{e:?}")})}));
   std::fs::write(&a[3],serde_json::to_string_pretty(&json!({"status":"running_candidate_only","identity_certified":false,"records":records})).unwrap()).unwrap();
  }
 }
 std::fs::write(&a[3],serde_json::to_string_pretty(&json!({"status":"completed_candidate_endpoint_diagnostic_only","identity_certified":false,"physical_endpoint_accepted":false,"ordinary_owner":"PreparedFrobenius::evaluate, match_values, physical_limit; numerical tolerance policy inherited","exponents":prepared.exponents().map(|x|x.to_string()).collect::<Vec<_>>(),"nonzero_conditions":prepared.nonzero_conditions().iter().map(|x|x.to_string()).collect::<Vec<_>>(),"records":records,"seconds":started.elapsed().as_secs_f64()})).unwrap()).unwrap();
}
