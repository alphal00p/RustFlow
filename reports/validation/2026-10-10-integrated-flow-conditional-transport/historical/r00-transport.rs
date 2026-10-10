//! Candidate-only transport through the existing ordinary RustFlow owner.
use serde_json::{json, Value};
use symbolica_amflow::{BoundaryData, ComplexFloat as C, DifferentialSystem, FlowOptions, Precision, RunContext, symbolica};
use symbolica::prelude::*;
fn atom(s:&str)->Atom{Atom::parse(s,"integrated_candidate_transport",Default::default()).unwrap()}
fn render(x:&C)->Value{json!({"re":x.re.to_string(),"im":x.im.to_string()})}
fn main(){
 let args:Vec<_>=std::env::args().collect(); assert_eq!(args.len(),4,"COMPANION COEFFICIENTS OUTPUT");
 let meta:Value=serde_json::from_str(&std::fs::read_to_string(&args[1]).unwrap()).unwrap();
 let data:Value=serde_json::from_str(&std::fs::read_to_string(&args[2]).unwrap()).unwrap()).unwrap();
 let mut records=vec![];
 let eta=symbol!("integrated_candidate_transport::eta");
 let system=DifferentialSystem{variable:eta,matrix:meta["matrix_eta"].as_array().unwrap().iter().map(|row|row.as_array().unwrap().iter().map(|x|atom(x.as_str().unwrap())).collect()).collect()};
 let start_time=std::time::Instant::now();
 for digits in [30,50] { for terms in [32,48,64] { for start in [16,32] {
   let p=Precision::decimal(digits+20).unwrap();
   let compiled=system.compile(p,&Default::default()).unwrap();
   let beta=p.eval(&atom(data["beta"].as_str().unwrap()),&Default::default()).unwrap();
   let t=p.div(&p.i(1),&p.i(start)); let pref=p.pow(&t,&beta);
   let mut values=vec![p.zero();5]; let mut tpower=p.i(1);
   for (n,a) in data["coefficients"].as_array().unwrap().iter().take(terms).enumerate(){
      let c=p.eval(&atom(a.as_str().unwrap()),&Default::default()).unwrap();
      let term=p.mul(&p.mul(&pref,&tpower),&c);
      let b=p.add(&p.i(n as i64),&beta);
      for k in 0..4 { values[k]=p.add(&values[k],&p.mul(&term,&p.powi(&b,k as i64))); }
      tpower=p.mul(&tpower,&t);
   }
   values[4]=pref;
   let boundary=BoundaryData{point:p.i(start),values};
   let options=FlowOptions{digits,guard_digits:20,series_order:64,max_steps:128,max_precision_attempts:1,max_boundary_attempts:1,..Default::default()};
   let mut current=boundary; let mut profile_steps=0usize;
   for end in [16,8,4,2,1] {
      let begin=std::time::Instant::now();
      let remaining=128usize.saturating_sub(profile_steps);
      if remaining==0 {records.push(json!({"digits":digits,"terms":terms,"start_eta":start,"end_eta":end,"error":"profile step cap128"}));break;}
      let result=compiled.transport(&current,&[p.i(end)],&FlowOptions{max_steps:remaining,..options.clone()},&RunContext::default());
      match result {
       Ok(out)=>{ profile_steps+=out.diagnostics.steps; records.push(json!({"digits":digits,"working_digits":digits+20,"boundary_terms":terms,"start_eta":start,"end_eta":end,"state":out.values.iter().map(render).collect::<Vec<_>>(),"steps":out.diagnostics.steps,"profile_steps":profile_steps,"seconds":begin.elapsed().as_secs_f64()})); current=BoundaryData{point:p.i(end),values:out.values}; }
       Err(e)=>{records.push(json!({"digits":digits,"boundary_terms":terms,"start_eta":start,"end_eta":end,"error":format!("{e:?}"),"seconds":begin.elapsed().as_secs_f64()}));break;}
      }
      std::fs::write(&args[3],serde_json::to_string_pretty(&json!({"status":"running_candidate_only","identity_certified":false,"records":records})).unwrap()).unwrap();
   }
 }}}
 std::fs::write(&args[3],serde_json::to_string_pretty(&json!({"status":"completed_candidate_only","identity_certified":false,"source_equation":"finite-prefix recurrence, not certified identity","target_index":0,"normalization":"native hard master and compact angular factors omitted","endpoint_evaluated":false,"reference_values_read":false,"seconds":start_time.elapsed().as_secs_f64(),"records":records})).unwrap()).unwrap();
}
