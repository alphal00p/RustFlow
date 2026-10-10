#![allow(dead_code,unused_imports)]
use symbolica_amflow::*;
use symbolica::prelude::*;
use std::collections::{BTreeMap,BTreeSet};
mod algebra;mod cut_regions;mod coefficient;mod moments;
use finite_density::boundary::OccupiedBoundaryDistribution;
use reduction::ReductionBackend;
fn str_a(a:&Atom)->String{a.to_canonical_string()}
fn family_key(f:&IntegralFamily)->String{serde_json::json!({"loops":f.loops.len(),"external":f.external,"gram":f.external_gram.iter().map(|r|r.iter().map(str_a).collect::<Vec<_>>()).collect::<Vec<_>>(),"physical":f.physical_propagators,"dimension":f.dimension,"epsilon":Atom::var(f.epsilon).to_canonical_string(),"propagators":f.propagators.iter().map(|p|serde_json::json!({"constant":str_a(&p.constant),"coefficients":p.scalar_products.iter().map(str_a).collect::<Vec<_>>()})).collect::<Vec<_>>()}).to_string()}
fn polynomial_zero(loops:usize,eps:Symbol)->Result<serde_json::Value>{
 use rustred::sector::{Mask,zero::{Analyzer,Decision}};
 let mut propagators=Vec::new();for i in 0..loops{let mut row=vec![0;loops];row[i]=1;propagators.push(Propagator::quadratic(&row,&[],Atom::zero(),&[])?)}for i in 0..loops{for j in i+1..loops{let mut row=vec![0;loops];row[i]=1;row[j]=-1;propagators.push(Propagator::quadratic(&row,&[],Atom::zero(),&[])?);}}
 let f=IntegralFamily{name:"unrestricted_polynomial_virtual_factor".into(),loops:(0..loops).map(|i|format!("k{i}")).collect(),external:vec![],external_gram:vec![],physical_propagators:propagators.len(),propagators,epsilon:eps,dimension:4};
 let ctx=rustred::algebra::CoefficientContext::try_new(["polynomial_vacuum_D"]).map_err(|e|Error::InvalidInput(e.to_string()))?;
 let c=|a:&Atom|a.try_to_rational_polynomial(&Q,&Z,Some(ctx.one().get_variables().clone())).map_err(|e|Error::InvalidInput(e.to_string()));
 let den=f.propagators.iter().map(|p|Ok(rustred::family::AffineDenominator::new(c(&p.constant)?,p.scalar_products.iter().map(c).collect::<Result<Vec<_>>>()?))).collect::<Result<Vec<_>>>()?;
 let native=rustred::family::IntegralFamily::new(f.name.clone(),f.loops.clone(),vec![],ctx.clone(),ctx.parameter("polynomial_vacuum_D").unwrap(),den,vec![],vec![ctx.zero();f.propagators.len()]).map_err(|e|Error::InvalidInput(e.to_string()))?;
 let a=Analyzer::try_unrestricted(&native).map_err(|e|Error::Reduction(e.to_string()))?;let m=Mask::try_new(vec![false;f.propagators.len()]).map_err(|e|Error::Reduction(e.to_string()))?;
 let Decision::ProvedZero(c)=a.analyze(&m).map_err(|e|Error::Reduction(e.to_string()))?else{return Err(Error::Unsupported("native polynomial vacuum not proved zero".into()))};
 Ok(serde_json::json!({"family":family_key(&f),"native_certificate":format!("{c:?}"),"scope":"all finite polynomial virtual tensor moments; compact coefficients untouched"}))
}
fn main()->Result<()>{
 let args=std::env::args().collect::<Vec<_>>();let out=std::path::PathBuf::from(&args[2]);std::fs::create_dir_all(&out).unwrap();let half_order:usize=args[3].parse().unwrap();if half_order>16{return Err(Error::Limit("initial half-order cap16".into()))}
 let input:finite_density::DensityInput=serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();let prepared=input.prepare()?;let family=prepared.occupied_cut(&[0,3],4096)?.at_physical_masses();let eps=symbol!("integrated_moment::epsilon");
 let shifted_slots=(0..family.physical_slots()).filter(|&s|family.roles()[s]==rustred::solver::guarded::IndexRole::Ordinary).collect::<Vec<_>>();let evidence=finite_density::massless_endpoint::MasslessFlowEvidence::new(&prepared,&family,&shifted_slots,Default::default())?;let labels=family.targets().iter().flat_map(|t|t.keys().cloned()).collect::<Vec<_>>();let audit=evidence.validate_labels(&family,&labels)?;
 let ordinary=family.region_family(eps,4)?;let compact=(0..family.loops()).map(|i|family.shells().iter().any(|s|s.loop_index==i)).collect::<Vec<_>>();let shifted=(0..family.input_slots()).map(|s|shifted_slots.contains(&s)).collect::<Vec<_>>();let context=RunContext::default();let regions=cut_regions::enumerate_compact(&ordinary,&compact,&context)?;
 let backend=RustRedBackend{bubble_subloops:false,max_targets:4096,..Default::default()};let mut hard_cache:BTreeMap<String,Vec<(String,Atom)>>=BTreeMap::new();let mut hard_reports=Vec::new();let mut conditions=BTreeSet::new();let mut totals:BTreeMap<(usize,String,String),Atom>=BTreeMap::new();let mut records=Vec::new();let mut zeros=Vec::new();let mut processed=0usize;
 for (ri,region)in regions.iter().enumerate(){
  let transformed=ordinary.transform_loops(&region.transformation)?;let virtual_soft=(0..family.loops()).filter(|&i|!compact[i]&&!region.hard[i]).count();
  for(ti,target)in family.targets().iter().enumerate(){for(original,weight)in target{
   let indices=Integral(original.0[..family.input_slots()].iter().enumerate().map(|(s,&n)|if family.roles()[s]==rustred::solver::guarded::IndexRole::RequiredCut{0}else{n}).collect());
   let ds=family.shells().iter().map(|s|OccupiedBoundaryDistribution{source_loop_index:s.loop_index,shell:finite_density::compact::CompactShell{mass_squared:Rational::zero(),chemical_potential:s.chemical_potential.clone()},cut_index:original.0[s.physical_slot],upper_index:original.0[s.upper_slot],lower_index:original.0[s.lower_slot]}).collect::<Vec<_>>();
   let ex=regions::expand_region(&ordinary,&indices,&shifted,region,half_order)?;
   for(grade,expr)in ex.coefficients.iter().enumerate(){let projected=integrand::projected_factor_region(expr,&ex.coordinates,&transformed,&region.hard,200000)?;let mut integrated=Vec::new();
    for term in &projected.terms{
     if virtual_soft>0{let count=moments::validate_polynomial(&term.soft,&projected.soft)?;zeros.push(serde_json::json!({"region":ri,"target":ti,"indices":original.0,"half_grade":grade,"monomials":count,"proof":polynomial_zero(virtual_soft,eps)?}));continue}
     let ds_ordered=projected.soft.source_loop_indices.iter().map(|i|ds.iter().find(|d|d.source_loop_index==*i).unwrap().clone()).collect::<Vec<_>>();
     let soft=moments::integrate(&term.soft,&projected.soft,&ds_ordered)?;if soft.is_zero(){continue}
     let cache_key=format!("{}:{}",family_key(&projected.hard.template),str_a(&term.hard));
     if !hard_cache.contains_key(&cache_key){let mut result=Vec::new();for ordinary_term in integrand::to_integrals(&term.hard,&projected.hard.coordinates,&projected.hard.template,200000)?{
      let r=backend.reduce(&ordinary_term.family,std::slice::from_ref(&ordinary_term.integral),&context)?;
      for c in &r.nonzero_conditions{conditions.insert(str_a(c));}
      let mapped=r.rules.get(&ordinary_term.integral).cloned().unwrap_or_else(||BTreeMap::from([(ordinary_term.integral.clone(),Atom::one())]));
      hard_reports.push(serde_json::json!({"family":family_key(&ordinary_term.family),"input":ordinary_term.integral.0,"input_coefficient":str_a(&ordinary_term.coefficient),"masters":r.residuals.iter().map(|i|&i.0).collect::<Vec<_>>(),"rule":mapped.iter().map(|(i,c)|serde_json::json!({"indices":i.0,"coefficient":str_a(c)})).collect::<Vec<_>>(),"conditions":r.nonzero_conditions.iter().map(str_a).collect::<Vec<_>>()}));
      for(i,c)in mapped{let channel=serde_json::json!({"family":family_key(&ordinary_term.family),"master":i.0}).to_string();result.push((channel,(&ordinary_term.coefficient*c).together().cancel()));}
     }hard_cache.insert(cache_key.clone(),result);}
     let exponent=(&ex.eta_power-Atom::num((grade as i64,2))).together().cancel();
     for(channel,hard)in &hard_cache[&cache_key]{let value=(weight*&soft*hard).together().cancel();*totals.entry((ti,str_a(&exponent),channel.clone())).or_default()+=&value;integrated.push(serde_json::json!({"master":channel,"coefficient":str_a(&value)}));}
     processed+=1;
    }
    records.push(serde_json::json!({"region":ri,"target":ti,"indices":original.0,"eta_power":str_a(&ex.eta_power),"half_grade":grade,"integrated":integrated}));
   }
  }}
 }
 let points=BTreeMap::from([(Atom::var(eps),Atom::num((-5,4)))]);
 let coefficients=totals.into_iter().map(|((target,power,master),c)|{let c=c.together().cancel();serde_json::json!({"target":target,"eta_exponent":power,"master":master,"coefficient":str_a(&c),"at_D_13_over_2":str_a(&family::substitute(&c,&points).together().cancel())})}).collect::<Vec<_>>();
 let report=serde_json::json!({"schema":"rustflow.integrated-moment-series.v1","scope":"Source-derived formal large-eta coefficients; exact ordinary RustRed plus HEPKit moments; no reference values; no weighted reduction","input_identity":prepared.identity(),"cut_slots":[0,3],"shifted_slots":shifted_slots,"half_order":half_order,"field":"Q(epsilon)","D":"4-2 epsilon","compact_normalization_per_shell":"mu^(D-2)/(sqrt(pi)*Gamma((D-1)/2))","mu":family.shells().iter().map(|s|s.chemical_potential.to_string()).collect::<Vec<_>>(),"routing_determinant":family.routing_determinant().to_string(),"all_region_count":regions.len(),"native_polynomial_zero_records":zeros,"ordinary_backend":backend.identity(),"hard_reductions":hard_reports,"original_label_audit":format!("{audit:?}"),"origin_identity":evidence.source_identity(),"nonzero_conditions":conditions,"term_records":records,"coefficients":coefficients,"processed_nonzero_products":processed,"hard_cache_entries":hard_cache.len(),"unknown_tail":true,"convergent_germ_proved":false});
 let bytes=serde_json::to_vec_pretty(&report).unwrap();if bytes.len()>128*1024*1024{return Err(Error::Limit("retained algebra report128MiB".into()))}std::fs::write(out.join("coefficients.json"),bytes).unwrap();println!("regions={} products={} hard_cache={} coefficients={}",regions.len(),processed,hard_cache.len(),report["coefficients"].as_array().unwrap().len());Ok(())
}
