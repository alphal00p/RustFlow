#![allow(dead_code,unused_imports)]
//! Isolated exact census; does not construct a weighted source or reduction.
use symbolica_amflow::*;
use symbolica::prelude::*;
use std::collections::BTreeMap;
mod algebra;
mod cut_regions;
fn strings(xs:&[Atom])->Vec<String>{xs.iter().map(AtomCore::to_canonical_string).collect()}
fn targets(f:&finite_density::geometry::OccupiedCutFamily)->serde_json::Value {
 serde_json::json!(f.targets().iter().map(|t|t.iter().map(|(i,c)|serde_json::json!({"indices":i.0,"coefficient":c.to_canonical_string()})).collect::<Vec<_>>()).collect::<Vec<_>>())
}
fn main()->Result<()> {
 let args=std::env::args().collect::<Vec<_>>();
 let input:finite_density::DensityInput=serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();
 let prepared=input.prepare()?; let independent=prepared.occupied_cut(&[0,3],4096)?;let family=independent.at_physical_masses();
 let eta=symbol!("integrated_moment::eta");let eps=symbol!("integrated_moment::epsilon");
 let shifted_slots=(0..family.physical_slots()).filter(|&s|family.roles()[s]==rustred::solver::guarded::IndexRole::Ordinary).collect::<Vec<_>>();
 let measure=family.deformed_measure::<16>(eta,&shifted_slots)?;
 let ordinary=family.region_family(eps,4)?;let compact=(0..family.loops()).map(|i|family.shells().iter().any(|s|s.loop_index==i)).collect::<Vec<_>>();
 let shifted=(0..family.input_slots()).map(|s|shifted_slots.contains(&s)).collect::<Vec<_>>();
 let regions=cut_regions::enumerate_compact(&ordinary,&compact,&RunContext::default())?;
 let exact_targets=family.targets().iter().map(|t|t.iter().fold(Atom::zero(),|sum,(i,c)|sum+i.0.iter().enumerate().filter(|(s,_)|*s<family.input_slots()&&family.roles()[*s]==rustred::solver::guarded::IndexRole::Ordinary).fold(c.clone(),|a,(s,&p)|a*measure.factors()[s].pow(-i64::from(p))))).collect::<Vec<_>>();
 let mut census=Vec::new();
 for (r,region) in regions.iter().enumerate(){
  let transformed=ordinary.transform_loops(&region.transformation)?;let mut labels=Vec::new();
  for (ti,target) in family.targets().iter().enumerate(){for (original,coefficient) in target{
   let indices=Integral(original.0[..family.input_slots()].iter().enumerate().map(|(s,&i)|if family.roles()[s]==rustred::solver::guarded::IndexRole::RequiredCut{0}else{i}).collect());
   let ex=regions::expand_region(&ordinary,&indices,&shifted,region,8)?;
   let mut pieces=Vec::new();
   for (grade,expr) in ex.coefficients.iter().enumerate(){
    let p=integrand::projected_factor_region(expr,&ex.coordinates,&transformed,&region.hard,200000)?;
    pieces.push(serde_json::json!({"half_grade":grade,"before_projection":expr.to_canonical_string(),"hard_coordinates":strings(&p.hard.coordinates),"soft_coordinates":strings(&p.soft.coordinates),"hard_source_loops":p.hard.source_loop_indices,"soft_source_loops":p.soft.source_loop_indices,"terms":p.terms.iter().map(|t|serde_json::json!({"hard":t.hard.to_canonical_string(),"soft":t.soft.to_canonical_string()})).collect::<Vec<_>>() }));
   }
   labels.push(serde_json::json!({"target":ti,"indices":original.0,"coefficient":coefficient.to_canonical_string(),"region_integral":indices.0,"eta_power":ex.eta_power.to_canonical_string(),"coefficients":pieces}));
  }}
  census.push(serde_json::json!({"id":r,"hard":region.hard,"hard_branches":region.hard_branches,"transformation":region.transformation.iter().map(|row|strings(row)).collect::<Vec<_>>(),"determinant":region.jacobian_determinant.to_canonical_string(),"labels":labels}));
 }
 let output=serde_json::json!({"scope":"Exact source-bound all-region census and half-grades0..8; no integrated coefficient/reference/closure claim","input_identity":prepared.identity(),"cut_slots":[0,3],"shifted_slots":shifted_slots,"coordinates":strings(family.coordinates()),"independent_factors":strings(independent.factors()),"physical_factors":strings(family.factors()),"eta_factors":strings(measure.factors()),"roles":family.roles().iter().map(|r|format!("{r:?}")).collect::<Vec<_>>(),"independent_targets":targets(&independent),"physical_targets":targets(&family),"off_shell_weighted_expressions_without_C_H":strings(&exact_targets),"routing":family.inverse_routing().iter().map(|row|strings(row)).collect::<Vec<_>>(),"routing_determinant":family.routing_determinant().to_string(),"shells":family.shells().iter().map(|s|serde_json::json!({"loop_index":s.loop_index,"physical_slot":s.physical_slot,"upper_slot":s.upper_slot,"lower_slot":s.lower_slot,"mass_squared":s.mass_squared.to_canonical_string(),"chemical_potential":s.chemical_potential.to_string()})).collect::<Vec<_>>(),"regions":census});
 std::fs::write(&args[2],serde_json::to_vec_pretty(&output).unwrap()).unwrap();
 println!("regions={} targets={} half_grades=9",regions.len(),family.targets().len()); Ok(())
}
