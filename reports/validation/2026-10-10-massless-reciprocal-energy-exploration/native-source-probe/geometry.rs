use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
fn main(){
 let args:Vec<_>=std::env::args().collect();let input:DensityInput=serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();
 let prepared=input.prepare().unwrap();let family=prepared.occupied_cut(&[0],4096).unwrap().at_physical_masses();
 let mut certified=Vec::new();let energy_offset=family.loops()*(family.loops()+1)/2;
 for slot in family.physical_slots()..family.input_slots(){for shell in family.shells(){
  let energy=&family.coordinates()[energy_offset+shell.loop_index];let AtomView::Var(v)=energy.as_view() else {panic!("energy coordinate")};
  let factor=&family.factors()[slot];let Ok(c)=Rational::try_from(factor.derivative(v.get_symbol()).as_view()) else{continue};
  if c.is_zero()||!(factor-Atom::num(c.clone())*energy).expand().together().cancel().is_zero(){continue}
  assert!(shell.mass_squared.is_zero());assert!(shell.chemical_potential>Rational::zero());
  certified.push(serde_json::json!({"slot":slot,"loop":shell.loop_index,"coefficient":c.to_string(),"energy_coordinate":energy.to_string(),"factor":factor.to_string(),"mass_squared":shell.mass_squared.to_string(),"mu":shell.chemical_potential.to_string(),"exact_factor_equality":true,"production_massive_certificate_absent":family.compact_energy_completion(slot).is_none()}));
 }}
 assert_eq!(certified.len(),1);assert_eq!(certified[0]["slot"],6);assert_eq!(certified[0]["loop"],0);assert_eq!(certified[0]["coefficient"],"1");
 let out=serde_json::json!({"scope":"Regenerate exact occupied geometry from input; algebraic energy descriptor only, not production source/boundary admission","input":args[1],"input_identity":prepared.identity(),"loops":family.loops(),"physical_slots":family.physical_slots(),"input_slots":family.input_slots(),"factors":family.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"coordinates":family.coordinates().iter().map(ToString::to_string).collect::<Vec<_>>(),"roles":family.roles().iter().map(|r|format!("{r:?}")).collect::<Vec<_>>(),"routing":family.inverse_routing().iter().map(|r|r.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),"determinant":family.routing_determinant().to_string(),"certified_completions":certified,"targets":family.targets().iter().map(|r|r.iter().map(|(i,c)|serde_json::json!({"indices":i.0,"coefficient":c.to_string()})).collect::<Vec<_>>()).collect::<Vec<_>>()});
 std::fs::write(&args[2],serde_json::to_vec_pretty(&out).unwrap()).unwrap();println!("exact occupied energy descriptor: slot6=E0; massless; original targets unchanged");
}
