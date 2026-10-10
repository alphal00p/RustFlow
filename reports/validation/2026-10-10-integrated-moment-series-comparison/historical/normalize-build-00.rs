//! Report-only saved-channel normalization. Never constructs moment coefficients.
use symbolica_amflow::{*, finite_density::normalization::native_measure_to_euclidean};
use symbolica::prelude::*;
use serde_json::{Value,json};
use std::{collections::BTreeSet,fs};
fn parse(s:&str)->Result<Atom>{Atom::parse(s,"integrated_moment_adapter",Default::default()).map_err(|e|Error::InvalidInput(e.to_string()))}
fn rational(s:&str)->Result<Rational>{Rational::try_from(parse(s)?.as_view()).map_err(|e|Error::InvalidInput(e.to_string()))}
fn out(c:&ComplexFloat)->Value{json!({"re":c.re.to_string(),"im":c.im.to_string()})}
fn main()->Result<()>{
 let args=std::env::args().collect::<Vec<_>>();
 if args.len()!=4{return Err(Error::InvalidInput("usage: normalize COEFFICIENTS INPUT OUTPUT".into()))}
 let bytes=fs::read(&args[1]).map_err(|e|Error::InvalidInput(e.to_string()))?;
 let data:Value=serde_json::from_slice(&bytes).map_err(|e|Error::InvalidInput(e.to_string()))?;
 if data["schema"]!="integrated-gaussian-rational-prefix.v1" || data["dimension"]!="13/2" || data["q"]!=1 || data["tail"]!="unknown"{return Err(Error::Unsupported("unsupported saved prefix schema; no missing data invented".into()))}
 let input:finite_density::DensityInput=serde_json::from_slice(&fs::read(&args[2]).map_err(|e|Error::InvalidInput(e.to_string()))?).map_err(|e|Error::InvalidInput(e.to_string()))?;
 let prepared=input.prepare()?;let family=prepared.occupied_cut(&[0,3],4096)?.at_physical_masses();
 if data["input_identity"]!=prepared.identity(){return Err(Error::InvalidInput("saved coefficient input identity mismatch".into()))}
 if data["compact_normalization"]!="prod_i mu_i^(D-2)/(sqrt(pi)*Gamma((D-1)/2))" { return Err(Error::Unsupported("unknown producer normalization contract".into())); }
 let D=Rational::from((13,2));let mut channels=BTreeSet::new();
 for row in data["coefficients"].as_array().ok_or_else(||Error::InvalidInput("missing coefficient records".into()))?{channels.insert(row["master"].as_str().ok_or_else(||Error::InvalidInput("missing source-bound channel".into()))?.to_owned());}
 if channels.is_empty(){return Err(Error::InvalidInput("empty channels".into()))}
 let mut profiles=Vec::new();
 for digits in [30,50]{
  let p=Precision::decimal(digits)?;
  let parameters=ahash::HashMap::from_iter([(parse("integrated_moment::{}::epsilon")?,p.rational(&Rational::from((-5,4))))]);
  let mut conditions=Vec::new();
  for condition in data["conditions"].as_array().ok_or_else(||Error::InvalidInput("missing retained conditions".into()))? { let raw=condition.as_str().ok_or_else(||Error::InvalidInput("condition is not an exact expression".into()))?;let value=p.eval(&parse(raw)?,&parameters)?;if value==p.zero(){return Err(Error::Numerical("retained producer condition vanished".into()))}conditions.push(json!({"condition":raw,"value":out(&value)})); }
  let det=family.routing_determinant();let abs_det=if det<&Rational::zero(){-det.clone()}else{det.clone()};
  let measure_atom=native_measure_to_euclidean(family.loops(),family.shells().len(),&Atom::num(D.clone()))?*Atom::num(abs_det).pow(Atom::num(-D.clone()));
  let measure=p.eval(&measure_atom,&ahash::HashMap::default())?;
  let sqrt_pi=p.eval(&Atom::var(Symbol::PI).pow(Atom::num((1,2))),&ahash::HashMap::default())?;
  let gamma=p.gamma_real(&p.rational(&((&D-Rational::one())/&Rational::from(2))).re)?;
  let mut compact=p.i(1);
  for shell in family.shells(){
   let mu=p.rational(&shell.chemical_potential);
   compact=p.mul(&compact,&p.div(&p.pow(&mu,&p.rational(&(&D-Rational::from(2)))),&p.mul(&sqrt_pi,&gamma)));
  }
  let mut masters=Vec::new();
  for channel in &channels{
   let key:Value=serde_json::from_str(channel).map_err(|e|Error::InvalidInput(e.to_string()))?;
   let f:Value=serde_json::from_str(key["family"].as_str().ok_or_else(||Error::InvalidInput("family identity absent".into()))?).map_err(|e|Error::InvalidInput(e.to_string()))?;
   let indices=key["master"].as_array().ok_or_else(||Error::InvalidInput("master label absent".into()))?;
   if f["loops"]!=1 || f["external"]!=json!([]) || f["physical"]!=1 || f["D"]!="13/2" || indices.len()!=1 || f["propagators"].as_array().map(Vec::len)!=Some(1){return Err(Error::Unsupported("this adapter only admits ordinary one-factor Gaussian masters; other channels require a declared native ordinary owner".into()))}
   let power=indices[0].as_u64().ok_or_else(||Error::Unsupported("nonpositive ordinary master index".into()))?;
   if power==0 || power>u32::MAX as u64{return Err(Error::InvalidInput("master power range".into()))}
   let prop=&f["propagators"][0];let coeff=prop["coefficients"].as_array().ok_or_else(||Error::InvalidInput("quadratic coefficient array absent".into()))?;
   if coeff.len()!=1{return Err(Error::Unsupported("one-factor master must have one virtual scalar product".into()))}
   let scale=rational(coeff[0].as_str().ok_or_else(||Error::InvalidInput("scale absent".into()))?)?;
   let constant=rational(prop["constant"].as_str().ok_or_else(||Error::InvalidInput("mass constant absent".into()))?)?;
   if scale<=Rational::zero(){return Err(Error::Unsupported("positive exact quadratic scale required".into()))}
   let mass=-constant/&scale;if mass<=Rational::zero(){return Err(Error::Unsupported("positive exact Gaussian mass required".into()))}
   let value=p.mul(&p.powi(&p.rational(&scale),-(power as i64)),&vacuum::tadpole(power as u32,&p.rational(&mass),&D,p)?);
   masters.push(json!({"channel":channel,"owner":"symbolica_amflow::vacuum::tadpole","power":power,"mass_squared":mass.to_string(),"quadratic_scale":scale.to_string(),"value":out(&value)}));
  }
  profiles.push(json!({"digits":digits,"bits":p.bits,"native_measure_to_euclidean_atom":measure_atom.to_canonical_string(),"native_measure_to_euclidean":out(&measure),"compact_normalization":out(&compact),"masters":masters,"checked_conditions":conditions}));
 }
 let output=json!({"schema":1,"scope":"Saved formal master channels evaluated only through existing native Gaussian and measure owners; no coefficient inference or reference read","input_identity":prepared.identity(),"cut_slots":[0,3],"dimension":"13/2","routing_determinant":family.routing_determinant().to_string(),"target_phases":"already present in saved coefficients; no extra Wick or per-cut sign","profiles":profiles});
 use std::io::Write;let mut f=fs::OpenOptions::new().write(true).create_new(true).open(&args[3]).map_err(|e|Error::InvalidInput(e.to_string()))?;f.write_all(serde_json::to_string_pretty(&output).unwrap().as_bytes()).map_err(|e|Error::InvalidInput(e.to_string()))?;Ok(())
}
