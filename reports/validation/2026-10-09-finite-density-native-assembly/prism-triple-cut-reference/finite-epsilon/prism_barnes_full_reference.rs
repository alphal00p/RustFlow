//! Validation-only Barnes quadrature. Uses native arithmetic and Gamma only;
//! imports no finite-density flow, AMF reduction, boundary or prediction.
use symbolica_amflow::{Precision, ComplexFloat as C, symbolica};
use symbolica::prelude::*;
use serde_json::{Value,json};
use std::collections::HashMap;
use std::time::Instant;
fn r(s:&str)->Rational{Rational::try_from(Atom::parse(s,"prism_reference",Default::default()).unwrap().as_view()).unwrap()}
fn q(v:&Value)->Rational{r(v.as_str().unwrap())}
fn conj(c:&C)->C{C::new(c.re.clone(),-c.im.clone())}
struct GammaCache{p:Precision,step:C,argument:Atom,expression:Atom,cache:HashMap<(String,i64),C>,evaluations:usize}
impl GammaCache{
 fn new(p:Precision,step:C)->Self{let a=Atom::parse("gamma_arg","prism_reference",Default::default()).unwrap();let e=symbolica::transcendental::gamma().call(&a);Self{p,step,argument:a,expression:e,cache:HashMap::new(),evaluations:0}}
 fn gamma(&mut self,offset:&Rational,imag_index:i64)->C{
  let p=self.p;let mut n=offset.to_multi_prec_float(53).to_f64().floor() as i64;
  let mut base=offset.clone()-Rational::from(n);
  while base<Rational::from(0){n-=1;base+=Rational::from(1);}
  while base>=Rational::from(1){n+=1;base-=Rational::from(1);}
  if base==Rational::from(0){base=Rational::from(1);n-=1;}
  let key=(base.to_string(),imag_index.abs());
  let base_value=if let Some(v)=self.cache.get(&key){v.clone()}else{
   let z=C::new(p.rational(&base).re,p.scale(&self.step,imag_index.abs(),1).re);
   let value=if imag_index==0{p.gamma_real(&z.re).unwrap()}else{p.eval(&self.expression,&[(self.argument.clone(),z)].into_iter().collect()).unwrap()};self.evaluations+=1;self.cache.insert(key,value.clone());value
  };
  let mut value=if imag_index<0{conj(&base_value)}else{base_value};
  let z=C::new(p.rational(&base).re,p.scale(&self.step,imag_index,1).re);
  if n>=0{for j in 0..n{value=p.mul(&value,&p.add(&z,&p.i(j)));}}
  else{for j in 1..=-n{value=p.div(&value,&p.sub(&z,&p.i(j)));}}
  value
 }
}
#[derive(Clone)]struct Factor{offset:Rational,slope:i64,sign:i64}
fn kernel(cache:&mut GammaCache,factors:&[Factor],node:i64)->C{
 let p=cache.p;let mut out=p.i(1);
 for f in factors{let value=cache.gamma(&f.offset,node*f.slope);out=if f.sign>0{p.mul(&out,&value)}else{p.div(&out,&value)};}
 out
}
fn integer_pole(a:&Rational)->bool{if a>&Rational::from(0){return false;}let n=a.to_multi_prec_float(53).to_f64().round() as i64;*a==Rational::from(n)}
fn leaf(cache:&mut GammaCache,source:&Value,n:i64,pi:&C)->C{
 let p=cache.p;let contour=source["contour"].as_object().unwrap();let mut axes:Vec<usize>=contour.keys().map(|s|s.parse().unwrap()).collect();axes.sort();
 let mut groups:HashMap<Vec<i64>,Vec<Factor>>=HashMap::new();let mut constants=Vec::new();
 for (name,sign)in[("num",1),("den",-1)]{for g in source[name].as_array().unwrap(){let g:Vec<Rational>=g.as_array().unwrap().iter().map(q).collect();let mut offset=g[0].clone()+g[1].clone();for &a in &axes{offset+=g[a].clone()*q(&contour[&a.to_string()]);}let coefficients:Vec<i64>=axes.iter().map(|&a|g[a].to_multi_prec_float(53).to_f64() as i64).collect();if coefficients.iter().all(|&x|x==0){if integer_pole(&offset){assert_eq!(sign,-1,"uncancelled numerator pole");return p.zero();}constants.push(Factor{offset,slope:0,sign});}else{let slope=*coefficients.iter().find(|&&x|x!=0).unwrap();assert!(slope.abs()==1);let key=coefficients.iter().map(|x|x/slope).collect();groups.entry(key).or_default().push(Factor{offset,slope,sign});}}}
 let constant=p.mul(&p.rational(&q(&source["coefficient"])),&kernel(cache,&constants,0));
 if axes.is_empty(){return constant;}
 let values=|cache:&mut GammaCache,fs:&[Factor],m:i64|->Vec<C>{(-m..=m).map(|j|kernel(cache,fs,j)).collect()};
 let weight=p.div(&cache.step,&p.scale(pi,2,1));
 if axes.len()==1{let fs=groups.remove(&vec![1]).unwrap_or_default();let vs=values(cache,&fs,n);let sum=vs.iter().fold(p.zero(),|s,x|p.add(&s,x));return p.mul(&constant,&p.mul(&sum,&weight));}
 assert_eq!(axes.len(),2);assert!(groups.keys().all(|k|[vec![1,0],vec![0,1],vec![1,1]].contains(k)));
 let a=values(cache,&groups.remove(&vec![1,0]).unwrap_or_default(),n);let b=values(cache,&groups.remove(&vec![0,1]).unwrap_or_default(),n);let c=values(cache,&groups.remove(&vec![1,1]).unwrap_or_default(),2*n);
 let mut sum=p.zero();for(i,av)in a.iter().enumerate(){let mut row=p.zero();for(j,bv)in b.iter().enumerate(){row=p.add(&row,&p.mul(bv,&c[i+j]));}sum=p.add(&sum,&p.mul(av,&row));}
 p.mul(&constant,&p.mul(&sum,&p.powi(&weight,2)))
}
fn outside(cache:&mut GammaCache,row:&Value)->C{let p=cache.p;let mut out=p.i(1);for(name,sign)in[("outside_num",1),("outside_den",-1)]{for g in row[name].as_array().unwrap(){let a=q(&g[0])+q(&g[1]);if integer_pole(&a){assert_eq!(sign,-1);return p.zero();}let v=cache.gamma(&a,0);out=if sign>0{p.mul(&out,&v)}else{p.div(&out,&v)};}}p.mul(&out,&p.pow(&p.i(4),&p.rational(&(q(&row["four_exponent"][0])+q(&row["four_exponent"][1])))))}
fn fp_row(cache:&mut GammaCache,row:&Value,n:i64,pi:&C)->(C,C,C){
 let p=cache.p;let mut total=p.zero();let mut pole=p.zero();let mut pole_scale=p.zero();let log4=p.log(&p.i(4));let s0=q(&row["four_exponent"][0])+q(&row["four_exponent"][1]);let s1=q(&row["four_exponent"][4])+q(&row["four_exponent"][5]);let four=p.pow(&p.i(4),&p.rational(&s0));
 let eg=p.eval(&Atom::var(symbolica::transcendental::euler_gamma()),&Default::default()).unwrap();
 let arg=Atom::parse("psi_arg","prism_reference",Default::default()).unwrap();let psi=symbolica::transcendental::polygamma().call((Atom::num(0),arg.clone()));
 for old in row["leaves"].as_array().unwrap(){
  let axes:Vec<usize>=old["contour"].as_object().unwrap().keys().map(|s|s.parse().unwrap()).collect();let mut power=0i64;let mut leading=q(&old["coefficient"]);let mut logs=Vec::new();let mut pole_logs=p.zero();let mut output=json!({"coefficient":"0","contour":old["contour"],"num":[],"den":[]});let mut identically_zero=false;
  for(name,outside,sign)in[("num","outside_num",1),("den","outside_den",-1)]{
   for g in old[name].as_array().unwrap().iter().chain(row[outside].as_array().unwrap()){
    let off=q(&g[0])+q(&g[1]);let slope=q(&g[4])+q(&g[5]);let has_variables=axes.iter().any(|&a|q(&g[a])!=Rational::from(0));
    if !has_variables&&integer_pole(&off){
     if slope==Rational::from(0){assert_eq!(sign,-1,"identical numerator pole");identically_zero=true;break;}
     power+=sign;let degree=(-off.to_multi_prec_float(53).to_f64())as i64;let mut fact=Rational::from(1);let mut harmonic=Rational::from(0);for j in 1..=degree{fact*=Rational::from(j);harmonic+=Rational::from((1,j));}let pref=Rational::from(if degree%2==0{1}else{-1})/(fact*slope.clone());leading=if sign>0{leading*pref}else{leading/pref};let log=p.mul(&p.rational(&slope),&p.sub(&p.rational(&harmonic),&eg));pole_logs=p.add(&pole_logs,&p.scale(&log,sign,1));
    }else{
     output[name].as_array_mut().unwrap().push(json!([off.to_string(),"0",g[2],g[3]]));if slope!=Rational::from(0){logs.push((off,slope,sign,has_variables));}
    }
   }
   if identically_zero{break;}
  }
  if identically_zero||power<0{continue;}assert!(power<=1,"higher regulator pole requires another derivative");output["coefficient"]=json!(leading.to_string());let base=leaf(cache,&output,n,pi);
  if power==0{total=p.add(&total,&base);}else{
   assert!(axes.is_empty(),"nonconstant regulator derivative not implemented");let mut log=p.add(&pole_logs,&p.mul(&p.rational(&s1),&log4));
   for(off,slope,sign,hasvar)in logs{assert!(!hasvar);let value=p.eval(&psi,&[(arg.clone(),p.rational(&off))].into_iter().collect()).unwrap();log=p.add(&log,&p.scale(&p.mul(&p.rational(&slope),&value),sign,1));}
   total=p.add(&total,&p.mul(&base,&log));pole=p.add(&pole,&base);pole_scale=p.add(&pole_scale,&C::new(p.norm(&base),p.real(0)));
  }
 }
 (p.mul(&four,&total),p.mul(&four,&pole),p.mul(&four,&pole_scale))
}

fn two_cut_reference(p:Precision,epsilon:&str)->Vec<C>{
 let parse=|s:&str|Atom::parse(s,"prism_two_cut_reference",Default::default()).unwrap();
 let mut parameters=HashMap::<String,C>::from([("eps".into(),p.rational(&r(epsilon)))]);
 let eval=|s:&str,values:&HashMap<String,C>|p.eval(&parse(s),&values.iter().map(|(k,v)|(parse(k),v.clone())).collect()).unwrap();
 let dimension=p.sub(&p.i(4),&p.scale(&parameters["eps"],2,1));parameters.insert("D".into(),dimension);
 let masters=vec![eval("gamma(eps)*gamma(1-eps)^2/gamma(2-2*eps)*gamma(2*eps)*gamma(1-2*eps)^2/gamma(2-3*eps)",&parameters),eval("-gamma(1-eps)^3*gamma(2*eps-1)/gamma(3-3*eps)",&parameters),eval("(gamma(eps)*gamma(1-eps)^2/gamma(2-2*eps))^2",&parameters)];
 let candidates:Value=serde_json::from_slice(&std::fs::read("fixtures/finite_density/prism_virtual_candidate_relations.json").unwrap()).unwrap();assert_eq!(candidates["master_indices"],json!([[0,1,0,1,1,1,0],[0,1,1,0,0,1,0],[1,1,1,1,0,0,0]]));let virtual_factor=eval("(4*pi)^(-D)",&parameters);
 for relation in candidates["relations"].as_array().unwrap(){let label=relation["label"].as_str().unwrap();let coefficients=relation["candidate_coefficients"].as_array().unwrap();assert_eq!(coefficients.len(),3);let mut value=p.zero();for(i,c)in coefficients.iter().enumerate(){value=p.add(&value,&p.mul(&eval(c.as_str().unwrap(),&parameters),&masters[i]));}value=p.mul(&virtual_factor,&value);if label=="CC"{value=p.neg(&value);}parameters.insert(label.into(),value);}
 for(name,expression)in[("Ad","2^(2-D)*pi^(-(D-1)/2)/gamma((D-1)/2)"),("alpha","(D-2)/2"),("b","D-5"),("n","2*D-7"),("Kb","gamma(alpha+b)*gamma(2*alpha)/(gamma(2*alpha+b)*gamma(alpha))"),("Kbm","gamma(alpha+b-1)*gamma(2*alpha)/(gamma(2*alpha+b-1)*gamma(alpha))")]{let value=eval(expression,&parameters);parameters.insert(name.into(),value);}
 let mut result=Vec::new();for suffix in ["A","B"]{parameters.insert("C0".into(),parameters[&format!("C0{suffix}")].clone());parameters.insert("C1".into(),parameters[&format!("C1{suffix}")].clone());result.push(eval("Ad^2/4*4^(b-1)*(2*C0*Kb*(n-1)/(n*(n-2))+(b*C0+C1)*Kbm/(n-1)^2-b*C0*Kbm/(n*(n-2)))",&parameters));}
 result.push(eval("Ad^2/4*4^(D-6)*Kbm/(2*D-8)^2*CC",&parameters));result
}

fn main(){let args:Vec<String>=std::env::args().collect();assert_eq!(args.len(),6,"FILE DIGITS STEP TRUNCATION OUTPUT");let data:Value=serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();let digits:u32=args[2].parse().unwrap();let p=Precision::decimal(digits).unwrap();let step=p.rational(&r(&args[3]));let truncation:f64=args[4].parse().unwrap();let n=(truncation/step.re.to_f64()).ceil()as i64;let pi=p.eval(&Atom::var(Symbol::PI),&Default::default()).unwrap();let mut cache=GammaCache::new(p,step);let started=Instant::now();let mut values=Vec::new();let mut pieces=HashMap::new();
 for(i,row)in data["rows"].as_array().unwrap().iter().enumerate(){let (value,residue,residue_scale)=fp_row(&mut cache,row,n,&pi);if p.norm(&residue_scale)>p.real(0){assert!(p.norm(&residue)<p.norm(&residue_scale)*p.tolerance(digits.saturating_sub(12)),"rho pole did not cancel in row {i}: {residue} / {residue_scale}");}println!("row {i}: {value}; gamma={} elapsed={:?}",cache.evaluations,started.elapsed());let sector=row["row"]["sector"].as_str().unwrap();let entry=pieces.entry(sector.to_string()).or_insert_with(||p.zero());*entry=p.add(entry,&value);values.push(value);}
 let dim=q(&data["rows"][0]["D"][0])+q(&data["rows"][0]["D"][1]);let d=dim.clone()-Rational::from(1);let ad=p.div(&p.pow(&p.i(2),&p.rational(&(Rational::from(1)-d.clone()))),&p.mul(&p.pow(&pi,&p.rational(&(d.clone()/Rational::from(2)))),&cache.gamma(&(d/Rational::from(2)),0)));let pref=p.div(&p.powi(&ad,3),&p.scale(&p.pow(&p.scale(&pi,4,1),&p.rational(&(dim/Rational::from(2)))),8,1));let total=p.mul(&pref,&values.iter().fold(p.zero(),|s,x|p.add(&s,x)));let pieces:HashMap<String,String>=pieces.into_iter().map(|(k,v)|(k,p.mul(&pref,&v).to_string())).collect();let two_cuts=two_cut_reference(p,data["epsilon"].as_str().unwrap());let complete=p.add(&total,&two_cuts.iter().fold(p.zero(),|s,x|p.add(&s,x)));let eps=p.rational(&r(data["epsilon"].as_str().unwrap()));let eg=p.eval(&Atom::var(symbolica::transcendental::euler_gamma()),&Default::default()).unwrap();let scale=p.mul(&p.powi(&p.scale(&pi,4,1),8),&p.exp(&p.mul(&p.scale(&eps,4,1),&p.sub(&eg,&p.log(&pi)))));let scaled=p.mul(&scale,&complete);let report=json!({"status":"validation-only rho finite-part three-cut sample; not a complete reference","epsilon":data["epsilon"],"analytic_regulator":"rho removed by exact finite-part rule; all per-monomial rho residues checked","digits":digits,"step":args[3],"truncation":truncation,"nodes_per_axis":2*n+1,"value":total.to_string(),"pieces":pieces,"two_cut_values":two_cuts.iter().map(ToString::to_string).collect::<Vec<_>>(),"full_raw_value":complete.to_string(),"full_msbar_reference_normalization":scaled.to_string(),"normalization":"raw times (4*pi)^8*(exp(EulerGamma)/pi)^(4*eps), mu=1","remaining_sector_prescription":"vacuum and singleton sectors dimensionally zero under the common massless regulated continuation; final proof binding reviewed separately","rows":values.iter().map(ToString::to_string).collect::<Vec<_>>(),"gamma_evaluations":cache.evaluations,"wall_seconds":started.elapsed().as_secs_f64(),"backend":"native Precision and Symbolica Gamma; no finite-density or AMF evaluation"});std::fs::write(&args[5],serde_json::to_string_pretty(&report).unwrap()+"\n").unwrap();println!("three-cut {total}; complete raw {complete}; normalized {scaled}");}
