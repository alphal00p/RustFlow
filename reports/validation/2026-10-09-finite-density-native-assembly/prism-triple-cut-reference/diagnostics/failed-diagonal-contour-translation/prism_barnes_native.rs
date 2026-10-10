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
fn main(){let args:Vec<String>=std::env::args().collect();assert_eq!(args.len(),6,"FILE DIGITS STEP TRUNCATION OUTPUT");let data:Value=serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();let digits:u32=args[2].parse().unwrap();let p=Precision::decimal(digits).unwrap();let step=p.rational(&r(&args[3]));let truncation:f64=args[4].parse().unwrap();let n=(truncation/step.re.to_f64()).ceil()as i64;let pi=p.eval(&Atom::var(Symbol::PI),&Default::default()).unwrap();let mut cache=GammaCache::new(p,step);let started=Instant::now();let mut values=Vec::new();let mut pieces=HashMap::new();
 for(i,row)in data["rows"].as_array().unwrap().iter().enumerate(){let mut sum=p.zero();for l in row["leaves"].as_array().unwrap(){sum=p.add(&sum,&leaf(&mut cache,l,n,&pi));}let value=p.mul(&outside(&mut cache,row),&sum);println!("row {i}: {value}; gamma={} elapsed={:?}",cache.evaluations,started.elapsed());let sector=row["row"]["sector"].as_str().unwrap();let entry=pieces.entry(sector.to_string()).or_insert_with(||p.zero());*entry=p.add(entry,&value);values.push(value);}
 let dim=q(&data["rows"][0]["D"][0])+q(&data["rows"][0]["D"][1]);let d=dim.clone()-Rational::from(1);let ad=p.div(&p.pow(&p.i(2),&p.rational(&(Rational::from(1)-d.clone()))),&p.mul(&p.pow(&pi,&p.rational(&(d.clone()/Rational::from(2)))),&cache.gamma(&(d/Rational::from(2)),0)));let pref=p.div(&p.powi(&ad,3),&p.scale(&p.pow(&p.scale(&pi,4,1),&p.rational(&(dim/Rational::from(2)))),8,1));let total=p.mul(&pref,&values.iter().fold(p.zero(),|s,x|p.add(&s,x)));let pieces:HashMap<String,String>=pieces.into_iter().map(|(k,v)|(k,p.mul(&pref,&v).to_string())).collect();let report=json!({"status":"regulated validation-only Barnes sample; not a complete reference","epsilon":data["epsilon"],"analytic_regulator":data["analytic_regulator"],"digits":digits,"step":args[3],"truncation":truncation,"nodes_per_axis":2*n+1,"value":total.to_string(),"pieces":pieces,"rows":values.iter().map(ToString::to_string).collect::<Vec<_>>(),"gamma_evaluations":cache.evaluations,"wall_seconds":started.elapsed().as_secs_f64(),"backend":"native Precision and Symbolica Gamma; no finite-density or AMF evaluation"});std::fs::write(&args[5],serde_json::to_string_pretty(&report).unwrap()+"\n").unwrap();println!("total {total}");}
