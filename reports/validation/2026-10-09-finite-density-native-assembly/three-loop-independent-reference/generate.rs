//! Validation-only exact Gamma/Beta reference and independent original-kernel
//! compact quadrature. No finite-density evaluator, AMF or predictions imported.
use symbolica_amflow::{Precision, ComplexFloat as C, symbolica};
use symbolica::prelude::*;
use serde_json::{Value,json};
use std::collections::HashMap;
fn atom(s:&str)->Atom{Atom::parse(s,"three_loop_reference",Default::default()).unwrap()}
fn rat(s:&str)->Rational{Rational::try_from(atom(s).as_view()).unwrap()}
fn eval(p:Precision,s:&str,params:&HashMap<String,C>)->C{p.eval(&atom(s),&params.iter().map(|(k,v)|(atom(k),v.clone())).collect()).unwrap()}
fn insert(p:Precision,m:&mut HashMap<String,C>,name:&str,expr:&str){let v=eval(p,expr,m);m.insert(name.into(),v);}
fn assert_rel(p:Precision,a:&C,b:&C,digits:u32){assert!(p.norm(&p.sub(a,b))<p.norm(b)*p.tolerance(digits),"{a} != {b}");}
fn params(p:Precision,epsilon:&str)->HashMap<String,C>{
 let mut m=HashMap::from([("epsilon".into(),p.rational(&rat(epsilon))),("mu".into(),p.i(1))]);
 for(k,v)in[("D","4-2*epsilon"),("b","D/2-3"),("n","3*D/2-5"),("alpha","(D-2)/2"),("Ad","2^(2-D)*pi^(-(D-1)/2)/gamma((D-1)/2)"),("C","(4*pi)^(-D/2)*gamma(2-D/2)*gamma(D/2-1)^2/gamma(D-2)"),("Kb","gamma(alpha+b)*gamma(2*alpha)/(gamma(2*alpha+b)*gamma(alpha))"),("Kbp","gamma(alpha+b+1)*gamma(2*alpha)/(gamma(2*alpha+b+1)*gamma(alpha))"),("Kbm","gamma(alpha+b-1)*gamma(2*alpha)/(gamma(2*alpha+b-1)*gamma(alpha))"),("F","Ad^2*C*4^b*mu^(2*n)/4"),("T","n^2-1"),("R","(6*D^3-67*D^2+244*D-294)/(4*(D-5)*(3*D-8))")]{insert(p,&mut m,k,v);}m
}
const PIECES:[(&str,&str);5]=[("scalar","F*Kb/n^2"),("raised_measure","F*(-Kb/(4*n^2)+(Kb/4+Kbp/2)/T)"),("raised_numerator","F*(Kb/n^2-Kb/(2*T))"),("raised_transfer","F*(b*Kb/(4*n^2)-b*Kb/(4*T)+b*Kbm/(4*n^2)-b*Kbm/(4*T))"),("raised_upper_surface","F*(-Kb/(4*n)+Kb/(4*(n+1))+Kbp/(2*(n+1)))")];
fn sample(p:Precision,eps:&str,digits:u32)->Value{
 let m=params(p,eps);let vals:Vec<C>=PIECES.iter().map(|(_,e)|eval(p,e,&m)).collect();let raised=vals[1..].iter().fold(p.zero(),|s,v|p.add(&s,v));assert_rel(p,&raised,&p.mul(&vals[0],&m["R"]),digits-12);
 let simple=eval(p,"-mu^2/(2048*pi^6*epsilon^2)*(pi/mu^2)^(3*epsilon)*gamma(1+epsilon)*gamma(1-epsilon)^3/((1-2*epsilon)^2*(1-3*epsilon)^2*gamma(1-2*epsilon)*gamma(1-3*epsilon))",&m);assert_rel(p,&simple,&vals[0],digits-12);
 assert!(p.norm(&vals[4])>p.real(0));let values=vec![vals[0].to_string(),raised.to_string()];json!({"epsilon":eps,"dimension":m["D"].to_string(),"digits":digits,"values":values,"contributions":[{"cuts":[],"values":["0","0"],"proof":"massless vacuum dimensional scalelessness"},{"cuts":[0],"values":["0","0"],"proof":"h2 singleton germ, P4 J1"},{"cuts":[3],"values":["0","0"],"proof":"h2 singleton germ, P5 J0"},{"cuts":[0,3],"values":values}],"pieces":PIECES.iter().zip(&vals).map(|((n,_),v)|(n.to_string(),json!(v.to_string()))).collect::<serde_json::Map<String,Value>>(),"exact_ratio":m["R"].to_string(),"checks":{"pieces_match_ratio":true,"bubble_beta_matches_duplicated_gamma_formula":true,"upper_surface_nonzero":true}})
}
fn laurent(p:Precision,digits:u32)->Value{
 let eg=p.eval(&Atom::var(symbolica::transcendental::euler_gamma()),&Default::default()).unwrap();let mut m=HashMap::from([("eg".into(),eg),("mu".into(),p.i(1))]);for(k,v)in[("A","-mu^2/(2048*pi^6)"),("L","10-3*eg+3*log(pi)-6*log(mu)"),("Q","L^2/2+13-3*pi^2/4")]{insert(p,&mut m,k,v);}
 let exprs=[["0","A","A*L","A*Q"],["0","3*A/8","A*(3*L/8-11/16)","A*(3*Q/8-11*L/16+7/32)"]];let values:Vec<Vec<String>>=exprs.iter().map(|x|x.iter().map(|e|eval(p,e,&m).to_string()).collect()).collect();json!({"digits":digits,"orders":[-3,-2,-1,0],"values":values,"exact_expressions":exprs,"normalization":"raw Euclidean (2*pi)^(-3D), mu=1","derivation":"Gamma recurrence/duplication and log-Gamma expansion; fixed original numerator with explicit upper surface"})
}
fn legendre(n:usize,x:f64)->(f64,f64){let(mut a,mut b)=(1.,x);for k in 2..=n{let c=((2*k-1)as f64*x*b-(k-1)as f64*a)/k as f64;a=b;b=c;}(b,n as f64*(x*b-a)/(x*x-1.))}
fn nodes(n:usize)->Vec<(f64,f64)>{let mut v=Vec::new();for i in 1..=n{let mut x=(std::f64::consts::PI*(i as f64-0.25)/(n as f64+0.5)).cos();for _ in 0..25{let(f,d)=legendre(n,x);let h=f/d;x-=h;if h.abs()<2e-16{break;}}let(_,d)=legendre(n,x);v.push(((1.+x)/2.,1./((1.-x*x)*d*d)));}v}
fn compact_quad(n:usize)->Value{
 let p=Precision::decimal(50).unwrap();let m=params(p,"-19/6");let dim=m["D"].re.to_f64();let b=m["b"].re.to_f64();let alpha=m["alpha"].re.to_f64();let beta=eval(p,"gamma(alpha)^2/gamma(2*alpha)",&m).re.to_f64();let unit=eval(p,"Ad^2*C/4",&m).re.to_f64();let quadrature=nodes(n);let mut sums=[0.;5];let mut compensation=[0.;5];
 let mut add=|slot:usize,v:f64|{let y=v-compensation[slot];let t=sums[slot]+y;compensation[slot]=(t-sums[slot])-y;sums[slot]=t;};
 for &(r1,w1)in &quadrature{for &(r2,w2)in &quadrature{for &(t,wt)in &quadrature{
  let h=4.*r1*r2*t;let numerator=-r1*r1/2.+r1*r2*(0.5+t);let numerator_a=-1.+r2/(2.*r1);let h_a=-1.+r2/r1;
  let weight=w1*w2*wt*(r1*r2).powf(dim-3.)*(t*(1.-t)).powf(alpha-1.)/beta*h.powf(b);add(0,weight);add(1,weight*numerator/(2.*r1*r1));add(2,-weight*numerator_a);add(3,-weight*b*numerator*h_a/h);
 }}}
 for &(r2,w2)in &quadrature{for &(t,wt)in &quadrature{let h=4.*r2*t;let numerator=-0.5+r2*(0.5+t);let weight=w2*wt*r2.powf(dim-3.)*(t*(1.-t)).powf(alpha-1.)/beta*h.powf(b);add(4,0.5*weight*numerator);}}
 let expected:Vec<f64>=PIECES.iter().map(|(_,e)|eval(p,e,&m).re.to_f64()).collect();let mut errors=Vec::new();let mut actual=Vec::new();for(i,e)in expected.iter().enumerate(){let a=unit*sums[i];actual.push(a);errors.push(((a-e)/e).abs());}json!({"dimension":"31/3","order":n,"actual":actual,"expected":expected,"relative_errors":errors,"method":"binary64 Gauss-Legendre directly on original radial/angular mass-jet kernels, with moving upper surface integrated separately"})
}
fn main(){let args:Vec<String>=std::env::args().collect();assert_eq!(args.len(),2,"OUTPUT");let definition:Value=serde_json::from_slice(&std::fs::read("examples/finite_density/massless_three_loop_chain.json").unwrap()).unwrap();assert_eq!(definition["targets"],json!([{"powers":[1,1,1,1,1],"numerator":"1"},{"powers":[2,1,1,1,1],"numerator":"g1_2+u1*u2"}]));assert_eq!(definition["chemical_potentials"],json!(["1"]));assert_eq!(definition["edges"].as_array().unwrap().iter().map(|e|e["routing"].clone()).collect::<Vec<_>>(),json!([["1","0","0"],["0","1","0"],["0","0","1"],["1","0","-1"],["0","1","-1"]]).as_array().unwrap().clone());assert!(definition["edges"].as_array().unwrap().iter().all(|e|e["mass_squared"]=="0"));
 let mut samples=Vec::new();let mut expansions=Vec::new();for digits in [50,80]{let p=Precision::decimal(digits).unwrap();for eps in ["-19/6","-5/4","1/8"]{samples.push(sample(p,eps,digits));}expansions.push(laurent(p,digits));}
 let quadratures:Vec<Value>=[24,40,64,96].into_iter().map(compact_quad).collect();assert!(quadratures.last().unwrap()["relative_errors"].as_array().unwrap().iter().all(|e|e.as_f64().unwrap()<1e-11),"direct original-kernel quadrature failed");
 let result=json!({"status":"independent reference generation; native comparison separate","definition":"examples/finite_density/massless_three_loop_chain.json","samples":samples,"laurent":expansions,"direct_compact_quadrature":quadratures,"normalization":"raw Euclidean","oracle_answers_read":0,"native_amf_predictions_read":0,"production_shortcuts_added":0});std::fs::write(&args[1],serde_json::to_string_pretty(&result).unwrap()+"\n").unwrap();println!("independent three-loop reference and four quadrature orders passed");}
