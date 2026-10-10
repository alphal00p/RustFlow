//! Independent reference algebra only: original mass jet and Beta rational identities.
use symbolica_amflow::symbolica;
use symbolica::prelude::*;
use serde_json::json;
fn atom(s:&str)->Atom{Atom::parse(s,"three_loop_reference_exact",Default::default()).unwrap()}
fn symbol(s:&str)->Symbol{match atom(s).as_view(){AtomView::Var(v)=>v.get_symbol(),_=>unreachable!()}}
fn zero(s:&str)->String{let value=atom(s).together().cancel();assert!(value.is_zero(),"{value}");value.to_canonical_string()}
fn main(){
 let args:Vec<String>=std::env::args().collect();assert_eq!(args.len(),2);
 let numerator=atom("(-a+2*E1*E2-2*r1*r2*z-a)/4-E1^2/2+E1*E2/2");
 let mass_jet=(numerator.derivative(symbol("a"))+numerator.derivative(symbol("E1"))*atom("1/(2*E1)")-atom("-1+E2/(2*E1)")).together().cancel();assert!(mass_jet.is_zero());
 let n="((3*D-10)/2)";let b="((D-6)/2)";let kp=format!("((D-4)/{n})");let km=format!("(({n}-1)/(D-5))");let t=format!("({n}^2-1)");
 let measure=format!("(-1/(4*{n}^2)+(1/4+{kp}/2)/{t})");
 let original_numerator=format!("(1/{n}^2-1/(2*{t}))");
 let transfer=format!("({b}/(4*{n}^2)-{b}/(4*{t})+{b}*{km}/(4*{n}^2)-{b}*{km}/(4*{t}))");
 let upper=format!("(-1/(4*{n})+1/(4*({n}+1))+{kp}/(2*({n}+1)))");
 let bulk=format!("({measure}+{original_numerator}+{transfer})*{n}^2");let surface=format!("{upper}*{n}^2");
 let ratio="(6*D^3-67*D^2+244*D-294)/(4*(D-5)*(3*D-8))";
 let cases=vec![
  ("bulk_beta_moments_to_rational",zero(&format!("{bulk}-(10*D^2-81*D+156)/(8*(D-5)*({n}+1))"))),
  ("upper_beta_moments_to_rational",zero(&format!("{surface}-{n}*(2*D-9)/(4*({n}+1))"))),
  ("complete_raised_to_scalar_ratio",zero(&format!("{bulk}+{surface}-({ratio})"))),
  ("epsilon_ratio_through_second_order",zero(&format!("({})-(3/8-11*epsilon/16+7*epsilon^2/32)+epsilon^3*(42*epsilon+53)/(64*(epsilon+1/2)*(3*epsilon-2))",ratio.replace("D","(4-2*epsilon)")))),
 ];
 let result=json!({"status":"passed","scope":"independent exact reference algebra; no native AMF evaluation or oracle values","mass_jet_difference":mass_jet.to_canonical_string(),"identities":cases,"method":"native Symbolica rational cancellation of the original mass derivative and unsimplified Beta moments","oracle_values_read":0,"native_predictions_read":0});
 std::fs::write(&args[1],serde_json::to_string_pretty(&result).unwrap()+"\n").unwrap();println!("five exact reference identities passed");
}
