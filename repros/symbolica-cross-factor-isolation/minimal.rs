use symbolica::prelude::*;
use std::time::Instant;
fn main(){
 let text=include_str!("polynomial.txt").trim();
 let a=Atom::parse(text,"np_w_long",Default::default()).unwrap();
 let f:RationalPolynomial<IntegerRing,u16>=a.try_to_rational_polynomial(&Q,&Z,None).unwrap();
 assert!(f.denominator.is_constant());
 let begin=Instant::now();
 if std::env::args().any(|a|a=="--factors") {
  let factors=f.numerator.factor();
  eprintln!("factorization {:.6}s count={}",begin.elapsed().as_secs_f64(),factors.len());
  for (i,(factor,multiplicity)) in factors.into_iter().enumerate(){
   if factor.is_constant(){continue;}
   eprintln!("factor {i} multiplicity={multiplicity} polynomial={}",factor.to_expression().to_canonical_string());
   let p=factor.to_univariate_from_univariate(0).map_coeff(|c|Rational::from(c.clone()),Q);
   let start=Instant::now();let roots=p.isolate_roots();
   eprintln!("factor {i} roots={} time={:.6}s",roots.len(),start.elapsed().as_secs_f64());
  }
 }else{
  let p=f.numerator.to_univariate_from_univariate(0).map_coeff(|c|Rational::from(c.clone()),Q);
  eprintln!("whole degree={} coefficients={:?}",p.degree(),p.coefficients());
  let roots=p.isolate_roots();eprintln!("whole roots={} time={:.6}s",roots.len(),begin.elapsed().as_secs_f64());
 }
}
