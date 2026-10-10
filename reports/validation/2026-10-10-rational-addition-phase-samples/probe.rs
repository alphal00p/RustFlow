use rustred::algebra::Coefficient;
use rustred::persistence::{CoefficientId,DecodedCoefficientTable};
use symbolica::prelude::*;
use symbolica::poly::gcd::PolynomialGCD;
use std::{borrow::Cow,collections::BTreeMap,time::Instant};
#[derive(Default)]struct Metric{calls:usize,seconds:f64}
#[derive(Default)]struct Stats{phases:BTreeMap<&'static str,Metric>,equal_denominator:usize,unequal_denominator:usize,denominator_one:usize,validation_seconds:f64,exact_comparisons:usize}
macro_rules! phase{($s:expr,$n:expr,$b:expr)=>{{let start=Instant::now();let value=$b;let seconds=start.elapsed().as_secs_f64();let m=$s.phases.entry($n).or_default();m.calls+=1;m.seconds+=seconds;value}}}
impl Stats{fn json(&self)->serde_json::Value{serde_json::json!({"phases":self.phases.iter().map(|(n,m)|(*n,serde_json::json!({"calls":m.calls,"seconds":m.seconds}))).collect::<BTreeMap<_,_>>(),"equal_denominator":self.equal_denominator,"unequal_denominator":self.unequal_denominator,"denominator_one":self.denominator_one,"validation_seconds":self.validation_seconds,"exact_native_add_comparisons":self.exact_comparisons})}}
// Literal arithmetic path of pinned Symbolica ed2374f RationalPolynomial::Add,
// with IntegerRing from_num_den(..., true) expanded solely to time its phases.
// No unreduced coefficient is returned or supplied to another arithmetic API.
fn add_phased(left:&Coefficient,right:&Coefficient,stats:&mut Stats)->Coefficient{
 assert_eq!(left.get_variables(),right.get_variables());
 let result=if left.denominator==right.denominator{
  stats.equal_denominator+=1;
  let mut num=phase!(stats,"polynomial-numerator-add",&left.numerator+&right.numerator);
  let mut den=phase!(stats,"denominator-clone",left.denominator.clone());
  phase!(stats,"constructor-variable-unification",num.unify_variables(&mut den));
  if den.is_one(){stats.denominator_one+=1;}else{
   let gcd=phase!(stats,"numerator-cancellation-gcd",num.gcd(&den));
   if !gcd.is_one(){num=phase!(stats,"numerator-cancellation-quotient",num/&gcd);den=phase!(stats,"denominator-cancellation-quotient",den/&gcd);}
   phase!(stats,"constructor-sign-normalization",{if den.lcoeff().is_negative(){num=-num;den=-den;}});
  }
  Coefficient{numerator:num,denominator:den}
 }else{
  stats.unequal_denominator+=1;
  let denom_gcd=phase!(stats,"denominator-gcd",left.denominator.gcd(&right.denominator));
  let mut a_red=Cow::Borrowed(&left.denominator);let mut b_red=Cow::Borrowed(&right.denominator);
  if !denom_gcd.is_one(){a_red=Cow::Owned(phase!(stats,"denominator-gcd-quotient",&left.denominator/&denom_gcd));b_red=Cow::Owned(phase!(stats,"denominator-gcd-quotient",&right.denominator/&denom_gcd));}
  let num1=phase!(stats,"polynomial-numerator-cross-product",&left.numerator*&b_red);
  let num2=phase!(stats,"polynomial-numerator-cross-product",&right.numerator*&a_red);
  let mut num=phase!(stats,"polynomial-numerator-add",num1+num2);
  let mut den=phase!(stats,"polynomial-denominator-product",if left.denominator.nterms()>right.denominator.nterms()&&left.denominator.nterms()>a_red.nterms(){b_red.as_ref()*&left.denominator}else{a_red.as_ref()*&right.denominator});
  let gcd=phase!(stats,"numerator-cancellation-gcd",num.gcd(&denom_gcd));
  if !gcd.is_one(){num=phase!(stats,"numerator-cancellation-quotient",num/&gcd);den=phase!(stats,"denominator-cancellation-quotient",den/&gcd);}
  Coefficient{numerator:num,denominator:den}
 };
 let start=Instant::now();let native=left+right;stats.validation_seconds+=start.elapsed().as_secs_f64();assert_eq!(result,native);assert_eq!(result.get_variables(),native.get_variables());stats.exact_comparisons+=1;result
}
fn sum_phased(mut values:Vec<Coefficient>,stats:&mut Stats)->Option<Coefficient>{values.retain(|c|!c.is_zero());while values.len()>1{let mut next=Vec::with_capacity(values.len().div_ceil(2));let mut it=values.into_iter();while let Some(a)=it.next(){if let Some(b)=it.next(){let sum=add_phased(&a,&b,stats);if !sum.is_zero(){next.push(sum);}}else{next.push(a);}}values=next;}values.pop()}
fn main(){let args=std::env::args().collect::<Vec<_>>();let input=std::path::PathBuf::from(&args[1]);let out=std::path::PathBuf::from(&args[2]);let metadata:serde_json::Value=serde_json::from_slice(&std::fs::read(input.join("samples.json")).unwrap()).unwrap();let state=std::fs::read(input.join("symbolica-state.bin")).unwrap();let bytes=std::fs::read(input.join("coefficients.bin")).unwrap();let table=DecodedCoefficientTable::import_generated_normalized(&state,&bytes,Default::default()).unwrap();assert!(metadata["records"].as_array().unwrap().len()<=8);let mut records=Vec::new();
 for(s,r)in metadata["records"].as_array().unwrap().iter().enumerate(){let values=r["ids"].as_array().unwrap().iter().map(|i|table.coefficient(CoefficientId::try_from_index(i.as_u64().unwrap()as usize).unwrap()).unwrap().clone()).collect::<Vec<_>>();assert!(values.len()<=32);let expected=r["expected"].as_u64().map(|i|table.coefficient(CoefficientId::try_from_index(i as usize).unwrap()).unwrap().clone());let mut stats=Stats::default();let start=Instant::now();let value=sum_phased(values,&mut stats);let wall=start.elapsed().as_secs_f64();assert_eq!(value,expected);records.push(serde_json::json!({"sample":s,"ordinal":r["ordinal"],"exact_captured_sum":true,"wall_seconds_including_native_validation":wall,"stats":stats.json()}));std::fs::write(&out,serde_json::to_vec_pretty(&serde_json::json!({"scope":"pinned native rational Add phase attribution on eight existing samples, no new integral reduction","completed_samples":records.len(),"records":records})).unwrap()).unwrap();println!("sample {s} phase check complete");}}
