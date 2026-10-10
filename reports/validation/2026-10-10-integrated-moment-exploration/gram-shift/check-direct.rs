//! Dimension shift tested directly through native Gaussian terminals.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{gaussian, *};
type Poly=BTreeMap<Vec<u16>,i64>;
fn determinant(h:usize)->Poly {
    fn perm(p:&mut Vec<usize>,h:usize,out:&mut Vec<Vec<usize>>) {
        if p.len()==h {out.push(p.clone());return;}
        for j in 0..h {if !p.contains(&j) {p.push(j);perm(p,h,out);p.pop();}}
    }
    let slots=(0..h).flat_map(|i|(i..h).map(move|j|(i,j))).collect::<Vec<_>>();
    let mut permutations=vec![];perm(&mut vec![],h,&mut permutations);
    let mut result=Poly::new();
    for p in permutations {
        let inversions=(0..h).map(|i|((i+1)..h).filter(|&j|p[i]>p[j]).count()).sum::<usize>();
        let mut powers=vec![0;slots.len()];
        for (i,&j) in p.iter().enumerate() {powers[slots.iter().position(|&p|p==(i.min(j),i.max(j))).unwrap()]+=1;}
        *result.entry(powers).or_insert(0)+=if inversions%2==0 {1}else{-1};
    }
    result.retain(|_,v|*v!=0);result
}
fn multiply(a:&Poly,b:&Poly)->Poly {
    let mut result=Poly::new();
    for (pa,ca) in a {for(pb,cb)in b {
        *result.entry(pa.iter().zip(pb).map(|(a,b)|a+b).collect()).or_insert(0)+=ca*cb;
    }}
    result.retain(|_,v|*v!=0);result
}
fn integrate(poly:&Poly,h:usize,dimension:i64,eps:&Rational,p:Precision)->Result<ComplexFloat> {
    let slots=(0..h).flat_map(|i|(i..h).map(move|j|(i,j))).collect::<Vec<_>>();
    let mut propagators=vec![Propagator{constant:Atom::num(-3),scalar_products:slots.iter().map(|&(i,j)|Atom::num(if i==j{3}else{2})).collect()}];
    for s in 0..slots.len() {propagators.push(Propagator{constant:Atom::zero(),scalar_products:(0..slots.len()).map(|j|Atom::num(if s==j{1}else{0})).collect()});}
    let family=IntegralFamily{name:"gram_shift_direct_gaussian_check".into(),loops:(0..h).map(|i|format!("k{i}")).collect(),external:vec![],external_gram:vec![],propagators,physical_propagators:1,epsilon:symbol!("gram_shift_check::epsilon"),dimension};
    let mut sum=p.zero();
    for(powers,c)in poly {
        let mut indices=vec![40];indices.extend(powers.iter().map(|&k|-(k as i16)));
        sum=p.add(&sum,&p.mul(&p.i(*c),&gaussian::terminal(&family,&Integral(indices),eps,p)?));
    }
    Ok(sum)
}
fn main()->Result<()> {
    let mut results=vec![];
    for digits in [50,80] {
        let p=Precision::decimal(digits)?;
        for(h,n)in[(1usize,1usize),(1,2),(1,3),(1,4),(2,1),(2,2),(3,1),(4,1)] {
            eprintln!("start loops={h} shifts={n} digits={digits}");
            let eps=Rational::from((1,7));let d=Rational::from(6)-&eps*&Rational::from(2);
            let one=Poly::from([(vec![0;h*(h+1)/2],1)]);let det=determinant(h);
            let mut insertion=one.clone();let mut factor=Rational::one();
            for r in 0..n {insertion=multiply(&insertion,&det);for j in 0..h {
                factor*=Rational::from(-2)/(d.clone()+Rational::from((2*r)as i64)-Rational::from(j as i64));
            }}
            let shifted=integrate(&one,h,6+(2*n)as i64,&eps,p)?;
            let inserted=integrate(&insertion,h,6,&eps,p)?;
            let bridged=p.mul(&p.rational(&factor),&inserted);
            let error=p.norm(&p.div(&p.sub(&bridged,&shifted),&shifted));
            let passed=error<p.parse("1e-35","0")?.re;
            results.push(serde_json::json!({"loops":h,"shift_count":n,"digits":digits,"base_dimension":"40/7","quadratic":"A=2I+11^T; mass_squared=3; denominator_power=40","native_bridge_factor":factor.to_string(),"gram_monomials":insertion.len(),"shifted_value":shifted.re.to_string(),"relative_error":error.to_string(),"passed":passed}));
            std::fs::write(std::env::args().nth(1).expect("output path"),serde_json::to_vec_pretty(&serde_json::json!({"scope":"Native Gaussian identity check, no physical finite-density prediction","cases":results,"complete":false})).unwrap()).unwrap();
            eprintln!("done passed={passed} relative_error={error}");
        }
    }
    let passed=results.iter().all(|r|r["passed"].as_bool()==Some(true));
    let output=serde_json::json!({"scope":"Generic vacuum dimension-shift normalization checks, no physical finite-density prediction","passed":passed,"complete":true,"cases":results});
    std::fs::write(std::env::args().nth(1).expect("output path"),serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    println!("checks={} passed={passed}",output["cases"].as_array().unwrap().len());assert!(passed);Ok(())
}
