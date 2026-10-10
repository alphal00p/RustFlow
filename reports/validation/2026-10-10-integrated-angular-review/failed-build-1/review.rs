//! Independent Cartesian sphere-coordinate expansion, not Wick pair enumeration.
#[path="/common/dev/rustflow_fermi/reports/validation/2026-10-10-integrated-moment-exploration/angular-aggregate/angular.rs"]
mod angular;
use std::collections::BTreeMap;
use symbolica::{prelude::*,domains::{Ring,RingOps}};
fn sphere(exponents:&[u16],d:usize)->Rational {
    if exponents.iter().any(|e|e%2!=0){return Q.zero();}
    let degree:usize=exponents.iter().map(|e|usize::from(*e)).sum();
    let numerator=exponents.iter().fold(1i64,|p,e|p*(1..usize::from(*e)).step_by(2).map(|k|k as i64).product::<i64>());
    let denominator=(0..degree/2).map(|k|(d+2*k)as i64).product::<i64>();
    Q.to_element(Integer::from(numerator),Integer::from(denominator),true)
}
fn cartesian(loops:usize,d:usize,edges:&[(usize,usize)],powers:&[u16])->Rational {
    // Expand each Gram dot product as sum_{coordinate=0}^{d-1} x_i x_j.
    // Independent sphere-coordinate moments follow by dividing the product
    // Gaussian moment by its radial moment; no graph pairings are used here.
    let mut polynomial=BTreeMap::from([(vec![0u16;loops*d],1u64)]);
    for (&(i,j),&power) in edges.iter().zip(powers) {
        for _ in 0..power {
            let mut next=BTreeMap::new();
            for (exponents,coefficient) in polynomial {
                for a in 0..d {
                    let mut exponents=exponents.clone();exponents[i*d+a]+=1;exponents[j*d+a]+=1;
                    *next.entry(exponents).or_insert(0u64)+=coefficient;
                }
            }
            polynomial=next;
        }
    }
    polynomial.iter().fold(Q.zero(),|sum,(exponents,coefficient)|{
        let v=(0..loops).fold(Rational::from(*coefficient as i64),|v,i|Q.mul(&v,&sphere(&exponents[i*d..(i+1)*d],d)));
        Q.add(&sum,&v)
    })
}
fn powers(out:&mut Vec<Vec<u16>>,at:usize,budget:u16,max:u16,vector:&mut [u16]){
    if at==vector.len(){out.push(vector.to_vec());return;}
    for x in 0..=budget.min(max){vector[at]=x;powers(out,at+1,budget-x,max,vector);}
}
#[test]
fn independent_cartesian_three_and_four_vector_moments(){
    let mut count=0;
    for (loops,max,total) in [(3,3,6),(4,2,5)] {
        let edges=(0..loops).flat_map(|i|(i+1..loops).map(move|j|(i,j))).collect::<Vec<_>>();
        let slots=(0..loops).flat_map(|i|(i..loops).map(move|j|(i,j))).collect::<Vec<_>>();
        let mut cases=Vec::new();powers(&mut cases,0,total,max,&mut vec![0;edges.len()]);
        for d in [1,2,3,5] {
            let mut averager=angular::AngularAverager::new(loops,Atom::num(d as i64),200000).unwrap();
            for case in &cases {
                let mut full=vec![0;slots.len()];
                for (edge,p) in edges.iter().zip(case){full[slots.iter().position(|x|x==edge).unwrap()]=*p;}
                let result=averager.monomial(&full).unwrap().1;
                let expected=Atom::num(cartesian(loops,d,&edges,case));
                assert!((result-expected).together().cancel().is_zero(),"loops={loops} d={d} powers={case:?}");count+=1;
            }
        }
    }
    println!("independent Cartesian cases: {count}");
}
#[test]
fn diagonal_powers_and_loop_permutations(){
    let loops=4;let slots=(0..loops).flat_map(|i|(i..loops).map(move|j|(i,j))).collect::<Vec<_>>();
    let original=vec![3,2,1,1,4,1,1,5,0,6];
    let mut ave=angular::AngularAverager::new(loops,parse!("11/2"),200000).unwrap();
    let(base_rad,base)=ave.monomial(&original).unwrap();assert!(!base.is_zero());
    for p in [[0,1,2,3],[3,2,1,0],[1,3,0,2],[2,0,3,1]] {
        let mut transformed=vec![0;slots.len()];
        for (s,&(i,j)) in slots.iter().enumerate(){let(a,b)=(p[i].min(p[j]),p[i].max(p[j]));transformed[slots.iter().position(|x|*x==(a,b)).unwrap()]=original[s];}
        let(rad,c)=ave.monomial(&transformed).unwrap();assert!((c-&base).together().cancel().is_zero());
        for i in 0..loops {assert_eq!(rad[p[i]],base_rad[i]);}
    }
}
