#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::Decode;
use rustred::persistence::{inspect_program, BinaryIoLimits, SectionTag, DecodedCoefficientTable, CoefficientId};
use rustred::solver::guarded::{GuardedSource, GuardedSourceSystem, GuardedProgram, IndexRole, IndexBounds, IndexDomain};
use rustred::solver::{Integral, Power, Term, SearchOptions};
use std::sync::Arc;
type Domain=Vec<(Option<i64>,Option<i64>)>;
type Label=Vec<(bool,i16)>;
#[derive(Decode)] struct Tr { integral: Label, coefficient:usize }
#[derive(Decode)] struct Sr { id:String, domain:Domain, conditions:Vec<usize>, terms:Vec<Tr> }
#[derive(Decode)] struct Sd { row:usize, integral:Label, shifts:Vec<i16> }
#[derive(Decode)] struct Rr { fixed:Vec<Option<i16>>,target:Label,rhs:Vec<Tr>,sources:Vec<Sd>,domain:Domain,discovery_domain:Domain,conditions:Vec<usize>,sector:Vec<bool>,permutation:Option<Vec<usize>> }
#[derive(Decode)] struct Rec { schema:String,measure:String,roles:Vec<u8>,indices:Vec<usize>,sources:Vec<Sr>,zero_domains:Vec<Domain>,rules:Vec<Rr>,terminals:Vec<Vec<i64>> }
fn domain<const N:usize>(v:&Domain)->IndexDomain<N>{IndexDomain::new(v.iter().map(|&(lo,hi)|IndexBounds::new(lo,hi).unwrap()).collect::<Vec<_>>().try_into().unwrap()).unwrap()}
fn label<const N:usize>(v:&Label)->Integral<N>{Integral::new(v.iter().map(|&(s,p)|Power::new(s,p).unwrap()).collect::<Vec<_>>().try_into().unwrap())}
fn bounds<const N:usize>(v:&IndexDomain<N>)->Domain{v.bounds().iter().map(|b|(b.lower(),b.upper())).collect()}

use symbolica::prelude::*;



use symbolica_amflow::finite_density::{DensityInput,guarded::{GuardedContext,GuardedMeasureIdentity}};
use symbolica_amflow::finite_density::guarded::GuardedIdentity;
use std::collections::BTreeMap;

fn shifted(mut row: GuardedIdentity<12>, indices: &[Symbol;12], slot: usize) -> GuardedIdentity<12> {
    let replacement=BTreeMap::from([(Atom::var(indices[slot]),Atom::var(indices[slot])+Atom::one())]);
    for term in &mut row.terms {
        term.shift[slot]=term.shift[slot].checked_add(1).unwrap();
        term.coefficient=symbolica_amflow::family::substitute(&term.coefficient,&replacement);
    }
    for c in &mut row.nonzero_conditions {*c=symbolica_amflow::family::substitute(c,&replacement)}
    let mut shift=[0i64;12];shift[slot]=1;
    row.domain=row.domain.pullback(&shift).unwrap();
    assert_eq!(row.domain.bounds()[slot],IndexBounds::new(None,Some(-1)).unwrap());
    row.id=format!("mixed-energy-polynomial-slot-{slot}/{}",row.id);
    // At the only boundary where a +2 image could be a positive completion,
    // its coefficient must vanish identically. This tests the old exponent0.
    for term in &row.terms {
        if term.shift[slot]>1 {
            let c=symbolica_amflow::family::substitute(&term.coefficient,&BTreeMap::from([(Atom::var(indices[slot]),Atom::num(-1))]));
            assert!(c.expand().together().cancel().is_zero(),"positive completion boundary image survives");
        }
    }
    row
}

fn main(){
    let args:Vec<_>=std::env::args().collect();let out=std::path::PathBuf::from(&args[4]);std::fs::create_dir_all(&out).unwrap();
    let mode=&args[5];assert!(["original","slot6","slot8","both","both-reverse","temporal-only","boost-only"].contains(&mode.as_str()));
    let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
    let(record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();
    assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v2");
    let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
    let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
    let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
    let original=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect::<Vec<_>>();
    assert_eq!(original.len(),54);let variables=original[0].row[0].coefficient.variables().clone();
    let symbols=variables.iter().map(|v|{let PolyVariable::Symbol(s)=v else{panic!()};*s}).collect::<Vec<_>>();
    let indices:[Symbol;12]=record.indices.iter().map(|&i|symbols[i]).collect::<Vec<_>>().try_into().unwrap();
    let parameters=symbols.iter().enumerate().filter(|(i,_)|!record.indices.contains(i)).map(|(_,s)|*s).collect::<Vec<_>>();
    let eta=symbol!("rustflow_occupied::eta");let epsilon=symbol!("rustflow_occupied::epsilon");assert!(parameters.contains(&eta)&&parameters.contains(&epsilon));
    let input:DensityInput=serde_json::from_slice(&std::fs::read(&args[3]).unwrap()).unwrap();let prepared=input.prepare().unwrap();
    assert!(record.measure.contains(prepared.identity()));let family=prepared.occupied_cut(&[3],4096).unwrap().at_physical_masses();
    assert_eq!(family.loops(),3);assert_eq!(family.shells().len(),1);assert_eq!(family.shells()[0].loop_index,0);
    let measure=family.deformed_measure::<12>(eta,&[0,1,2,4]).unwrap();assert_eq!(measure.roles(),&roles);
    let q=family.coordinates();let pairs=[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)];
    let energy=|i:usize|q[6+i].clone();let gram=|i:usize,j:usize|q[pairs.iter().position(|&(a,b)|(a,b)==(i.min(j),i.max(j))).unwrap()].clone();
    assert!((&measure.factors()[6]-&measure.factors()[8]-energy(0)).expand().together().cancel().is_zero());
    for s in [6,8]{assert_eq!(roles[s],IndexRole::Ordinary);assert!(s>=family.physical_slots()&&s<family.input_slots())}
    let mut admitted_bounds=*IndexDomain::for_roles(&roles).bounds();for b in &mut admitted_bounds[5..9]{*b=IndexBounds::new(None,Some(0)).unwrap()}admitted_bounds[11]=IndexBounds::fixed(0);let admitted=IndexDomain::new(admitted_bounds).unwrap();
    let mut u=pairs.iter().map(|&(i,j)|Atom::num(i32::from(i==0))*energy(j)+Atom::num(i32::from(j==0))*energy(i)).collect::<Vec<_>>();u.extend([Atom::one(),Atom::zero(),Atom::zero()]);
    let mut boost=vec![Atom::zero();6];boost.extend((0..3).map(|i|(energy(i)*energy(0)-gram(i,0)).expand()));
    let divergence=(Atom::num(3)-Atom::num(2)*Atom::var(epsilon))*energy(0);
    for factor in &measure.factors()[..5]{let action=q.iter().zip(&boost).fold(Atom::zero(),|sum,(coordinate,action)|{let AtomView::Var(v)=coordinate.as_view()else{panic!()};sum+factor.derivative(v.get_symbol())*action});assert!(action.expand().together().cancel().is_zero())}
    let mut temporal=measure.ibp("mixed-temporal-U/0",&indices,&u,&Atom::zero(),16).unwrap();
    let mut common_boost=measure.ibp("mixed-common-boost/0",&indices,&boost,&divergence,16).unwrap();
    assert_eq!(temporal.len(),4);assert_eq!(common_boost.len(),4);
    for row in temporal.iter_mut().chain(common_boost.iter_mut()){row.domain=row.domain.intersection(&admitted).unwrap()}
    let slots:Vec<usize>=match mode.as_str(){"original"=>vec![],"slot6"=>vec![6],"slot8"=>vec![8],"both-reverse"=>vec![8,6],_=>vec![6,8]};
    let mut added=Vec::new();
    for (rows,enabled) in [(&temporal,mode!="boost-only"),(&common_boost,mode!="temporal-only")]{if enabled{for &s in &slots{for row in rows{added.push(shifted(row.clone(),&indices,s))}}}}
    let audit=added.iter().map(|s|serde_json::json!({"id":s.id,"domain":bounds(&s.domain),"conditions":s.nonzero_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),"terms":s.terms.iter().map(|t|serde_json::json!({"shift":t.shift,"coefficient":t.coefficient.to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
    let mut rows=original;
    if !added.is_empty(){
        let ctx=GuardedContext::new_with_physical_arity(GuardedMeasureIdentity{measure:format!("mixed completion polynomial source presentation; original={}",record.measure),support:"old polynomial source guard pulled back separately along each ordinary completion".into(),orientation:"occupied future loop0; temporal U0 and simultaneous common boost".into(),normalization:"homogeneous source identities unchanged; no period values".into(),branch:"same sealed origin, same zeros; no inverse-energy or domain widening".into(),deformation:"original uncut0,1,2,4 D-eta".into()},roles,indices,parameters,added,11).unwrap();
        assert_eq!(ctx.sources().native_sources().coefficient_variables(),variables.as_ref().as_slice());
        rows.extend(ctx.sources().sources().iter().zip(ctx.sources().native_sources().rows()).map(|(m,r)|GuardedSource::new(m.id.clone(),r.clone(),m.domain.clone()).with_nonzero_conditions(m.nonzero_conditions.clone())));
    }
    let sources=Arc::new(GuardedSourceSystem::new(format!("{};exact-polynomial-mixed-completion-presentation-v1;mode={mode};occupiedE=F6-F8;no-domain-widening",record.measure),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
    let points:Vec<[i64;12]>=serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();assert_eq!(points.len(),15);assert!(points.iter().all(|p|admitted.contains(p)));
    let mut rules=Vec::new();let mut gaps=Vec::new();
    for p in &points{let found=sources.solve_domains(vec![IndexDomain::new(p.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},1).unwrap();gaps.push(found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>());rules.extend(found.rules)}
    let program=GuardedProgram::new(sources.clone(),rules,[]).unwrap();let encoded=program.encode_native(limits).unwrap();std::fs::write(out.join("program.bin"),&encoded).unwrap();let replay=GuardedProgram::decode_generated(&encoded,sources.clone(),limits).unwrap();
    let probes=points.iter().map(|p|{let x=program.apply(p).unwrap();let y=replay.apply(p).unwrap();assert_eq!(x.status,y.status);assert_eq!(x.terms,y.terms);assert_eq!(x.nonzero_conditions,y.nonzero_conditions);assert!(x.terms.keys().all(|i|admitted.contains(i)));serde_json::json!({"point":p,"status":format!("{:?}",x.status),"terms":x.terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":x.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})}).collect::<Vec<_>>();
    let output=serde_json::json!({"scope":"Validation-only first15 actual cut3 requests, native pointdepth3/cap1; original54 rows retained; added exact whole-row reindexing only; no closure/numerical claim","mode":mode,"source_count":sources.sources().len(),"rule_count":program.rules().len(),"input_identity":prepared.identity(),"factor_basis":measure.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"coordinates":q.iter().map(ToString::to_string).collect::<Vec<_>>(),"routing":family.inverse_routing().iter().map(|row|row.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),"occupied_energy":energy(0).to_string(),"energy_expansion":[[6,1],[8,-1]],"added_sources":audit,"same_zero_domains":true,"old_rules_imported":false,"encoded_decoded_replay":true,"all_rhs_labels_admitted":true,"source_boundary_positive_completion_terms_vanish":true,"probes":probes,"gaps":gaps,"proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>(),"shifts":s.seed.shifts})).collect::<Vec<_>>()).collect::<Vec<_>>()});
    std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&output).unwrap()).unwrap();println!("mode={mode} sources={} rules={} points={} replay=pass",sources.sources().len(),program.rules().len(),points.len());
}
