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



use symbolica_amflow::finite_density::{DensityInput,guarded::{GuardedContext,GuardedMeasureIdentity,GuardedIdentity,GuardedIdentityTerm}};
fn main(){
 let args:Vec<_>=std::env::args().collect();let mode=&args[5];let out=std::path::PathBuf::from(&args[4]);std::fs::create_dir_all(&out).unwrap();assert!(["original","normal-ward","ward-normal"].contains(&mode.as_str()));
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let original=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect::<Vec<_>>();
 assert_eq!(original.len(),52);let variables=original[0].row[0].coefficient.variables().clone();
 let symbols=variables.iter().map(|v|{let PolyVariable::Symbol(s)=v else{panic!("non-symbol coefficient coordinate")};*s}).collect::<Vec<_>>();
 let indices:[Symbol;12]=record.indices.iter().map(|&i|symbols[i]).collect::<Vec<_>>().try_into().unwrap();let parameters=symbols.iter().enumerate().filter(|(i,_)|!record.indices.contains(i)).map(|(_,s)|*s).collect::<Vec<_>>();
 let parameter=|name:&str|*parameters.iter().find(|&&s|Atom::var(s).to_string()==name).expect("original parameter name");let eta=parameter("eta");let epsilon=parameter("epsilon");
 let input:DensityInput=serde_json::from_slice(&std::fs::read(&args[3]).unwrap()).unwrap();let prepared=input.prepare().unwrap();let family=prepared.occupied_cut(&[0],4096).unwrap().at_physical_masses();assert!(record.measure.contains(prepared.identity()));assert_eq!(family.loops(),3);assert_eq!(family.shells().len(),1);
 let old_factors=record.measure.split("; factors=[").nth(1).unwrap().split("]; massless").next().unwrap().replace("rustflow_occupied::","");assert_eq!(old_factors,family.factors().iter().map(ToString::to_string).collect::<Vec<_>>().join(", "));
 let measure=family.deformed_measure::<12>(eta,&[1,2,3,4]).unwrap();assert_eq!(measure.roles(),&roles);let coordinates=family.coordinates();let loops=family.loops();let pairs=(0..loops).flat_map(|i|(i..loops).map(move|j|(i,j))).collect::<Vec<_>>();let gram=|i:usize,j:usize|coordinates[pairs.iter().position(|&(a,b)|(a,b)==(i.min(j),i.max(j))).unwrap()].clone();let energy=|i:usize|coordinates[pairs.len()+i].clone();let occupied=family.shells()[0].loop_index;
 let mut direction=vec![Atom::zero();pairs.len()];direction.extend((0..loops).map(|i|(energy(i)*energy(occupied)-gram(i,occupied)).expand()));let divergence=(Atom::num(3)-Atom::num(2)*Atom::var(epsilon))*energy(occupied);
 let action=|f:&Atom|coordinates.iter().zip(&direction).fold(Atom::zero(),|a,(q,v)|{let AtomView::Var(q)=q.as_view() else{panic!()};a+f.derivative(q.get_symbol())*v}).expand().together().cancel();
 for f in &measure.factors()[..family.physical_slots()]{assert!(action(f).is_zero(),"physical factor not invariant")}
 for g in &coordinates[..pairs.len()]{assert!(action(g).is_zero())}
 let mut bounds=*IndexDomain::for_roles(&roles).bounds();bounds[11]=IndexBounds::fixed(0);for b in &mut bounds[5..9]{*b=IndexBounds::new(None,Some(0)).unwrap()}let admitted=IndexDomain::new(bounds).unwrap();
 // Independently derived polynomial-only singleton shell Ward identity.
 // No reciprocal-energy source seed or widened factor domain is used.
 assert_eq!(family.shells()[0].physical_slot,0);assert!(family.shells()[0].mass_squared.is_zero());
 assert_eq!(family.shells()[0].upper_slot,9);assert_eq!(family.shells()[0].lower_slot,10);
 for i in 0..loops {assert_eq!(family.factors()[6+i],energy(i));assert_eq!(roles[6+i],IndexRole::Ordinary);}
 let mu=Atom::num(family.shells()[0].chemical_potential.clone());let dminus2=Atom::num(2)-Atom::num(2)*Atom::var(epsilon);let upper=Atom::var(indices[9]);
 let recenter=true;let mut identities=Vec::new();
 for positive in [false,true] {
   let mut domain_bounds=*admitted.bounds();domain_bounds[0]=IndexBounds::fixed(1);for b in &mut domain_bounds[6..9]{*b=IndexBounds::fixed(0)}domain_bounds[10]=IndexBounds::fixed(0);
   let mut previous=[0i16;12];let mut next=[0i16;12];
   let (left,right)=if recenter {
     previous[9]=-1;domain_bounds[9]=if positive{IndexBounds::new(Some(2),None).unwrap()}else{IndexBounds::fixed(1)};
     if positive{(dminus2.clone()+Atom::num(1)-upper.clone(),(upper.clone()-Atom::num(1))*mu.clone())}else{(dminus2.clone(),-mu.clone())}
   } else {
     next[9]=1;domain_bounds[9]=if positive{IndexBounds::new(Some(1),None).unwrap()}else{IndexBounds::fixed(0)};
     if positive{(dminus2.clone()-upper.clone(),upper.clone()*mu.clone())}else{(dminus2.clone(),-mu.clone())}
   };
   identities.push(GuardedIdentity{id:format!("polynomial-singleton-normalized-boost/{}{}",if positive{"positive-upper"}else{"bulk"},if recenter{"/recentered"}else{""}),terms:vec![GuardedIdentityTerm{shift:previous,coefficient:left},GuardedIdentityTerm{shift:next,coefficient:right}],domain:IndexDomain::new(domain_bounds).unwrap(),nonzero_conditions:vec![]});
 }
 assert_eq!(identities.len(),2);
 let source_audit=identities.iter().map(|s|serde_json::json!({"id":s.id,"domain":bounds_of(&s.domain),"conditions":s.nonzero_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),"terms":s.terms.iter().map(|t|serde_json::json!({"shift":t.shift,"coefficient":t.coefficient.to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
 let ctx=GuardedContext::new_with_physical_arity(GuardedMeasureIdentity{measure:format!("validation-only independently derived polynomial singleton normalized boost; original source context={}",record.measure),support:"C1 massless singleton; all energy completions0; lower occupation0; arbitrary scalar Gram polynomial and original virtual powers; all original support unchanged".into(),orientation:"cutoff joint Lorentz boost on future null shell; independent tangent shell derivation, not inverse-energy seed replay".into(),normalization:"homogeneous total derivative; original normalization unchanged".into(),branch:"polynomial-singleton-normalized-boost-v1; original fixed eta joint high-D massless origin cutoff removed; no inverse-energy integrals".into(),deformation:"only original uncut physical factors1..4 D-eta, all scalar factors invariant".into()},roles,indices,parameters,identities,11).unwrap();
 assert_eq!(ctx.sources().native_sources().coefficient_variables(),variables.as_ref().as_slice());
 let boost=ctx.sources().sources().iter().zip(ctx.sources().native_sources().rows()).map(|(m,r)|GuardedSource::new(m.id.clone(),r.clone(),m.domain.clone()).with_nonzero_conditions(m.nonzero_conditions.clone())).collect::<Vec<_>>();
 let strict_sources=Arc::new(GuardedSourceSystem::new(format!("polynomial-singleton-normalized-boost-source-only-v1/{mode}"),roles,record.indices.clone().try_into().unwrap(),boost.clone()).unwrap());
 let base_point:[i64;12]=[1,1,1,1,1,0,0,0,0,1,0,0];let mut validation_points=vec![base_point];for h in [2,3]{let mut p=base_point;p[9]=h;validation_points.push(p)}
 for (axis,power) in [(0,2),(6,-1),(7,-1),(8,-1),(10,1)]{let mut p=base_point;p[axis]=power;validation_points.push(p)}
 let mut guard_validation=Vec::new();
 for (number,point) in validation_points.iter().enumerate(){let found=strict_sources.solve_domains(vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},1).unwrap();let p=GuardedProgram::new(strict_sources.clone(),found.rules,[]).unwrap();let result=p.apply(point).unwrap();assert_eq!(!result.terms.is_empty(),number<3,"strict source guard/application case {number}");let encoded=p.encode_native(limits).unwrap();let decoded=GuardedProgram::decode_generated(&encoded,strict_sources.clone(),limits).unwrap();assert_eq!(decoded.apply(point).unwrap().terms,result.terms);guard_validation.push(serde_json::json!({"point":point,"expect_recurrence":number<3,"status":format!("{:?}",result.status),"terms":result.terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"decoded_replay":true}));}
 let mut normal=Vec::new();
 for r in original.iter().filter(|r|r.id.starts_with("lorentz/0/3/")){
   let mut shift=[0i64;12];shift[6]=1;
   let row=r.row.iter().map(|t|{let mut powers=*t.integral.powers();assert!(powers[6].is_symbolic());powers[6]=powers[6].shifted(1).unwrap();Term{integral:Integral::new(powers),coefficient:t.coefficient.shift_var(record.indices[6],&Integer::one())}}).collect();
   let domain=r.domain.pullback(&shift).unwrap();assert_eq!(domain.bounds()[6],IndexBounds::new(None,Some(-1)).unwrap());
   normal.push(GuardedSource::new(format!("shifted-polynomial-U/{}",r.id),row,domain).with_nonzero_conditions(r.nonzero_conditions.iter().map(|c|c.shift_var(record.indices[6],&Integer::one())).collect()));
 }
 assert_eq!(normal.len(),4);let normal_audit=normal.iter().map(|s|serde_json::json!({"id":s.id,"domain":bounds_of(&s.domain),"terms":s.row.iter().map(|t|serde_json::json!({"powers":t.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"coefficient":t.coefficient.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":s.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})).collect::<Vec<_>>();
 let mut rows=original;
 if mode=="normal-ward"{rows.extend(normal);rows.extend(boost)}else if mode=="ward-normal"{rows.extend(boost);rows.extend(normal)}
 let measure_id=format!("{};source-presentation=combined-polynomial-U-shift-and-singleton-Ward-v1;mode={mode};U-original-seed-a6+1-pullback-a6<=-1;guard-C1-Eindices0-lower0;cutoff-shell-Euler-Dminus2;bulk-(D-2)I-muH1;positive-(D-2-s)I+s*muHnext;original-polynomial-domain-unchanged",record.measure);
 let sources=Arc::new(GuardedSourceSystem::new(measure_id,roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let points:Vec<[i64;12]>=serde_json::from_slice::<Vec<Vec<i64>>>(&std::fs::read(&args[2]).unwrap()).unwrap().into_iter().map(|i|i.try_into().unwrap()).collect();let mut rules=Vec::new();let mut gaps=Vec::new();
 for point in &points{let found=sources.solve_domains(vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},1).unwrap();gaps.push(found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds_of(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>());rules.extend(found.rules)}
 let program=GuardedProgram::new(sources.clone(),rules,[]).unwrap();let encoded=program.encode_native(limits).unwrap();std::fs::write(out.join("program.bin"),&encoded).unwrap();let replay=GuardedProgram::decode_generated(&encoded,sources.clone(),limits).unwrap();
 let probes=points.iter().map(|point|{let x=program.apply(point).unwrap();let y=replay.apply(point).unwrap();assert_eq!(x.status,y.status);assert_eq!(x.terms,y.terms);assert_eq!(x.nonzero_conditions,y.nonzero_conditions);serde_json::json!({"point":point,"status":format!("{:?}",x.status),"terms":x.terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":x.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})}).collect::<Vec<_>>();
 let result=serde_json::json!({"scope":"Matched combination of four original temporal-U shifted-polynomial rows and two independently derived strict C1/Eindices0/lower0 Ward rows; original domain, zeros and targets unchanged; native replay only, no inverse-energy admission or closure claim","mode":mode,"source_count":sources.sources().len(),"rules":program.rules().len(),"input_identity":prepared.identity(),"factors":measure.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"coordinates":coordinates.iter().map(ToString::to_string).collect::<Vec<_>>(),"direction":direction.iter().map(ToString::to_string).collect::<Vec<_>>(),"divergence":divergence.to_string(),"all_scalar_and_physical_factor_actions_exact_zero":true,"all_boost_physical_shifts_nonpositive":true,"additional_sources":source_audit,"shifted_polynomial_U_sources":normal_audit,"source_recentered":recenter,"encoded_decoded_replay":true,"old_rules_imported":false,"polynomial_completion_domains_unchanged":true,"zero_domains_unchanged":true,"strict_source_guard_validation":guard_validation,"probes":probes,"gaps":gaps,"proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>(),"shifts":s.seed.shifts})).collect::<Vec<_>>()).collect::<Vec<_>>()});std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();println!("mode={mode} points={} native_rules={} replay=pass",points.len(),program.rules().len());
}
fn bounds_of<const N:usize>(d:&IndexDomain<N>)->Domain{bounds(d)}
