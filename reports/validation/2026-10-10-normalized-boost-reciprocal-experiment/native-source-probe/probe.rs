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
fn main(){
 let args:Vec<_>=std::env::args().collect();let mode=&args[5];let out=std::path::PathBuf::from(&args[4]);std::fs::create_dir_all(&out).unwrap();assert!(["original","polynomial-shifted","widened","widened-raw","normalized-first","normalized-last"].contains(&mode.as_str()));
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let widened=!["original","polynomial-shifted"].contains(&mode.as_str());
 let extend=|v:&Domain|{let mut d=v.clone();assert_eq!(d[6],(None,Some(0)));if widened{d[6]=(None,None)}domain::<12>(&d)};
 let original=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),extend(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect::<Vec<_>>();
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
 let mut bounds=*IndexDomain::for_roles(&roles).bounds();bounds[11]=IndexBounds::fixed(0);for b in &mut bounds[5..9]{*b=IndexBounds::new(None,Some(0)).unwrap()}if widened{bounds[6]=IndexBounds::unbounded()}let admitted=IndexDomain::new(bounds).unwrap();
 assert!(family.shells()[0].mass_squared.is_zero());assert_eq!(family.factors()[6],energy(occupied));assert_eq!(roles[6],IndexRole::Ordinary);assert_eq!(occupied,0);assert_eq!(family.shells()[0].chemical_potential,Rational::one());
 let mut identities=measure.ibp("global-boost/occupied-0",&indices,&direction,&divergence,16).unwrap();assert_eq!(identities.len(),4);
 for source in &mut identities{source.domain=source.domain.intersection(&admitted).unwrap();assert!(source.terms.iter().all(|t|t.shift[..5].iter().all(|&s|s<=0)),"boost raised a physical factor");assert!(source.terms.iter().all(|t|t.shift[11]==0));}
 let source_audit=identities.iter().map(|s|serde_json::json!({"id":s.id,"domain":bounds_of(&s.domain),"conditions":s.nonzero_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),"terms":s.terms.iter().map(|t|serde_json::json!({"shift":t.shift,"coefficient":t.coefficient.to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
 let ctx=GuardedContext::new_with_physical_arity(GuardedMeasureIdentity{measure:format!("validation-only correlated boost; original source context={}",record.measure),support:format!("original completion/support except exact future compact E0 all finite integers iff widened={widened}; no other domain extension"),orientation:"all loops simultaneous v_i=E_i P0-g_0i U, occupied future P0".into(),normalization:"homogeneous total derivative; original normalization unchanged".into(),branch:"finite-label joint high-D radial valuation; lower cutoff and independent line regulators removed at fixed positive eta,T; real-FermiT0 then meromorphicD; exact E0 reciprocal extension only when source epoch enables it".into(),deformation:"only original uncut physical factors1..4 D-eta, all scalar factors invariant".into()},roles,indices,parameters,identities,11).unwrap();
 assert_eq!(ctx.sources().native_sources().coefficient_variables(),variables.as_ref().as_slice());
 let boost=ctx.sources().sources().iter().zip(ctx.sources().native_sources().rows()).map(|(m,r)|GuardedSource::new(m.id.clone(),r.clone(),m.domain.clone()).with_nonzero_conditions(m.nonzero_conditions.clone())).collect::<Vec<_>>();
 let mut normal=Vec::new();
 for r in &boost {
   let mut shift=[0i64;12];shift[6]=1;
   let row=r.row.iter().map(|t|{let mut powers=*t.integral.powers();assert!(powers[6].is_symbolic());powers[6]=powers[6].shifted(1).unwrap();Term{integral:Integral::new(powers),coefficient:t.coefficient.shift_var(record.indices[6],&Integer::one())}}).collect();
   let pulled=r.domain.pullback(&shift).unwrap();assert_eq!(pulled.bounds()[6],if widened{IndexBounds::unbounded()}else{IndexBounds::new(None,Some(-1)).unwrap()});
   normal.push(GuardedSource::new(format!("normalized-global-boost/{}",r.id),row,pulled).with_nonzero_conditions(r.nonzero_conditions.iter().map(|c|c.shift_var(record.indices[6],&Integer::one())).collect()));
 }
 let normalized_audit=normal.iter().map(|r|serde_json::json!({"id":r.id,"domain":bounds_of(&r.domain),"conditions":r.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"terms":r.row.iter().map(|t|serde_json::json!({"powers":t.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"coefficient":t.coefficient.to_expression().to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
 let mut rows=Vec::new();if mode=="normalized-first"{rows.extend(normal.clone())}rows.extend(original);if mode=="widened-raw"{rows.extend(boost)}else if mode=="normalized-last"||mode=="polynomial-shifted"{rows.extend(normal)}
 let measure_id=format!("{};source-presentation=normalized-global-boost-attribution-v1;mode={mode};widened={widened};only-factor6=E0-future-massless-mu1;finite-label-joint-highD-radial-origin;fixed-positive-eta-T-before-real-Fermi-T0-and-meromorphicD;lower-contact-and-free-virtual-zero-boxes-widened-only-E0;all-original-targets-unchanged;shifted-entire-boost-row-a6+1-with-exact-pullback",record.measure);
 let sources=Arc::new(GuardedSourceSystem::new(measure_id,roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(extend).collect()).unwrap());
 let mut points:Vec<[i64;12]>=serde_json::from_slice::<Vec<Vec<i64>>>(&std::fs::read(&args[2]).unwrap()).unwrap().into_iter().map(|i|i.try_into().unwrap()).collect();let polynomial_count=points.len();for power in [1,2]{let mut point=[1,1,1,1,1,0,0,0,0,0,0,0];point[6]=power;points.push(point)}let mut rules=Vec::new();let mut gaps=Vec::new();
 for point in &points{let found=sources.solve_domains(vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},1).unwrap();gaps.push(found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds_of(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>());rules.extend(found.rules)}
 let program=GuardedProgram::new(sources.clone(),rules,[]).unwrap();let encoded=program.encode_native(limits).unwrap();std::fs::write(out.join("program.bin"),&encoded).unwrap();let replay=GuardedProgram::decode_generated(&encoded,sources.clone(),limits).unwrap();
 let probes=points.iter().map(|point|{let x=program.apply(point).unwrap();let y=replay.apply(point).unwrap();assert_eq!(x.status,y.status);assert_eq!(x.terms,y.terms);assert_eq!(x.nonzero_conditions,y.nonzero_conditions);serde_json::json!({"point":point,"status":format!("{:?}",x.status),"terms":x.terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":x.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})}).collect::<Vec<_>>();
 let result=serde_json::json!({"scope":"Matched normalized-global-boost attribution with separately derived finite-label reciprocal radial continuation; no production admission, numerical or closure claim","mode":mode,"source_count":sources.sources().len(),"rules":program.rules().len(),"input_identity":prepared.identity(),"factors":measure.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"coordinates":coordinates.iter().map(ToString::to_string).collect::<Vec<_>>(),"direction":direction.iter().map(ToString::to_string).collect::<Vec<_>>(),"divergence":divergence.to_string(),"all_scalar_and_physical_factor_actions_exact_zero":true,"all_boost_physical_shifts_nonpositive":true,"boost_sources":source_audit,"normalized_boost_sources":normalized_audit,"reciprocal_domain_enabled":widened,"polynomial_points":polynomial_count,"encoded_decoded_replay":true,"old_rules_imported":false,"polynomial_completion_domains_unchanged":!widened,"zero_domains_unchanged":!widened,"only_widened_axis":if widened{Some(6)}else{None},"new_zero_boxes_added":false,"probes":probes,"gaps":gaps,"proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>(),"shifts":s.seed.shifts})).collect::<Vec<_>>()).collect::<Vec<_>>()});std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();println!("mode={mode} points={} native_rules={} replay=pass",points.len(),program.rules().len());
}
fn bounds_of<const N:usize>(d:&IndexDomain<N>)->Domain{bounds(d)}
