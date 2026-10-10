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
use std::collections::BTreeMap;
use symbolica_amflow::finite_density::{DensityInput,guarded::GuardedIdentity};
fn atom_row(id:&GuardedIdentity<12>)->BTreeMap<[i16;12],Atom>{id.terms.iter().map(|t|(t.shift,t.coefficient.clone().together().cancel())).collect()}
fn equal_rows(a:&BTreeMap<[i16;12],Atom>,b:&BTreeMap<[i16;12],Atom>)->bool {
 let keys=a.keys().chain(b.keys()).collect::<std::collections::BTreeSet<_>>();
 keys.into_iter().all(|k|(a.get(k).cloned().unwrap_or_default()-b.get(k).cloned().unwrap_or_default()).together().cancel().is_zero())
}

fn main() {
 let args:Vec<_>=std::env::args().collect();
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let (record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r {0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let indices:[usize;12]=record.indices.clone().try_into().unwrap();
 let mode=&args[4];
 let make_original=||record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect::<Vec<_>>();
 let original=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,indices,make_original()).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let replay=GuardedProgram::decode_generated(&bytes,original.clone(),limits).unwrap();
 let sample=polynomial(record.sources[0].terms[0].coefficient);
 let variables=sample.variables().clone();
 let symbols=variables.iter().map(|v|{let PolyVariable::Symbol(s)=v else{panic!("non-symbol source variable")};*s}).collect::<Vec<_>>();
 let index_symbols:[Symbol;12]=indices.map(|i|symbols[i]);
 let eta=*symbols.iter().find(|&&s|Atom::var(s).to_string()=="eta").unwrap();
 let epsilon=*symbols.iter().find(|&&s|Atom::var(s).to_string()=="epsilon").unwrap();
 let input_path="examples/finite_density/massless_three_loop_chain.json";
 let input:DensityInput=serde_json::from_slice(&std::fs::read(input_path).unwrap()).unwrap();
 let input=input.prepare().unwrap();let family=input.occupied_cut(&[0],4096).unwrap().at_physical_masses();
 let shifted=(0..family.physical_slots()).filter(|&i|family.roles()[i]==IndexRole::Ordinary).collect::<Vec<_>>();
 let measure=family.deformed_measure::<12>(eta,&shifted).unwrap();assert_eq!(measure.roles(),&roles);assert_eq!(family.loops(),3);
 let mut admitted=*IndexDomain::for_roles(&roles).bounds();
 for b in &mut admitted[family.physical_slots()..family.input_slots()]{*b=IndexBounds::new(None,Some(0)).unwrap();}
 admitted[family.factors().len()..].fill(IndexBounds::fixed(0));let admitted=IndexDomain::new(admitted).unwrap();
 let standard=measure.lorentz_ibps(3,&(Atom::num(4)-Atom::num(2)*Atom::var(epsilon)),&index_symbols,4096).unwrap();
 let raw_atom=|s:&Sr|s.terms.iter().map(|t|(t.integral.iter().map(|&(symbolic,v)|{assert!(symbolic);v}).collect::<Vec<_>>().try_into().unwrap(),polynomial(t.coefficient).to_expression())).collect::<BTreeMap<[i16;12],Atom>>();
 // Match every preserved legacy row against independently regenerated geometry.
 // Legacy partitions all occupations; active partitions may have broader guards.
 for old in &record.sources {
   let name=old.id.split("/occupation-case-").next().unwrap();let old_domain=domain::<12>(&old.domain);let row=raw_atom(old);
   assert!(standard.iter().any(|new|new.id.split("/occupation-case-").next().unwrap()==name && new.domain.intersection(&old_domain).as_ref()==Some(&old_domain) && equal_rows(&atom_row(new),&row)),"geometry mismatch for {}",old.id);
 }
 let coords=family.coordinates();let pairs=(0..3).flat_map(|i|(i..3).map(move|j|(i,j))).collect::<Vec<_>>();let offset=pairs.len();
 let mut extra=Vec::new();let mut proof=Vec::new();
 for i in 0..3 {
  let energy=coords[offset+i].clone();
  let mut dir=pairs.iter().map(|&(a,b)| &energy*(Atom::num(i32::from(a==i))*&coords[offset+b]+Atom::num(i32::from(b==i))*&coords[offset+a])).collect::<Vec<_>>();
  dir.extend((0..3).map(|a|Atom::num(i32::from(a==i))*&energy));
  let generated=measure.ibp(&format!("temporal-euler/{i}"),&index_symbols,&dir,&Atom::num(1),4096).unwrap();
  let (slot,scale)=(family.physical_slots()..family.input_slots()).find_map(|slot|{
   let ratio=(&family.factors()[slot]/&energy).together().cancel();let scale=Rational::try_from(ratio.as_view()).ok()?;(!scale.is_zero()).then_some((slot,scale))
  }).expect("diagnostic fixture has no pure energy completion");
  for mut new in generated {
   new.domain=new.domain.intersection(&admitted).unwrap();let row=atom_row(&new);let mut audited=0;
   for old in record.sources.iter().filter(|s|s.id.starts_with(&format!("lorentz/{i}/3/"))) {
    let mut pulled=*domain::<12>(&old.domain).bounds();let b=pulled[slot];pulled[slot]=IndexBounds::new(b.lower().map(|v|v+1),b.upper().map(|v|v+1)).unwrap();
    let Some(overlap)=new.domain.intersection(&IndexDomain::new(pulled).unwrap())else{continue};
    let mut shifted_row=BTreeMap::new();
    for (mut shift,coeff) in raw_atom(old) {shift[slot]-=1;let coeff=coeff.replace(Atom::var(index_symbols[slot])).with(Atom::var(index_symbols[slot])-Atom::num(1))/Atom::num(scale.clone());*shifted_row.entry(shift).or_insert_with(Atom::new)+=coeff;}
    assert!(equal_rows(&row,&shifted_row),"Euler/shifted-U source mismatch for {} against {}",new.id,old.id);
    proof.push(serde_json::json!({"new_source":new.id,"original_source":old.id,"completion_slot":slot,"factor_over_energy":scale.to_string(),"pulled_guard_intersection":bounds(&overlap),"exact_coefficients_equal":true}));audited+=1;
   }
   assert!(audited>0);extra.push(new);
  }
 }
 let make_extra=||extra.iter().map(|id|GuardedSource::new(id.id.clone(),id.terms.iter().map(|t|{let c:rustred::algebra::Coefficient=t.coefficient.try_to_rational_polynomial(&Q,&Z,Some(variables.clone())).unwrap();assert!(c.denominator.is_one());Term{integral:Integral::symbolic(t.shift).unwrap(),coefficient:c.numerator}}).collect(),id.domain.clone())).collect::<Vec<_>>();
 let mut rows=make_original();match mode.as_str(){"original"=>{},"appended"=>rows.extend(make_extra()),"prepended"=>{let mut first=make_extra();first.extend(rows);rows=first;},_=>panic!("unknown mode")};
 let sources=Arc::new(GuardedSourceSystem::new(format!("{};generic-temporal-euler-presentation={mode};geometry={}",record.measure,input.identity()),roles,indices,rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let raw:Vec<Vec<i64>>=serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();
 let points:Vec<[i64;12]>=raw.into_iter().map(|p|p.try_into().unwrap()).collect();
 let mut probes=Vec::new();
 for target in points{
  let started=std::time::Instant::now();
  let found=sources.solve_domains(vec![IndexDomain::new(target.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},1).unwrap();
  let gaps=found.unresolved.iter().map(|g|serde_json::json!({"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
  let program=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();
  let applied=program.apply(&target).unwrap();
  probes.push(serde_json::json!({"target":target.to_vec(),"mode":mode,"seconds":started.elapsed().as_secs_f64(),"status":format!("{:?}",applied.status),"rhs":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"discovery_gaps":gaps,"rule_sources":program.rules().iter().flat_map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"ordinal":s.basis_row,"seed":(0..12).map(|i|s.seed.integral[i].value()).collect::<Vec<_>>()}))).collect::<Vec<_>>() }));
 }
 let sources_json=record.sources.iter().enumerate().map(|(i,s)|{serde_json::json!({"original_ordinal":i,"source_id":s.id,"domain":s.domain,"conditions":s.conditions.iter().map(|&id|polynomial(id).to_expression().to_string()).collect::<Vec<_>>(),"terms":s.terms.iter().map(|t|serde_json::json!({"label":t.integral,"coefficient":polynomial(t.coefficient).to_expression().to_string()})).collect::<Vec<_>>()})}).collect::<Vec<_>>();
 let output=serde_json::json!({"scope":"Bounded native fixed-point discovery and original-source replay; presentation comparison only, no analytic zero or amplitude value", "mode":mode,"source_program":args[1],"point_fixture":args[2],"points":probes.len(),"max_domains_per_point":1,"max_depth":3,"full_corpus_preserved":true,"geometry_original_row_checks":record.sources.len(),"shifted_source_span_checks":proof,"geometry":input.identity(),"additional_source_count":extra.len(),"historical_program_replayed_rules":replay.rules().len(),"original_source_count":record.sources.len(),"source_count":sources.sources().len(),"source_term_count":sources.native_sources().rows().iter().map(|r|r.len()).sum::<usize>(),"probes":probes,"original_source_rows":sources_json});
 std::fs::write(&args[3],serde_json::to_vec_pretty(&output).unwrap()).unwrap();
 println!("{}",serde_json::json!({"mode":mode,"points":output["points"],"source_count":output["source_count"],"applied":output["probes"].as_array().unwrap().iter().filter(|p|p["status"].as_str().unwrap().starts_with("Applied")).count(),"output":args[3]}));
}
