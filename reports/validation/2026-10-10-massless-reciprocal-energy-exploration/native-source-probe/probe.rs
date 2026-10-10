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


fn main(){
 let args:Vec<_>=std::env::args().collect();let mode=&args[5];let out=std::path::PathBuf::from(&args[4]);std::fs::create_dir_all(&out).unwrap();
 assert!(["polynomial","reciprocal","normal-first","normal-last","shifted-polynomial"].contains(&mode.as_str()));
 let widened=mode!="polynomial"&&mode!="shifted-polynomial";
 let geometry:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[3]).unwrap()).unwrap();
 assert_eq!(geometry["certified_completions"][0]["slot"],6);assert_eq!(geometry["certified_completions"][0]["coefficient"],"1");assert_eq!(geometry["certified_completions"][0]["mass_squared"],"0");
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 assert!(record.measure.contains(geometry["input_identity"].as_str().unwrap()));
 let old_factors=record.measure.split("; factors=[").nth(1).unwrap().split("]; massless").next().unwrap().replace("rustflow_occupied::","");
 let factors=geometry["factors"].as_array().unwrap().iter().map(|v|v.as_str().unwrap()).collect::<Vec<_>>().join(", ");assert_eq!(old_factors,factors);
 assert_eq!(record.sources.len(),52);assert_eq!(record.roles[6],0);assert!(record.sources.iter().all(|s|s.domain[6]==(None,Some(0))));assert!(record.zero_domains.iter().all(|s|s[6]==(None,Some(0))));
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let extend=|d:&Domain|{let mut v=d.clone();if widened{v[6]=(None,None)}assert_eq!(v[11],(Some(0),Some(0)));domain::<12>(&v)};
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),extend(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect::<Vec<_>>();
 let mut normal=Vec::new();if mode.starts_with("normal-")||mode=="shifted-polynomial"{
 for r in rows.iter().filter(|r|r.id.starts_with("lorentz/0/3/")){
  let mut shift=[0i64;12];shift[6]=1;
  let row=r.row.iter().map(|t|{let mut integral=*t.integral.powers();assert!(integral[6].is_symbolic());integral[6]=integral[6].shifted(1).unwrap();Term{integral:Integral::new(integral),coefficient:t.coefficient.shift_var(record.indices[6],&Integer::one())}}).collect();
  normal.push(GuardedSource::new(format!("reciprocal-normal/{}",r.id),row,r.domain.pullback(&shift).unwrap()).with_nonzero_conditions(r.nonzero_conditions.iter().map(|c|c.shift_var(record.indices[6],&Integer::one())).collect()));
 }assert_eq!(normal.len(),4);
 }
 let normal_count=normal.len();let mut all=Vec::new();if mode=="normal-first"{all.append(&mut normal)}all.extend(rows);all.append(&mut normal);
 let extended_identity=serde_json::json!({"epoch":"exploratory-massless-compact-energy-joint-high-D-v1","mode":mode,"original_source_context":record.measure,"geometry_descriptor":geometry,"new_admission":if !widened{"original polynomial completions"}else{"only slot6=E0 has all finite integer indices; all other completion/occupation/tail constraints unchanged"},"origin_protocol":"joint on-shell distribution; for each finite shifted row choose common sufficiently high ReD, remove lower cutoff and independent line regulators at fixed positive eta,T, then real-Fermi T0 and meromorphic D; no new uncut energy pole","zero_policy":if !widened{"original zero domains"}else{"same lower-contact and free-virtual zero boxes widened only in exact compact E0 under finite-label high-D radial valuation; virtual null direction unchanged"},"normal_policy":"twice the normal IBP: original temporal-U IBP with a6 replaced by a6+1 and every term shifted +1 in slot6; explicit upper/lower terms retained"}).to_string();
 let sources=Arc::new(GuardedSourceSystem::new(extended_identity,roles,record.indices.clone().try_into().unwrap(),all).unwrap().with_zero_domains(record.zero_domains.iter().map(extend).collect()).unwrap());
 let mut points:Vec<[i64;12]>=serde_json::from_slice::<Vec<Vec<i64>>>(&std::fs::read(&args[2]).unwrap()).unwrap().into_iter().map(|i|i.try_into().unwrap()).collect();
 for row in geometry["targets"].as_array().unwrap(){for t in row.as_array().unwrap(){let mut i=[0;12];for(j,v)in t["indices"].as_array().unwrap().iter().enumerate(){i[j]=v.as_i64().unwrap()}if !points.contains(&i){points.push(i)}}}
 let polynomial_count=points.len();let p=[1,1,1,1,1,0,0,0,0,0,0,0];for power in [1,2]{let mut i=p;i[6]=power;points.push(i)}
 let mut rules=Vec::new();let mut discoveries=Vec::new();let depth=std::env::var("PROBE_DEPTH").ok().map(|v|v.parse().unwrap()).unwrap_or(3);let budget=std::env::var("PROBE_DOMAINS").ok().map(|v|v.parse().unwrap()).unwrap_or(1);
 for (n,point)in points.iter().enumerate(){let start=std::time::Instant::now();let found=sources.solve_domains(vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],SearchOptions{max_depth:Some(depth),sample_seed:0,..Default::default()},budget).unwrap();
 discoveries.push(serde_json::json!({"point":point,"physical_polynomial_probe":n<polynomial_count,"seconds":start.elapsed().as_secs_f64(),"rules":found.rules.len(),"gaps":found.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>() }));rules.extend(found.rules);
 }
 let program=GuardedProgram::new(sources.clone(),rules,[]).unwrap();let encoded=program.encode_native(limits).unwrap();std::fs::write(out.join("program.bin"),&encoded).unwrap();let replay=GuardedProgram::decode_generated(&encoded,sources.clone(),limits).unwrap();
 let mut probes=Vec::new();for point in &points{let x=program.apply(point).unwrap();let y=replay.apply(point).unwrap();assert_eq!(x.status,y.status);assert_eq!(x.terms,y.terms);assert_eq!(x.nonzero_conditions,y.nonzero_conditions);
 probes.push(serde_json::json!({"point":point,"status":format!("{:?}",x.status),"terms":x.terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":x.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"rhs_inverse_energy_terms":x.terms.keys().filter(|i|i[6]>0).count()}));}
 let num=|i:[i64;12]|Integral::numeric(i.map(|x|x.try_into().unwrap())).unwrap();let order=rustred::solver::IntegralOrder::new(p.map(|v|v>0),roles.map(|r|r==IndexRole::RequiredCut)).with_roles(roles).unwrap();let mut inverse=p;inverse[6]=1;let mut pinched=inverse;pinched[1]=0;
 let result=serde_json::json!({"scope":"Fresh native exact replay on regenerated-energy-certified source extension; no historical rules imported, no numerical evaluation or production admission","mode":mode,"input_source_program":args[1],"geometry":args[3],"sources":sources.sources().len(),"normal_sources":normal_count,"reciprocal_domain_enabled":widened,"rules":program.rules().len(),"depth":depth,"domains_per_point":budget,"encoded_decoded_replay":true,"old_rules_imported":false,"measure_id":sources.measure_id(),"polynomial_points":polynomial_count,"discoveries":discoveries,"probes":probes,"rule_proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>(),"shifts":s.seed.shifts})).collect::<Vec<_>>()).collect::<Vec<_>>(),"order":{"convention":"Less is harder; valid RHS is easier","same_physical_support_inverse_vs_polynomial":format!("{:?}",order.compare(&num(inverse),&num(p))),"one_physical_pinch_inverse_vs_full_polynomial":format!("{:?}",order.compare(&num(pinched),&num(p)))},"source_guards":sources.sources().iter().map(|s|serde_json::json!({"id":s.id,"domain":bounds(&s.domain),"conditions":s.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})).collect::<Vec<_>>()});
 std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();println!("mode={mode} points={} rules={} replay=pass",points.len(),program.rules().len());
}
