#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::{Decode,Encode};
use rustred::persistence::{inspect_program, BinaryIoLimits, SectionTag, DecodedCoefficientTable, CoefficientId};
use rustred::solver::guarded::{GuardedSource, GuardedSourceSystem, GuardedProgram, IndexRole, IndexBounds, IndexDomain};
use rustred::solver::{Integral, Power, Term, SearchOptions};
use std::sync::Arc;
type Domain=Vec<(Option<i64>,Option<i64>)>;
type Label=Vec<(bool,i16)>;
#[derive(Decode,Encode)] struct Tr { integral: Label, coefficient:usize }
#[derive(Decode,Encode)] struct Sr { id:String, domain:Domain, conditions:Vec<usize>, terms:Vec<Tr> }
#[derive(Decode,Encode)] struct Sd { row:usize, integral:Label, shifts:Vec<i16> }
#[derive(Decode,Encode)] struct Rr { fixed:Vec<Option<i16>>,target:Label,rhs:Vec<Tr>,sources:Vec<Sd>,domain:Domain,discovery_domain:Domain,conditions:Vec<usize>,sector:Vec<bool>,permutation:Option<Vec<usize>> }
#[derive(Decode,Encode)] struct Rec { schema:String,measure:String,roles:Vec<u8>,indices:Vec<usize>,sources:Vec<Sr>,zero_domains:Vec<Domain>,rules:Vec<Rr>,terminals:Vec<Vec<i64>> }
fn domain<const N:usize>(v:&Domain)->IndexDomain<N>{IndexDomain::new(v.iter().map(|&(lo,hi)|IndexBounds::new(lo,hi).unwrap()).collect::<Vec<_>>().try_into().unwrap()).unwrap()}
fn label<const N:usize>(v:&Label)->Integral<N>{Integral::new(v.iter().map(|&(s,p)|Power::new(s,p).unwrap()).collect::<Vec<_>>().try_into().unwrap())}
fn bounds<const N:usize>(v:&IndexDomain<N>)->Domain{v.bounds().iter().map(|b|(b.lower(),b.upper())).collect()}

fn main(){
 let a=std::env::args().collect::<Vec<_>>();let output=std::path::Path::new(&a[2]);std::fs::create_dir_all(output).unwrap();
 let bytes=std::fs::read(&a[1]).unwrap();let limits=BinaryIoLimits::default();let env=inspect_program(&bytes,limits).unwrap();
 let(mut record,used):(Rec,usize)=bincode::decode_from_slice(env.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,env.section(SectionTag::PROGRAM).unwrap().len());
 let table=DecodedCoefficientTable::import_generated_normalized(env.section(SectionTag::SYMBOLICA_STATE).unwrap(),env.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;16]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let original=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();assert!(original.terminals().is_empty());
 let source_dump=serde_json::json!({"scope":"original polynomial sources, not additional identities or generated rules","measure_id":sources.measure_id(),"index_variable_positions":sources.native_sources().index_variables().to_vec(),"coefficient_variables":sources.native_sources().coefficient_variables().iter().map(|x|format!("{x}")).collect::<Vec<_>>(),"roles":roles.iter().map(|r|format!("{r:?}")).collect::<Vec<_>>(),"zero_domains":sources.zero_domains().iter().map(bounds).collect::<Vec<_>>(),"sources":sources.sources().iter().zip(sources.native_sources().rows()).enumerate().filter(|(_, (info,_))|info.id.starts_with("lorentz/2/0/")||info.id.starts_with("lorentz/2/1/")).map(|(ordinal,(info,row))|serde_json::json!({"ordinal":ordinal,"id":info.id,"domain":bounds(&info.domain),"conditions":info.nonzero_conditions.iter().map(|x|x.to_expression().to_string()).collect::<Vec<_>>(),"terms":row.iter().map(|t|serde_json::json!({"integral":t.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"coefficient":t.coefficient.to_expression().to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>()});
 std::fs::write(output.join("source-rows.json"),serde_json::to_vec_pretty(&source_dump).unwrap()).unwrap();
 let targets:Vec<[i64;16]>=vec![
 [1,2,2,1,0,-1,0,0,0,0,0,0,0,0,0,0], [1,3,2,1,0,-1,0,0,0,0,0,0,0,0,0,0],
 [1,1,2,1,0,-1,0,0,0,0,0,0,0,0,0,0], [1,1,1,1,0,-1,0,0,0,0,0,0,0,0,0,0],
 [1,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,0], [1,2,0,1,0,-2,0,0,0,0,0,0,0,0,0,0],
 [1,1,0,1,0,-2,0,0,0,0,0,0,0,0,0,0],
 [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0], [2,1,1,1,1,-1,0,0,0,0,0,0,0,0,0,0], [2,1,1,1,1,0,-1,-1,0,0,0,0,0,0,0,0]];
 let audit=|program:&GuardedProgram<16>|->serde_json::Value{
   serde_json::Value::Array(targets.iter().map(|target|{
    let matching=program.rules().iter().enumerate().filter(|(_,rule)|rule.domain().contains(target)&&rule.candidate().target.powers().iter().zip(target).all(|(p,n)|p.is_symbolic()||i64::from(p.value())==*n)).map(|(i,rule)|serde_json::json!({"ordinal":i,"rhs_terms":rule.candidate().rhs.len(),"conditions":rule.nonzero_conditions().iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})).collect::<Vec<_>>();
    let first=program.apply(target).unwrap();let result=program.reduce(*target,Default::default()).unwrap();
    serde_json::json!({"target":target.to_vec(),"structural_matches":matching,"first_application":format!("{:?}",first.status),"first_terms":first.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"rule_applications":result.rule_applications,"unresolved":result.unresolved.iter().map(|r|serde_json::json!({"indices":r.integral.to_vec(),"coefficient":r.coefficient.to_expression().to_string(),"reason":format!("{:?}",r.reason)})).collect::<Vec<_>>(),"conditions":result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})
   }).collect())
 };
 let mut results=vec![serde_json::json!({"policy":"original","audit":audit(&original)})];
 for mode in ["zero-first","shortest-rhs"]{
  if mode=="zero-first"{record.rules.sort_by_key(|r|!r.rhs.is_empty());}else{record.rules.sort_by_key(|r|r.rhs.len());}
  let structure=bincode::encode_to_vec(&record,bincode::config::standard()).unwrap();
  let generated=rustred::persistence::encode_program(env.kind(),&[rustred::persistence::BinarySection{tag:SectionTag::SYMBOLICA_STATE,bytes:env.section(SectionTag::SYMBOLICA_STATE).unwrap()},rustred::persistence::BinarySection{tag:SectionTag::COEFFICIENTS,bytes:env.section(SectionTag::COEFFICIENTS).unwrap()},rustred::persistence::BinarySection{tag:SectionTag::PROGRAM,bytes:&structure}],limits).unwrap();
  let replayed=GuardedProgram::decode_generated(&generated,sources.clone(),limits).unwrap();std::fs::write(output.join(format!("{mode}.bin")),generated).unwrap();
  results.push(serde_json::json!({"policy":mode,"audit":audit(&replayed)}));
 }
 let report=serde_json::json!({"scope":"original source-bound program independently replayed; only stable rule precedence changed via generated transport and every reordered rule replayed again; no discovery or physics changes","source_count":sources.sources().len(),"rules":original.rules().len(),"zero_rhs_rules":original.rules().iter().filter(|r|r.candidate().rhs.is_empty()).count(),"policies":results});std::fs::write(output.join("result.json"),serde_json::to_vec_pretty(&report).unwrap()).unwrap();println!("rule-choice audit complete, {} rules",original.rules().len());
}
