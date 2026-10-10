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
use std::collections::{BTreeMap,BTreeSet,VecDeque};

type Point=[i64;12];type Row=BTreeMap<Point,Atom>;
struct Reduced { leaves:Row,terms:Row,conditions:BTreeSet<String>,failures:Vec<serde_json::Value>,applications:usize }
fn add(row:&mut Row,i:Point,c:Atom){*row.entry(i).or_default()+=c;}
fn clean(row:&mut Row){row.retain(|_,c|{*c=c.together().cancel();!c.is_zero()});}
fn reduce(program:&GuardedProgram<12>,row:&Row,admitted:&IndexDomain<12>)->Reduced{
 let mut r=Reduced{leaves:Row::new(),terms:Row::new(),conditions:BTreeSet::new(),failures:Vec::new(),applications:0};
 for (i,c) in row {if c.is_zero(){continue;}assert!(admitted.contains(i));let x=program.reduce(*i,Default::default()).unwrap();r.applications+=x.rule_applications;
  r.conditions.extend(x.nonzero_conditions.iter().map(|c|c.to_expression().to_string()));
  for(j,w)in x.terms{assert!(admitted.contains(&j));add(&mut r.terms,j,c*w.to_expression());}
  for t in x.unresolved{assert!(admitted.contains(&t.integral));if t.reason==rustred::solver::guarded::GuardedApplicationFailure::NoApplicableRule{add(&mut r.leaves,t.integral,c*t.coefficient.to_expression());}else{r.failures.push(serde_json::json!({"indices":t.integral,"coefficient":(c*t.coefficient.to_expression()).to_string(),"reason":format!("{:?}",t.reason)}));}}
 }
 clean(&mut r.leaves);clean(&mut r.terms);r
}
fn derivative(i:Point,shifted:&[usize])->Row{let mut row=Row::new();for&slot in shifted{if i[slot]!=0{let mut j=i;j[slot]+=1;add(&mut row,j,Atom::num(i[slot]));}}clean(&mut row);row}
fn display(row:&Row)->Vec<serde_json::Value>{row.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_string()})).collect()}
fn setting(name:&str,default:usize)->usize{std::env::var(name).ok().map(|s|s.parse().unwrap()).unwrap_or(default)}
fn main(){
 let args:Vec<_>=std::env::args().collect();let out=std::path::PathBuf::from(&args[2]);std::fs::create_dir_all(&out).unwrap();
 let bytes=std::fs::read(&args[1]).unwrap();let limits=BinaryIoLimits::default();let envelope=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(record.schema,"rustred.guarded-source-program.v1");
 let table=DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;12]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());

 let program_bytes=std::fs::read(&args[3]).unwrap();let program=GuardedProgram::decode_generated(&program_bytes,sources.clone(),limits).unwrap();
 let report:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[4]).unwrap()).unwrap();
 let points=report["attempted_exact_points"].as_array().unwrap().iter().map(|p|serde_json::from_value::<Point>(p.clone()).unwrap()).collect::<Vec<_>>();
 let occupations=roles.iter().enumerate().filter_map(|(i,r)|(*r==IndexRole::Occupation).then_some(i)).collect::<Vec<_>>();
 let mut choices=Vec::new();
 for point in &points{
  let applied=program.apply(point).unwrap();let rustred::solver::guarded::GuardedApplicationStatus::Applied{rule}=applied.status else{continue};
  let n: i64=occupations.iter().map(|&i|point[i]).sum();
  for (term,coefficient) in &applied.terms{
   let m:i64=occupations.iter().map(|&i|term[i]).sum();if m<=n{continue;}
   let r=&program.rules()[rule];let order=r.ordering();
   let a=Integral::numeric(point.map(|x|i16::try_from(x).unwrap())).unwrap();let b=Integral::numeric(term.map(|x|i16::try_from(x).unwrap())).unwrap();
   assert_eq!(order.compare(&a,&b),std::cmp::Ordering::Less);
   choices.push(((r.candidate().sources.len(),applied.terms.len(),point.iter().map(|v|v.abs()).sum::<i64>(),rule),serde_json::json!({"requested_point":point,"rhs_child":term,"coefficient":coefficient.to_expression().to_string(),"native_rule":rule,"native_harder_first_compare":"Less","parent_occupation_degree":n,"child_occupation_degree":m,"occupation_first_grade_would_compare":"Greater","applied_conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"rule_domain":bounds(r.domain()),"rule_conditions":r.nonzero_conditions().iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"native_rhs":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"proof_sources":r.candidate().sources.iter().map(|seed|{let s=&record.sources[seed.basis_row];serde_json::json!({"source_ordinal":seed.basis_row,"source_id":s.id,"source_guard":s.domain,"source_conditions":s.conditions.iter().map(|&id|polynomial(id).to_expression().to_string()).collect::<Vec<_>>(),"seed_integral":seed.seed.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"seed_shifts":seed.seed.shifts,"original_source_row":s.terms.iter().map(|t|serde_json::json!({"integral":t.integral,"coefficient":polynomial(t.coefficient).to_expression().to_string()})).collect::<Vec<_>>()})}).collect::<Vec<_>>()})));
  }
 }
 choices.sort_by_key(|x|x.0);let count=choices.len();let witness=choices.first().map(|x|x.1.clone());
 let result=serde_json::json!({"scope":"Native replayed rule witness, not alternate-order discovery and not a closure-obstruction proof","source_program":args[1],"replayed_program":args[3],"requested_points":points.len(),"replayed_rules":program.rules().len(),"occupation_raising_edges":count,"minimal_score":"proof source count, RHS count, requested-index absolute sum, rule index","witness":witness,"limitation":"Guarded native roles cannot currently be combined with a general compiled ordering program. Occupation-first grade is an illustrative comparison only; no alternative rule or period is installed."});
 std::fs::write(out.join("ordering-witness.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();println!("rules={} points={} raising_edges={} witness={}",program.rules().len(),points.len(),count,result["witness"]["native_rule"]);
}
