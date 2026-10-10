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


use rustred::algebra::{Coefficient,CoefficientPolynomial};
use rustred::solver::guarded::{GuardedReduction,GuardedReductionLimits,GuardedUnresolvedTerm,GuardedApplicationFailure,GuardedApplicationStatus};
use std::collections::{BTreeMap,BTreeSet};
fn add<const N:usize>(m:&mut BTreeMap<[i64;N],Coefficient>,key:[i64;N],value:Coefficient){
 if value.is_zero(){return} match m.entry(key){std::collections::btree_map::Entry::Vacant(e)=>{e.insert(value);},std::collections::btree_map::Entry::Occupied(mut e)=>{let sum=e.get()+&value;if sum.is_zero(){e.remove();}else{*e.get_mut()=sum;}}}
}
fn numeric<const N:usize>(label:&[i64;N])->Option<Integral<N>>{let mut values=[0i16;N];for(i,&n)in label.iter().enumerate(){values[i]=i16::try_from(n).ok()?;}Integral::numeric(values).ok()}
fn scheduled<const N:usize>(p:&GuardedProgram<N>,target:[i64;N],limits:GuardedReductionLimits,native_order:bool)->(GuardedReduction<N>,usize,usize){
 let one:Coefficient=p.sources().native_sources().rows().iter().flatten().next().unwrap().coefficient.one().into();let mut result=GuardedReduction{terms:BTreeMap::new(),unresolved:Vec::new(),nonzero_conditions:Vec::new(),rule_applications:0};
 if limits.max_pending_integrals==0{result.unresolved.push(GuardedUnresolvedTerm{integral:target,coefficient:one,reason:GuardedApplicationFailure::WorkLimit});return(result,0,0)}
 let mut pending=BTreeMap::from([(target,one)]);let mut uncovered=BTreeMap::new();let mut expanded=BTreeSet::new();let mut repeat_expansions=0;let mut peak_pending=1;
 while !pending.is_empty(){
  let key=if native_order && !p.rules().is_empty(){let order=p.rules()[0].ordering();*pending.keys().min_by(|a,b|match(numeric(a),numeric(b)){(Some(a),Some(b))=>order.compare(&a,&b),(None,Some(_))=>std::cmp::Ordering::Less,(Some(_),None)=>std::cmp::Ordering::Greater,(None,None)=>a.cmp(b)}).unwrap()}else{*pending.first_key_value().unwrap().0};
  let coefficient=pending.remove(&key).unwrap();
  let zero=p.sources().roles().iter().zip(&key).any(|(role,n)|*role==IndexRole::RequiredCut&&*n<=0)||p.sources().zero_domains().iter().any(|d|d.contains(&key));
  if result.rule_applications>=limits.max_rule_applications&&!p.terminals().contains(&key)&&!zero{
   result.unresolved.push(GuardedUnresolvedTerm{integral:key,coefficient,reason:GuardedApplicationFailure::WorkLimit});result.unresolved.extend(pending.into_iter().map(|(integral,coefficient)|GuardedUnresolvedTerm{integral,coefficient,reason:GuardedApplicationFailure::WorkLimit}));break;
  }
  let applied=p.apply(&key).unwrap();match applied.status{
   GuardedApplicationStatus::Zero=>{},GuardedApplicationStatus::Terminal=>add(&mut result.terms,key,coefficient),
   GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)=>add(&mut uncovered,key,coefficient),
   GuardedApplicationStatus::Unresolved(reason)=>result.unresolved.push(GuardedUnresolvedTerm{integral:key,coefficient,reason}),
   GuardedApplicationStatus::Applied{..}=>{result.rule_applications+=1;if !expanded.insert(key){repeat_expansions+=1}
    for condition in applied.nonzero_conditions{if !condition.is_constant()&&!result.nonzero_conditions.contains(&condition){result.nonzero_conditions.push(condition)}}
    for(child,value)in applied.terms{let weighted=&coefficient*&value;if weighted.is_zero(){continue}if pending.contains_key(&child)||pending.len()<limits.max_pending_integrals{add(&mut pending,child,weighted);}else{result.unresolved.push(GuardedUnresolvedTerm{integral:child,coefficient:weighted,reason:GuardedApplicationFailure::WorkLimit});}}
   }
  }peak_pending=peak_pending.max(pending.len());
 }
 result.unresolved.extend(uncovered.into_iter().map(|(integral,coefficient)|GuardedUnresolvedTerm{integral,coefficient,reason:GuardedApplicationFailure::NoApplicableRule}));(result,repeat_expansions,peak_pending)
}
fn value<const N:usize>(r:&GuardedReduction<N>)->serde_json::Value{serde_json::json!({"terms":r.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"unresolved":r.unresolved.iter().map(|t|serde_json::json!({"indices":t.integral.to_vec(),"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>(),"conditions":r.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"rule_applications":r.rule_applications})}
fn main(){
 let a=std::env::args().collect::<Vec<_>>();let output=std::path::Path::new(&a[3]);std::fs::create_dir_all(output).unwrap();
 let bytes=std::fs::read(&a[1]).unwrap();let limits=BinaryIoLimits::default();let env=inspect_program(&bytes,limits).unwrap();
 let(record,used):(Rec,usize)=bincode::decode_from_slice(env.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,env.section(SectionTag::PROGRAM).unwrap().len());
 let table=DecodedCoefficientTable::import_generated_normalized(env.section(SectionTag::SYMBOLICA_STATE).unwrap(),env.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap();
 let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
 let roles:[IndexRole;16]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
 let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
 let sources=Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
 let program=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();let targets:Vec<[i64;16]>=serde_json::from_slice(&std::fs::read(&a[2]).unwrap()).unwrap();
 let audit_lex=std::env::var("PROBE_COMPARE_LEX").unwrap_or_else(|_|"true".into())=="true";let mut results=Vec::new();let mut unique=BTreeSet::new();let mut count=0;
 for target in targets{
  let start=std::time::Instant::now();let(native,repeats,peak)=scheduled(&program,target,Default::default(),true);let native_seconds=start.elapsed().as_secs_f64();for term in &native.unresolved{if term.reason==GuardedApplicationFailure::NoApplicableRule{unique.insert(term.integral);}}
  let mut row=serde_json::json!({"target":target.to_vec(),"native_order":value(&native),"native_seconds":native_seconds,"native_repeat_expansions":repeats,"native_peak_pending":peak});
  if audit_lex{let start=std::time::Instant::now();let(lex,repeats,peak)=scheduled(&program,target,Default::default(),false);let lex_seconds=start.elapsed().as_secs_f64();let actual=program.reduce(target,Default::default()).unwrap();assert_eq!(value(&lex),value(&actual),"diagnostic lex loop differs from native API");
   let no_fail=lex.unresolved.iter().all(|r|r.reason==GuardedApplicationFailure::NoApplicableRule)&&native.unresolved.iter().all(|r|r.reason==GuardedApplicationFailure::NoApplicableRule);
   let condition_equal=lex.nonzero_conditions.iter().all(|c|native.nonzero_conditions.contains(c))&&native.nonzero_conditions.iter().all(|c|lex.nonzero_conditions.contains(c));
   let lex_terms=lex.unresolved.iter().filter(|r|r.reason==GuardedApplicationFailure::NoApplicableRule).map(|r|(r.integral,r.coefficient.clone())).collect::<BTreeMap<_,_>>();let native_terms=native.unresolved.iter().filter(|r|r.reason==GuardedApplicationFailure::NoApplicableRule).map(|r|(r.integral,r.coefficient.clone())).collect::<BTreeMap<_,_>>();
   row["lex"]=value(&lex);row["lex_seconds"]=serde_json::json!(lex_seconds);row["lex_repeat_expansions"]=serde_json::json!(repeats);row["lex_peak_pending"]=serde_json::json!(peak);row["condition_sets_equal"]=serde_json::json!(condition_equal);row["complete_residual_coefficients_equal"]=serde_json::json!(if no_fail{Some(lex_terms==native_terms&&lex.terms==native.terms)}else{None});
  }
  count+=1;results.push(row);if count%32==0{eprintln!("audited {count} targets, {} native-order residual labels",unique.len());std::fs::write(output.join("progress.json"),serde_json::to_vec_pretty(&serde_json::json!({"audited":count,"unique_no_rule":unique.len()})).unwrap()).unwrap();}
 }
 let report=serde_json::json!({"scope":"native apply only; diagnostic mirrors current lifecycle bookkeeping and budgets, changing only pending-key selection to common native numeric order. No discovery, rule changes or elimination","program_rules":program.rules().len(),"limits":{"max_rule_applications":100000,"max_pending_integrals":100000},"targets":results,"unique_no_rule_labels":unique.len(),"remaining_labels":unique.into_iter().map(|i|i.to_vec()).collect::<Vec<_>>()});std::fs::write(output.join("result.json"),serde_json::to_vec_pretty(&report).unwrap()).unwrap();eprintln!("native-order diagnostic completed {count} targets");
}
