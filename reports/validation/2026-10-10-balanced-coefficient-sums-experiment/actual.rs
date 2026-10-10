#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::Decode;
use rustred::persistence::{
    BinaryIoLimits, CoefficientId, DecodedCoefficientTable, SectionTag, inspect_program,
};
use rustred::solver::guarded::{
    GuardedProgram, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole,
};
use rustred::solver::{Integral, Power, SearchOptions, Term};
use std::sync::Arc;
type Domain = Vec<(Option<i64>, Option<i64>)>;
type Label = Vec<(bool, i16)>;
#[derive(Decode)]
struct Tr {
    integral: Label,
    coefficient: usize,
}
#[derive(Decode)]
struct Sr {
    id: String,
    domain: Domain,
    conditions: Vec<usize>,
    terms: Vec<Tr>,
}
#[derive(Decode)]
struct Sd {
    row: usize,
    integral: Label,
    shifts: Vec<i16>,
}
#[derive(Decode)]
struct Rr {
    fixed: Vec<Option<i16>>,
    target: Label,
    rhs: Vec<Tr>,
    sources: Vec<Sd>,
    domain: Domain,
    discovery_domain: Domain,
    conditions: Vec<usize>,
    sector: Vec<bool>,
    permutation: Option<Vec<usize>>,
}
#[derive(Decode)]
struct Rec {
    schema: String,
    measure: String,
    roles: Vec<u8>,
    indices: Vec<usize>,
    sources: Vec<Sr>,
    zero_domains: Vec<Domain>,
    rules: Vec<Rr>,
    terminals: Vec<Vec<i64>>,
}
fn domain<const N: usize>(v: &Domain) -> IndexDomain<N> {
    IndexDomain::new(
        v.iter()
            .map(|&(lo, hi)| IndexBounds::new(lo, hi).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
    .unwrap()
}
fn label<const N: usize>(v: &Label) -> Integral<N> {
    Integral::new(
        v.iter()
            .map(|&(s, p)| Power::new(s, p).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
}
fn bounds<const N: usize>(v: &IndexDomain<N>) -> Domain {
    v.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}

use std::collections::BTreeMap;
use symbolica::prelude::*;
fn experiment_stats<const N: usize>(
    rule: &rustred::solver::guarded::GuardedRule<N>,
) -> serde_json::Value {
    #[cfg(source_portfolio)]
    {
        rule.candidate().stats.guarded_portfolio.map(|s| serde_json::json!({
        "baseline_rhs_terms":s.baseline_rhs_terms,"baseline_rows":s.baseline_rows,
        "baseline_seeds":s.baseline_seeds,"baseline_independent_rows":s.baseline_independent_rows,
        "baseline_exact_trace_rows":s.baseline_exact_trace_rows,"baseline_direct_hit":s.baseline_direct_hit,
        "baseline_elapsed":s.baseline_elapsed.as_secs_f64(),"baseline_exact_materialization":s.baseline_exact_materialization.as_secs_f64(),
        "selected_rhs_terms":s.selected_rhs_terms,"selected_arm":s.selected_arm,
        "preseal_calls":s.preseal_calls,"preseal_elapsed":s.preseal_elapsed.as_secs_f64(),
        "policy_elapsed":s.policy_elapsed.as_secs_f64(),"total_elapsed":s.total_elapsed.as_secs_f64(),"completion":s.completion,
        "trials":s.trials.iter().map(|t|serde_json::json!({
            "selector":t.selector,"selected_source_rows":t.selected_source_rows,
            "attempted_rows":t.attempted_rows,"accepted_rows":t.accepted_rows,"seeds":t.seeds,
            "independent_rows":t.independent_rows,"exact_materialization":t.exact_materialization.as_secs_f64(),
            "guard_rejected_rows":t.guard_rejected_rows,"empty_rows":t.empty_rows,
            "exact_trace_rows":t.exact_trace_rows,"exact_trace_terms":t.exact_trace_terms,
            "rhs_terms":t.rhs_terms,"nonzero_conditions":t.nonzero_conditions,
            "completion":t.completion,"search_elapsed":t.search_elapsed.as_secs_f64(),
        })).collect::<Vec<_>>()
    })).unwrap_or(serde_json::Value::Null)
    }
    #[cfg(not(source_portfolio))]
    {
        let _ = rule;
        serde_json::Value::Null
    }
}



use std::cmp::Ordering;
use rustred::algebra::{Coefficient,CoefficientContext};
use rustred::solver::IntegralOrder;
use rustred::solver::guarded::{GuardedReduction,GuardedUnresolvedTerm,GuardedApplicationStatus,GuardedApplicationFailure,GuardedReductionLimits};
struct PendingKey<'a, const N: usize> {
    integral: [i64; N],
    numeric: Option<Integral<N>>,
    order: Option<&'a IntegralOrder<N>>,
}

impl<'a, const N: usize> PendingKey<'a, N> {
    fn new(integral: [i64; N], order: Option<&'a IntegralOrder<N>>) -> Self {
        let mut values = [0_i16; N];
        let numeric = integral
            .iter()
            .zip(&mut values)
            .try_for_each(|(&value, slot)| {
                *slot = i16::try_from(value).ok()?;
                Some(())
            })
            .and_then(|()| Integral::numeric(values).ok());
        Self {
            integral,
            numeric,
            order,
        }
    }
}

impl<const N: usize> PartialEq for PendingKey<'_, N> {
    fn eq(&self, other: &Self) -> bool {
        self.integral == other.integral
    }
}

impl<const N: usize> Eq for PendingKey<'_, N> {}

impl<const N: usize> PartialOrd for PendingKey<'_, N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for PendingKey<'_, N> {
    fn cmp(&self, other: &Self) -> Ordering {
        let ordering = match (self.order, &self.numeric, &other.numeric) {
            (Some(order), Some(left), Some(right)) => order.compare(left, right),
            (Some(_), None, Some(_)) => Ordering::Less,
            (Some(_), Some(_), None) => Ordering::Greater,
            _ => Ordering::Equal,
        };
        ordering.then_with(|| self.integral.cmp(&other.integral))
    }
}


#[derive(Default)]struct Stats{additions:usize,addition_seconds:f64,multiplies:usize,multiply_seconds:f64,apply_seconds:f64,buffer_flushes:usize,capacity_compactions:usize,popped_zero_keys:usize,peak_buffer:usize,peak_pending:usize}
impl Stats{fn json(&self)->serde_json::Value{serde_json::json!({"additions":self.additions,"addition_seconds":self.addition_seconds,"multiplies":self.multiplies,"multiply_seconds":self.multiply_seconds,"apply_seconds":self.apply_seconds,"buffer_flushes":self.buffer_flushes,"capacity_compactions":self.capacity_compactions,"popped_zero_keys":self.popped_zero_keys,"peak_buffer":self.peak_buffer,"peak_pending":self.peak_pending})}}
fn balanced(mut values:Vec<Coefficient>,stats:&mut Stats)->Option<Coefficient>{
 values.retain(|c|!c.is_zero());while values.len()>1{let mut next=Vec::with_capacity(values.len().div_ceil(2));let mut iter=values.into_iter();while let Some(a)=iter.next(){if let Some(b)=iter.next(){let start=std::time::Instant::now();let sum=&a+&b;stats.addition_seconds+=start.elapsed().as_secs_f64();stats.additions+=1;if !sum.is_zero(){next.push(sum);}}else{next.push(a);}}values=next;}values.pop()
}
fn push<K:Ord>(map:&mut BTreeMap<K,Vec<Coefficient>>,key:K,value:Coefficient,buffer:usize,stats:&mut Stats){
 if value.is_zero(){return;}match map.entry(key){std::collections::btree_map::Entry::Vacant(e)=>{e.insert(vec![value]);},std::collections::btree_map::Entry::Occupied(mut e)=>{e.get_mut().push(value);stats.peak_buffer=stats.peak_buffer.max(e.get().len());if e.get().len()>=buffer{stats.buffer_flushes+=1;let entries=std::mem::take(e.get_mut());if let Some(sum)=balanced(entries,stats){*e.get_mut()=vec![sum];}else{e.remove();}}}}
}
fn compact<K:Ord>(map:&mut BTreeMap<K,Vec<Coefficient>>,stats:&mut Stats){map.retain(|_,v|{let Some(c)=balanced(std::mem::take(v),stats)else{return false;};*v=vec![c];true});}
fn accumulate<K:Ord>(map:&mut BTreeMap<K,Coefficient>,key:K,value:Coefficient){if value.is_zero(){return;}match map.entry(key){std::collections::btree_map::Entry::Vacant(e)=>{e.insert(value);},std::collections::btree_map::Entry::Occupied(mut e)=>{let s=e.get()+&value;if s.is_zero(){e.remove();}else{*e.get_mut()=s;}}}}
fn reduce_balanced<const N:usize>(program:&GuardedProgram<N>,target:[i64;N],limits:GuardedReductionLimits,buffer:usize,progress:Option<&std::path::Path>)->(GuardedReduction<N>,Stats){
 assert!(buffer>=2);let mut stats=Stats::default();let start=std::time::Instant::now();
 let one:Coefficient=program.sources().native_sources().rows().iter().flatten().next().unwrap().coefficient.one().into();
 let mut result=GuardedReduction{terms:BTreeMap::new(),unresolved:Vec::new(),nonzero_conditions:Vec::new(),rule_applications:0};
 if limits.max_pending_integrals==0{result.unresolved.push(GuardedUnresolvedTerm{integral:target,coefficient:one,reason:GuardedApplicationFailure::WorkLimit});return(result,stats);}
 let order=program.rules().first().map(|r|r.ordering());let mut pending=BTreeMap::from([(PendingKey::new(target,order),vec![one])]);let mut uncovered=BTreeMap::new();let mut last_report=std::time::Instant::now();
 while let Some((key,values))=pending.pop_first(){
  let Some(coefficient)=balanced(values,&mut stats)else{stats.popped_zero_keys+=1;continue;};let integral=key.integral;
  if result.rule_applications>=limits.max_rule_applications&&!program.terminals().contains(&integral)&&!program.sources().is_zero(&integral){
   result.unresolved.push(GuardedUnresolvedTerm{integral,coefficient,reason:GuardedApplicationFailure::WorkLimit});
   for(key,values)in pending{if let Some(coefficient)=balanced(values,&mut stats){result.unresolved.push(GuardedUnresolvedTerm{integral:key.integral,coefficient,reason:GuardedApplicationFailure::WorkLimit});}}
   break;
  }
  let t=std::time::Instant::now();let applied=program.apply(&integral).unwrap();stats.apply_seconds+=t.elapsed().as_secs_f64();
  match applied.status{
   GuardedApplicationStatus::Zero=>{},GuardedApplicationStatus::Terminal=>accumulate(&mut result.terms,integral,coefficient),
   GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)=>accumulate(&mut uncovered,integral,coefficient),
   GuardedApplicationStatus::Unresolved(reason)=>result.unresolved.push(GuardedUnresolvedTerm{integral,coefficient,reason}),
   GuardedApplicationStatus::Applied{..}=>{result.rule_applications+=1;for condition in applied.nonzero_conditions{if !condition.is_constant()&&!result.nonzero_conditions.contains(&condition){result.nonzero_conditions.push(condition);}}
    for(child,value)in applied.terms{let t=std::time::Instant::now();let weighted=&coefficient*&value;stats.multiply_seconds+=t.elapsed().as_secs_f64();stats.multiplies+=1;if weighted.is_zero(){continue;}let key=PendingKey::new(child,order);
     if !pending.contains_key(&key)&&pending.len()>=limits.max_pending_integrals{stats.capacity_compactions+=1;compact(&mut pending,&mut stats);}
     if pending.contains_key(&key)||pending.len()<limits.max_pending_integrals{push(&mut pending,key,weighted,buffer,&mut stats);stats.peak_pending=stats.peak_pending.max(pending.len());}else{result.unresolved.push(GuardedUnresolvedTerm{integral:child,coefficient:weighted,reason:GuardedApplicationFailure::WorkLimit});}
    }
   }
  }
  if last_report.elapsed().as_secs()>=1{if let Some(p)=progress{let j=serde_json::json!({"applications":result.rule_applications,"pending":pending.len(),"elapsed_seconds":start.elapsed().as_secs_f64(),"stats":stats.json()});let temp=p.with_extension("json.part");std::fs::write(&temp,serde_json::to_vec_pretty(&j).unwrap()).unwrap();std::fs::rename(temp,p).unwrap();}last_report=std::time::Instant::now();}
 }
 result.unresolved.extend(uncovered.into_iter().map(|(integral,coefficient)|GuardedUnresolvedTerm{integral,coefficient,reason:GuardedApplicationFailure::NoApplicableRule}));(result,stats)
}
fn assert_equal<const N:usize>(a:&GuardedReduction<N>,b:&GuardedReduction<N>){assert_eq!(a.terms,b.terms);assert_eq!(a.rule_applications,b.rule_applications);assert_eq!(a.nonzero_conditions,b.nonzero_conditions);assert_eq!(a.unresolved.len(),b.unresolved.len());for(a,b)in a.unresolved.iter().zip(&b.unresolved){assert_eq!(a.integral,b.integral);assert_eq!(a.coefficient,b.coefficient);assert_eq!(a.reason,b.reason);}}

#[cfg(test)]mod tests{use super::*;use rustred::solver::{RuleCandidate,Case,SeedSource,Seed,SearchStats};
 fn discover(sources:Arc<GuardedSourceSystem<1>>,points:impl IntoIterator<Item=i64>,terminals:impl IntoIterator<Item=[i64;1]>)->GuardedProgram<1>{
  let mut rules=Vec::new();for point in points{let found=sources.solve_domains(vec![IndexDomain::new([IndexBounds::fixed(point)]).unwrap()],SearchOptions{max_depth:Some(0),sample_seed:0,..Default::default()},1).unwrap();assert!(!found.rules.is_empty(),"missing direct row at {point}: {:?}",found.unresolved);rules.extend(found.rules);}
  GuardedProgram::new(sources,rules,terminals).unwrap()
 }

    fn diamond_program(cancel: bool) -> GuardedProgram<1> {
        let context = CoefficientContext::try_new(["guarded_diamond_n", "guarded_diamond_x"]).unwrap();
        let x = context.parameter("guarded_diamond_x").unwrap();
        let target = Integral::symbolic([0]).unwrap();
        let branches = [
            (4, vec![(-1, context.one()), (-2, context.one())]),
            (2, vec![(-1, x.clone())]),
            (3, vec![(-2, if cancel { -x.clone() } else { x.clone() })]),
        ];
        let sources = Arc::new(
            GuardedSourceSystem::new(
                "replayed-uncovered-diamond",
                [IndexRole::Occupation],
                [0],
                branches
                    .iter()
                    .enumerate()
                    .map(|(ordinal, (point, rhs))| {
                        let mut row = vec![Term {
                            integral: target,
                            coefficient: context.one().numerator,
                        }];
                        row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                            integral: Integral::symbolic([*shift]).unwrap(),
                            coefficient: (-coefficient.clone()).numerator,
                        }));
                        GuardedSource::new(
                            format!("diamond-{ordinal}"),
                            row,
                            IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                        )
                        .with_nonzero_conditions(vec![x.numerator.clone()])
                    })
                    .collect(),
            )
            .unwrap(),
        );
        discover(sources, branches.iter().map(|(point, ..)| *point), [])
    }



    fn guarded_pending_diamond(cancel: bool) -> GuardedProgram<1> {
        let context = CoefficientContext::try_new(["pending_n", "pending_x", "pending_y"]).unwrap();
        let x = context.parameter("pending_x").unwrap();
        let y = context.parameter("pending_y").unwrap();
        let target = Integral::symbolic([0]).unwrap();
        let branches = [
            (
                6,
                context.one(),
                vec![(-1, context.one()), (-4, context.one())],
            ),
            (5, context.one(), vec![(-1, context.one())]),
            (
                4,
                context.one(),
                vec![(
                    -2,
                    if cancel {
                        -context.one()
                    } else {
                        context.one()
                    },
                )],
            ),
            (2, y.clone(), vec![(-1, &context.one() / &y)]),
        ];
        let sources = Arc::new(
            GuardedSourceSystem::new(
                "pending-native-order-guarded-diamond",
                [IndexRole::Occupation],
                [0],
                branches
                    .iter()
                    .enumerate()
                    .map(|(ordinal, (point, pivot, rhs))| {
                        let mut row = vec![Term {
                            integral: target,
                            coefficient: pivot.numerator.clone(),
                        }];
                        row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                            integral: Integral::symbolic([*shift]).unwrap(),
                            coefficient: (-(coefficient * pivot)).numerator,
                        }));
                        let source = GuardedSource::new(
                            format!("pending-{ordinal}"),
                            row,
                            IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                        );
                        if *point == 5 {
                            source.with_nonzero_conditions(vec![x.numerator.clone()])
                        } else {
                            source
                        }
                    })
                    .collect(),
            )
            .unwrap(),
        );
        discover(sources, branches.iter().map(|(point, ..)| *point), [[1]])
    }




 fn capacity_fixture()->GuardedProgram<1>{
  let context=CoefficientContext::try_new(["balanced_capacity_n","balanced_capacity_x"]).unwrap();let target=Integral::symbolic([0]).unwrap();let branches=[(8,vec![(-1,1),(-2,1)]),(7,vec![(-6,1)]),(6,vec![(-5,-1),(-4,1),(-3,1)])];
  let sources=Arc::new(GuardedSourceSystem::new("balanced-budget-cancellation",[IndexRole::Occupation],[0],branches.iter().enumerate().map(|(ordinal,(point,rhs))|{let mut row=vec![Term{integral:target,coefficient:context.one().numerator}];row.extend(rhs.iter().map(|(shift,c)|Term{integral:Integral::symbolic([*shift]).unwrap(),coefficient:context.integer(-i64::from(*c)).numerator}));GuardedSource::new(format!("budget-{ordinal}"),row,IndexDomain::new([IndexBounds::fixed(*point)]).unwrap())}).collect()).unwrap());
  discover(sources,branches.iter().map(|(point,_)|*point),[])
 }
 #[test]fn cancelled_pending_keys_do_not_create_false_capacity_failures(){let p=capacity_fixture();let limits=GuardedReductionLimits{max_rule_applications:100,max_pending_integrals:2};let native=p.reduce([8],limits).unwrap();let(lazy,s)=reduce_balanced(&p,[8],limits,32,None);assert_equal(&lazy,&native);assert!(lazy.unresolved.iter().all(|r|r.reason==GuardedApplicationFailure::NoApplicableRule));assert_eq!(lazy.unresolved.len(),2);assert_eq!(s.capacity_compactions,1);assert!(s.peak_pending<=2);}
 #[test]fn complete_native_budget_grid_preserves_all_conditions_failures_and_values(){let programs=[diamond_program(true),diamond_program(false),guarded_pending_diamond(true),guarded_pending_diamond(false),capacity_fixture()];for(p,start)in programs.iter().zip([4,4,6,6,8]){for applications in 0..=8{for pending in 0..=5{let limits=GuardedReductionLimits{max_rule_applications:applications,max_pending_integrals:pending};let native=p.reduce([start],limits).unwrap();for buffer in [2,4,32]{let(lazy,s)=reduce_balanced(p,[start],limits,buffer,None);assert_equal(&lazy,&native);assert!(s.peak_buffer<=buffer);assert!(s.peak_pending<=pending);}}}}}
 #[test]fn unused_guard_cancels_but_used_guard_survives_after_native_reload(){let p=guarded_pending_diamond(true);let p=GuardedProgram::decode_generated(&p.encode_native(Default::default()).unwrap(),p.sources().clone(),Default::default()).unwrap();let(a,s)=reduce_balanced(&p,[6],Default::default(),32,None);assert_equal(&a,&p.reduce([6],Default::default()).unwrap());let c=CoefficientContext::try_new(["pending_n","pending_x","pending_y"]).unwrap();assert!(a.nonzero_conditions.contains(&c.parameter("pending_x").unwrap().numerator));assert!(!a.nonzero_conditions.contains(&c.parameter("pending_y").unwrap().numerator));assert!(s.popped_zero_keys>0);}
 #[test]fn buffer_flush_is_bounded_and_balanced_rational_sum_is_exact(){let c=CoefficientContext::try_new(["balanced_sum_x"]).unwrap();let x=c.parameter("balanced_sum_x").unwrap();let mut map=BTreeMap::new();let mut stats=Stats::default();let mut expected=c.zero();for i in 1..=257{let term=&c.integer(if i%2==0{1}else{-1})/&(&x+&c.integer(i));expected=&expected+&term;push(&mut map,0,term,32,&mut stats);assert!(map[&0].len()<32);}let actual=balanced(map.remove(&0).unwrap(),&mut stats).unwrap();assert_eq!(actual,expected);assert!(stats.buffer_flushes>=8);assert_eq!(stats.peak_buffer,32);}

 #[test]fn thirty_one_thirty_two_thirty_three_and_cross_buffer_cancellation(){let c=CoefficientContext::try_new(["balanced_bound_x"]).unwrap();let x=c.parameter("balanced_bound_x").unwrap();for count in [31,32,33,65]{let mut map=BTreeMap::new();let mut stats=Stats::default();for _ in 0..count{push(&mut map,0,x.clone(),32,&mut stats);}for _ in 0..count{push(&mut map,0,-x.clone(),32,&mut stats);}compact(&mut map,&mut stats);assert!(map.is_empty());assert!(stats.peak_buffer<=32);}}
 #[test]fn coupled_index_guard_failure_remains_explicit(){
 let c=CoefficientContext::try_new(["balanced_guard_n","balanced_guard_m"]).unwrap();let n=c.parameter("balanced_guard_n").unwrap();let m=c.parameter("balanced_guard_m").unwrap();let d=IndexDomain::new([IndexBounds::new(Some(1),None).unwrap(),IndexBounds::new(Some(1),None).unwrap()]).unwrap();
 let sources=Arc::new(GuardedSourceSystem::new("balanced-coupled-guard",[IndexRole::Occupation;2],[0,1],vec![GuardedSource::new("coupled",vec![Term{integral:Integral::symbolic([0,0]).unwrap(),coefficient:(&n-&m).numerator},Term{integral:Integral::symbolic([-1,0]).unwrap(),coefficient:(-c.one()).numerator}],d.clone())]).unwrap());
 let found=sources.solve_domains(vec![d],SearchOptions{max_depth:Some(1),sample_seed:0,..Default::default()},1).unwrap();assert!(!found.rules.is_empty());let p=GuardedProgram::new(sources,found.rules,[]).unwrap();let a=p.apply(&[2,2]).unwrap();assert!(matches!(a.status,GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::ConditionVanished{..})));let expected=p.reduce([2,2],Default::default()).unwrap();let(got,_)=reduce_balanced(&p,[2,2],Default::default(),32,None);assert_equal(&got,&expected);assert_eq!(got.unresolved.len(),1);
 }
 #[test]fn invalid_and_unsupported_labels_are_explicit(){let p=diamond_program(true);for point in [-1,40000]{let limits=GuardedReductionLimits::default();assert_equal(&reduce_balanced(&p,[point],limits,32,None).0,&p.reduce([point],limits).unwrap());}}
}

fn process_cpu_ticks()->(u64,u64){
    let text=std::fs::read_to_string("/proc/self/stat").unwrap();let close=text.rfind(')').unwrap();let f=text[close+2..].split_whitespace().collect::<Vec<_>>();(f[11].parse().unwrap(),f[12].parse().unwrap())
}
fn progress(out:&std::path::Path, phase:&str, timings:&[serde_json::Value]){
    std::fs::write(out.join("progress.json"),serde_json::to_vec_pretty(&serde_json::json!({"phase":phase,"completed_timings":timings})).unwrap()).unwrap();println!("begin {phase}");
}
fn stamp(name:&str,start:std::time::Instant,cpu:(u64,u64))->serde_json::Value{
    let wall=start.elapsed().as_secs_f64();let now=process_cpu_ticks();println!("end {name} {wall:.6}s");serde_json::json!({"phase":name,"wall_seconds":wall,"user_cpu_ticks":now.0-cpu.0,"system_cpu_ticks":now.1-cpu.1})
}
fn main(){
    use std::time::Instant;
    use rustred::solver::guarded::GuardedReductionLimits;
    let args=std::env::args().collect::<Vec<_>>();let out=std::path::PathBuf::from(&args[3]);std::fs::create_dir_all(&out).unwrap();let mut timings=Vec::new();
    macro_rules! timed {($name:expr,$body:expr)=>{{progress(&out,$name,&timings);let cpu=process_cpu_ticks();let t=Instant::now();let value=$body;timings.push(stamp($name,t,cpu));value}}}
    let bytes=timed!("read_program_bytes",std::fs::read(&args[1]).unwrap());let limits=BinaryIoLimits::default();
    let envelope=timed!("inspect_envelope",inspect_program(&bytes,limits).unwrap());
    let record:Rec=timed!("decode_program_record",{let(r,used):(Rec,usize)=bincode::decode_from_slice(envelope.section(SectionTag::PROGRAM).unwrap(),bincode::config::standard()).unwrap();assert_eq!(used,envelope.section(SectionTag::PROGRAM).unwrap().len());assert_eq!(r.schema,"rustred.guarded-source-program.v2");r});
    let table=timed!("import_symbolica_and_coefficient_table",DecodedCoefficientTable::import_generated_normalized(envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),envelope.section(SectionTag::COEFFICIENTS).unwrap(),limits).unwrap());
    let polynomial=|id|{let c=table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone();assert!(c.denominator.is_one());c.numerator};
    let sources=timed!("construct_exact_source_context",{
        let roles:[IndexRole;16]=record.roles.iter().map(|r|match r{0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
        let rows=record.sources.iter().map(|s|GuardedSource::new(s.id.clone(),s.terms.iter().map(|t|Term{integral:label(&t.integral),coefficient:polynomial(t.coefficient)}).collect(),domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id|polynomial(id)).collect())).collect();
        Arc::new(GuardedSourceSystem::new(record.measure.clone(),roles,record.indices.clone().try_into().unwrap(),rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap())});
    let program=timed!("decode_generated_and_replay_all_rules",GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap());assert!(program.terminals().is_empty());assert_eq!(program.rules().len(),record.rules.len());
    let selections:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();
    let terminals:Vec<[i64;16]>=serde_json::from_value(selections["terminals"].clone()).unwrap();


    let count=program.rules().len();
    let program=timed!("bind_frontier_terminals",program.with_terminals_verified(terminals,count).unwrap());
    let target:[i64;16]=serde_json::from_value(selections["points"][25]["point"].clone()).unwrap();
    let start=Instant::now();let(result,stats)=reduce_balanced(&program,target,GuardedReductionLimits::default(),32,Some(&out.join("progress-balanced.json")));let elapsed=start.elapsed().as_secs_f64();
    let actual=serde_json::json!({"rule_applications":result.rule_applications,"terms":result.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"unresolved":result.unresolved.iter().map(|t|serde_json::json!({"indices":t.integral.to_vec(),"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>(),"conditions":result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()});
    // Read the prior native answer only after the complete independent reduction.
    let prior:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[4]).unwrap()).unwrap();
    let equal=actual==prior["points"][25]["reduction"];
    let final_report=serde_json::json!({"scope":"isolated balanced coefficient accumulation via unchanged native apply; no source search or period answer","target":target.to_vec(),"rules":count,"term_buffer":32,"native_limits":{"max_rule_applications":100000,"max_pending_integrals":100000},"elapsed_reduction_seconds":elapsed,"stats":stats.json(),"matches_prior_uninstrumented_native_exact_result":equal,"timings":timings,"reduction":actual});
    std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&final_report).unwrap()).unwrap();assert!(equal);println!("balanced target complete {elapsed:.6}s, {} applications, exact result match",result.rule_applications);
}

