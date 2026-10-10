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



#[derive(Clone)]struct Sample{ordinal:usize,values:Vec<Coefficient>,result:Option<Coefficient>,seconds:f64,group_sizes:Vec<usize>}
#[derive(Default)]struct Stats{additions:usize,addition_seconds:f64,multiplies:usize,multiply_seconds:f64,apply_seconds:f64,buffer_flushes:usize,capacity_compactions:usize,popped_zero_keys:usize,peak_buffer:usize,peak_pending:usize,buffers:usize,eligible_buffers:usize,eligible_seconds:f64,ineligible_seconds:f64,grouping_inspection_seconds:f64,buffer_histogram:BTreeMap<usize,usize>,eligible_histogram:BTreeMap<usize,usize>,samples:Vec<Sample>}
impl Stats{fn json(&self)->serde_json::Value{serde_json::json!({"additions":self.additions,"addition_seconds":self.addition_seconds,"multiplies":self.multiplies,"multiply_seconds":self.multiply_seconds,"apply_seconds":self.apply_seconds,"buffer_flushes":self.buffer_flushes,"capacity_compactions":self.capacity_compactions,"popped_zero_keys":self.popped_zero_keys,"peak_buffer":self.peak_buffer,"peak_pending":self.peak_pending,"buffers":self.buffers,"eligible_buffers":self.eligible_buffers,"eligible_seconds":self.eligible_seconds,"ineligible_seconds":self.ineligible_seconds,"grouping_inspection_seconds":self.grouping_inspection_seconds,"buffer_histogram":self.buffer_histogram,"eligible_histogram":self.eligible_histogram,"retained_samples":self.samples.len()})}}
fn denominator_groups(values:&[Coefficient])->Vec<Vec<usize>>{
 let mut groups:Vec<Vec<usize>>=Vec::new();for(i,c)in values.iter().enumerate(){assert_eq!(c.get_variables(),values[0].get_variables());if let Some(g)=groups.iter_mut().find(|g|values[g[0]].denominator==c.denominator){g.push(i);}else{groups.push(vec![i]);}}groups
}
fn balanced(mut values:Vec<Coefficient>,stats:&mut Stats)->Option<Coefficient>{
 values.retain(|c|!c.is_zero());stats.buffers+=1;*stats.buffer_histogram.entry(values.len()).or_default()+=1;
 let inspect=std::time::Instant::now();let groups=denominator_groups(&values);let eligible=groups.iter().any(|g|g.len()>1);stats.grouping_inspection_seconds+=inspect.elapsed().as_secs_f64();let saved=if eligible{Some(values.clone())}else{None};let start=std::time::Instant::now();
 while values.len()>1{let mut next=Vec::with_capacity(values.len().div_ceil(2));let mut iter=values.into_iter();while let Some(a)=iter.next(){if let Some(b)=iter.next(){let start=std::time::Instant::now();let sum=&a+&b;stats.addition_seconds+=start.elapsed().as_secs_f64();stats.additions+=1;if !sum.is_zero(){next.push(sum);}}else{next.push(a);}}values=next;}
 let result=values.pop();let seconds=start.elapsed().as_secs_f64();if let Some(original)=saved{stats.eligible_buffers+=1;stats.eligible_seconds+=seconds;*stats.eligible_histogram.entry(original.len()).or_default()+=1;let candidate=stats.samples.len()<8||seconds>stats.samples.last().unwrap().seconds;if candidate{stats.samples.push(Sample{ordinal:stats.buffers,values:original,result:result.clone(),seconds,group_sizes:groups.iter().map(Vec::len).collect()});stats.samples.sort_by(|a,b|b.seconds.total_cmp(&a.seconds).then(a.ordinal.cmp(&b.ordinal)));stats.samples.truncate(8);}}else{stats.ineligible_seconds+=seconds;}result
}
fn sum_balanced_plain(mut values:Vec<Coefficient>)->Option<Coefficient>{values.retain(|v|!v.is_zero());while values.len()>1{let mut next=Vec::with_capacity(values.len().div_ceil(2));let mut it=values.into_iter();while let Some(a)=it.next(){if let Some(b)=it.next(){let sum=&a+&b;if !sum.is_zero(){next.push(sum);}}else{next.push(a);}}values=next;}values.pop()}
fn sum_sequential(values:&[Coefficient])->Option<Coefficient>{let mut out:Option<Coefficient>=None;for value in values{if value.is_zero(){continue;}out=Some(match out{Some(c)=>&c+value,None=>value.clone()});}out.filter(|c|!c.is_zero())}
fn sum_grouped(values:&[Coefficient])->Option<Coefficient>{
 use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;
 let mut totals=Vec::new();for group in denominator_groups(values){if group.len()==1{totals.push(values[group[0]].clone());continue;}
 let first=&values[group[0]];let mut nums=group.iter().map(|&i|values[i].numerator.clone()).collect::<Vec<_>>();while nums.len()>1{let mut next=Vec::with_capacity(nums.len().div_ceil(2));let mut it=nums.into_iter();while let Some(a)=it.next(){next.push(match it.next(){Some(b)=>&a+&b,None=>a});}nums=next;}
 let c=Coefficient::from_num_den(nums.pop().unwrap(),first.denominator.clone(),first.numerator.ring(),true);if !c.is_zero(){totals.push(c);}}
 sum_balanced_plain(totals)
}
fn save_samples(out:&std::path::Path,stats:&Stats){
 let mut table=rustred::persistence::CoefficientTableBuilder::new(Default::default());let records=stats.samples.iter().map(|s|serde_json::json!({"ordinal":s.ordinal,"baseline_balanced_seconds":s.seconds,"group_sizes":s.group_sizes,"ids":s.values.iter().map(|c|table.intern(c).unwrap().index()).collect::<Vec<_>>(),"expected":s.result.as_ref().map(|c|table.intern(c).unwrap().index()),"input_shapes":s.values.iter().map(|c|serde_json::json!({"numerator_terms":c.numerator.nterms(),"denominator_terms":c.denominator.nterms(),"denominator_leading_coefficient":c.denominator.lcoeff().to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();let encoded=table.finish().unwrap();
 std::fs::create_dir_all(out).unwrap();std::fs::write(out.join("symbolica-state.bin"),&encoded.state).unwrap();std::fs::write(out.join("coefficients.bin"),&encoded.atoms).unwrap();let decoded=DecodedCoefficientTable::import_generated_normalized(&encoded.state,&encoded.atoms,Default::default()).unwrap();for(s,r)in stats.samples.iter().zip(&records){for(c,id)in s.values.iter().zip(r["ids"].as_array().unwrap()){assert_eq!(c,decoded.coefficient(CoefficientId::try_from_index(id.as_u64().unwrap()as usize).unwrap()).unwrap());}}
 std::fs::write(out.join("samples.json"),serde_json::to_vec_pretty(&serde_json::json!({"scope":"at most eight slowest observed balanced input buffers with exact repeated denominators; no answer-based selection","records":records,"stats":stats.json(),"native_table_roundtrip":true})).unwrap()).unwrap();
}
#[cfg(test)]mod grouping_tests{use super::*;use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;
 #[test]fn grouped_cancels_new_common_factors_and_signs(){let c=CoefficientContext::try_new(["group_x"]).unwrap();let x=c.parameter("group_x").unwrap();let a=&c.one()/&x;let b=&(&x-&c.one())/&x;assert_eq!(a.denominator,b.denominator);assert_eq!(sum_grouped(&[a.clone(),b.clone()]),Some(c.one()));assert_eq!(sum_grouped(&[a.clone(),-a.clone()]),None);let negative=Coefficient::from_num_den(-a.numerator.clone(),-a.denominator.clone(),a.numerator.ring(),true);assert_eq!(negative,a);assert_eq!(sum_grouped(&[negative,b]),Some(c.one()));}
 #[test]fn grouped_boundary_interleavings_match_canonical_add(){let c=CoefficientContext::try_new(["group_grid_x"]).unwrap();let x=c.parameter("group_grid_x").unwrap();for count in [0,1,2,31,32,33,65]{let v=(0..count).map(|i|&c.integer(if i%2==0{1}else{-1})/&(&x+&c.integer((i%3)+1))).collect::<Vec<_>>();assert_eq!(sum_grouped(&v),sum_sequential(&v));assert_eq!(sum_grouped(&v),sum_balanced_plain(v.clone()));for buffer in [2,4,32]{let groups=v.chunks(buffer).filter_map(sum_grouped).collect::<Vec<_>>();assert_eq!(sum_balanced_plain(groups),sum_sequential(&v));}}}
 #[test]fn grouped_polynomial_denominator_one_and_zero(){let c=CoefficientContext::try_new(["group_poly_x"]).unwrap();let x=c.parameter("group_poly_x").unwrap();let v=vec![x.clone(),-x,c.zero(),c.integer(2)];assert_eq!(sum_grouped(&v),Some(c.integer(2)));}
 #[test]fn incompatible_variable_maps_are_rejected(){let a=CoefficientContext::try_new(["group_map_x","group_map_y"]).unwrap();let b=CoefficientContext::try_new(["group_map_y","group_map_x"]).unwrap();let x=a.parameter("group_map_x").unwrap();let y=b.parameter("group_map_y").unwrap();assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(||sum_grouped(&[x,y]))).is_err());}
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


fn main(){
 let args=std::env::args().collect::<Vec<_>>();let inp=std::path::PathBuf::from(&args[1]);let out=std::path::PathBuf::from(&args[2]);let metadata:serde_json::Value=serde_json::from_slice(&std::fs::read(inp.join("samples.json")).unwrap()).unwrap();let table=DecodedCoefficientTable::import_generated_normalized(&std::fs::read(inp.join("symbolica-state.bin")).unwrap(),&std::fs::read(inp.join("coefficients.bin")).unwrap(),Default::default()).unwrap();let mut records=Vec::new();
 for(s,record)in metadata["records"].as_array().unwrap().iter().enumerate(){let values=record["ids"].as_array().unwrap().iter().map(|i|table.coefficient(CoefficientId::try_from_index(i.as_u64().unwrap()as usize).unwrap()).unwrap().clone()).collect::<Vec<_>>();assert!(values.len()<=32);let expected=record["expected"].as_u64().map(|i|table.coefficient(CoefficientId::try_from_index(i as usize).unwrap()).unwrap().clone());let mut runs=Vec::new();
 for repeat in 0..3{for method in [repeat%3,(repeat+1)%3,(repeat+2)%3]{let start=std::time::Instant::now();let value=match method{0=>sum_sequential(&values),1=>sum_balanced_plain(values.clone()),2=>sum_grouped(&values),_=>unreachable!()};let seconds=start.elapsed().as_secs_f64();assert_eq!(value,expected,"sample {s} method {method}");runs.push(serde_json::json!({"repeat":repeat,"method":(["sequential","balanced","equal-denominator-grouped"][method]),"seconds":seconds,"exact_match":true}));}}
 let mut medians=serde_json::Map::new();for name in ["sequential","balanced","equal-denominator-grouped"]{let mut times=runs.iter().filter(|r|r["method"]==name).map(|r|r["seconds"].as_f64().unwrap()).collect::<Vec<_>>();times.sort_by(f64::total_cmp);medians.insert(name.to_string(),serde_json::json!(times[1]));}records.push(serde_json::json!({"sample":s,"ordinal":record["ordinal"],"group_sizes":record["group_sizes"],"coefficient_count":values.len(),"runs":runs,"medians":medians}));std::fs::write(&out,serde_json::to_vec_pretty(&serde_json::json!({"scope":"matched exact sample summation; alternating method order, three repetitions; no native reduction or physical values","completed_samples":records.len(),"records":records})).unwrap()).unwrap();println!("sample {s} complete");}
}
