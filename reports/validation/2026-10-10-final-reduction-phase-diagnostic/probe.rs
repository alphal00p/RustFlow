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
    let mut observations=Vec::new();
    for point in &terminals { let a=program.apply(point).unwrap();observations.push(serde_json::json!({"point":point.to_vec(),"status":format!("{:?}",a.status),"terms":a.terms.len()})); }
    std::fs::write(out.join("terminal-observations.json"),serde_json::to_vec_pretty(&observations).unwrap()).unwrap();
    assert!(observations.iter().all(|v|v["status"]=="Unresolved(NoApplicableRule)"));
    let program=timed!("bind_frontier_terminals",program.with_terminals_verified(terminals,count).unwrap());
    let selected:Vec<serde_json::Value>=serde_json::from_value(selections["points"].clone()).unwrap();
    for i in [24,25] {
        let p:[i64;16]=serde_json::from_value(selected[i]["point"].clone()).unwrap();
        let sub=out.join(format!("point-{i}"));std::fs::create_dir_all(&sub).unwrap();
        let result=traced_reduce(&program,p,&sub);
        if let Some(result)=result {
            if i==25 {
                let saved:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[4]).unwrap()).unwrap();
                let expected=&saved["points"][25]["reduction"];
                let actual=serde_json::json!({"rule_applications":result.rule_applications,"terms":result.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"unresolved":result.unresolved.iter().map(|t|serde_json::json!({"indices":t.integral.to_vec(),"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>(),"conditions":result.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()});
                assert_eq!(&actual,expected);
                std::fs::write(sub.join("native-equivalence.json"),serde_json::to_vec_pretty(&serde_json::json!({"matches_prior_uninstrumented_native_result_exactly":true,"scope":"same ordered indices, rendered exact rational coefficients, conditions, statuses and application count","reduction":actual})).unwrap()).unwrap();continue;
            }
            let reference=program.reduce(p,GuardedReductionLimits::default()).unwrap();
            assert_eq!(result.terms,reference.terms);assert_eq!(result.rule_applications,reference.rule_applications);assert_eq!(result.nonzero_conditions,reference.nonzero_conditions);
            assert_eq!(result.unresolved.len(),reference.unresolved.len());
            for (a,b)in result.unresolved.iter().zip(&reference.unresolved){assert_eq!(a.integral,b.integral);assert_eq!(a.coefficient,b.coefficient);assert_eq!(a.reason,b.reason);}
            std::fs::write(sub.join("native-equivalence.json"),serde_json::to_vec_pretty(&serde_json::json!({"exact_native_equivalence":true,"rule_applications":result.rule_applications,"terms":result.terms.len(),"unresolved":result.unresolved.len(),"conditions":result.nonzero_conditions.len()})).unwrap()).unwrap();
        } else {break;}
    }
}

use std::cmp::Ordering;
use rustred::solver::IntegralOrder;
use rustred::algebra::{Coefficient,CoefficientPolynomial};
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


struct Trace {directory:std::path::PathBuf,start:std::time::Instant,sequence:usize,totals:BTreeMap<String,(usize,f64)>,iterations:usize,applications:usize,pending:usize,conditions:usize}
impl Trace {
 fn phase<T>(&mut self,name:&str,label:[i64;16],sizes:serde_json::Value,f:impl FnOnce()->T)->T {
  self.sequence+=1;let report=serde_json::json!({"phase":name,"sequence":self.sequence,"elapsed_seconds":self.start.elapsed().as_secs_f64(),"label":label.to_vec(),"coefficient_sizes":sizes,"iterations":self.iterations,"applications":self.applications,"pending":self.pending,"conditions":self.conditions,"completed_phase_totals":self.totals});
  let temporary=self.directory.join("progress.json.part");std::fs::write(&temporary,serde_json::to_vec_pretty(&report).unwrap()).unwrap();std::fs::rename(temporary,self.directory.join("progress.json")).unwrap();
  let t=std::time::Instant::now();let result=f();let elapsed=t.elapsed().as_secs_f64();let entry=self.totals.entry(name.into()).or_default();entry.0+=1;entry.1+=elapsed;result
 }
}
fn size(c:&Coefficient)->serde_json::Value {serde_json::json!({"numerator_terms":c.numerator.nterms(),"denominator_terms":c.denominator.nterms()})}
fn add<K:Ord>(map:&mut BTreeMap<K,Coefficient>,key:K,value:Coefficient,trace:&mut Trace,phase:&str,label:[i64;16]){
 if value.is_zero(){return;}match map.entry(key){std::collections::btree_map::Entry::Vacant(e)=>{e.insert(value);},std::collections::btree_map::Entry::Occupied(mut e)=>{let sizes=serde_json::json!({"left":size(e.get()),"right":size(&value)});let sum=trace.phase(phase,label,sizes,||e.get()+&value);if sum.is_zero(){e.remove();}else{*e.get_mut()=sum;}}}
}
fn traced_reduce(program:&GuardedProgram<16>,target:[i64;16],directory:&std::path::Path)->Option<GuardedReduction<16>>{
 let mut trace=Trace{directory:directory.into(),start:std::time::Instant::now(),sequence:0,totals:BTreeMap::new(),iterations:0,applications:0,pending:0,conditions:0};
 let one:Coefficient=program.sources().native_sources().rows().iter().flatten().next().unwrap().coefficient.one().into();
 let mut result=GuardedReduction{terms:BTreeMap::new(),unresolved:Vec::new(),nonzero_conditions:Vec::new(),rule_applications:0};
 let order=program.rules().first().map(|r|r.ordering());let mut pending=BTreeMap::from([(PendingKey::new(target,order),one)]);let mut uncovered=BTreeMap::new();let limits=GuardedReductionLimits::default();
 while let Some((key,coefficient))=pending.pop_first(){
  let integral=key.integral;trace.iterations+=1;trace.pending=pending.len();trace.applications=result.rule_applications;trace.conditions=result.nonzero_conditions.len();
  if trace.start.elapsed().as_secs()>480 || result.rule_applications>=limits.max_rule_applications || pending.len()>limits.max_pending_integrals{std::fs::write(directory.join("bounded-stop.json"),serde_json::to_vec_pretty(&serde_json::json!({"pending":pending.len(),"rule_applications":result.rule_applications,"reason":"diagnostic wall/work bound between native operations; no completed reduction claim"})).unwrap()).unwrap();return None;}
  let applied=trace.phase("native-apply",integral,size(&coefficient),||program.apply(&integral).unwrap());
  match applied.status {
   GuardedApplicationStatus::Zero=>{},
   GuardedApplicationStatus::Terminal=>add(&mut result.terms,integral,coefficient,&mut trace,"terminal-coefficient-add",integral),
   GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)=>add(&mut uncovered,integral,coefficient,&mut trace,"uncovered-coefficient-add",integral),
   GuardedApplicationStatus::Unresolved(reason)=>result.unresolved.push(GuardedUnresolvedTerm{integral,coefficient,reason}),
   GuardedApplicationStatus::Applied{..}=>{result.rule_applications+=1;
    for condition in applied.nonzero_conditions{trace.phase("condition-retention",integral,serde_json::Value::Null,||{if !condition.is_constant()&&!result.nonzero_conditions.contains(&condition){result.nonzero_conditions.push(condition);}});}
    for(child,value)in applied.terms{
     let weighted=trace.phase("coefficient-multiply",child,serde_json::json!({"left":size(&coefficient),"right":size(&value)}),||&coefficient*&value);
     if weighted.is_zero(){continue;}let key=PendingKey::new(child,order);
     if pending.contains_key(&key)||pending.len()<limits.max_pending_integrals {add(&mut pending,key,weighted,&mut trace,"pending-coefficient-add",child);}else{result.unresolved.push(GuardedUnresolvedTerm{integral:child,coefficient:weighted,reason:GuardedApplicationFailure::WorkLimit});}
    }
   }
  }
 }
 result.unresolved.extend(uncovered.into_iter().map(|(integral,coefficient)|GuardedUnresolvedTerm{integral,coefficient,reason:GuardedApplicationFailure::NoApplicableRule}));
 std::fs::write(directory.join("complete.json"),serde_json::to_vec_pretty(&serde_json::json!({"iterations":trace.iterations,"rule_applications":result.rule_applications,"wall_seconds":trace.start.elapsed().as_secs_f64(),"phase_totals":trace.totals})).unwrap()).unwrap();Some(result)
}
