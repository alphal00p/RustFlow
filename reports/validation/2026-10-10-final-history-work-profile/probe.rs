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
    let program=timed!("bind_frontier_terminals_after_full_cold_replay",program.with_terminals_verified(terminals,count).unwrap());
    let work=GuardedReductionLimits::default();
    let selected:Vec<serde_json::Value>=serde_json::from_value(selections["points"].clone()).unwrap();assert!(selected.len()<=64);
    let mut points=Vec::new();
    for(i,selection)in selected.into_iter().enumerate(){
        let p:[i64;16]=serde_json::from_value(selection["point"].clone()).unwrap();
        let reduced=timed!(&format!("point-{i:02}-reduce-native"),program.reduce(p,work).unwrap());
        let red_json=timed!(&format!("point-{i:02}-reduce-expression-json"),serde_json::json!({"rule_applications":reduced.rule_applications,"terms":reduced.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"unresolved":reduced.unresolved.iter().map(|t|serde_json::json!({"indices":t.integral.to_vec(),"coefficient":t.coefficient.to_expression().to_string(),"reason":format!("{:?}",t.reason)})).collect::<Vec<_>>(),"conditions":reduced.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()}));
        points.push(serde_json::json!({"selection":selection,"reduction":red_json}));std::fs::write(out.join("points-progress.json"),serde_json::to_vec_pretty(&points).unwrap()).unwrap();
    }
    let result=serde_json::json!({"scope":"Read-only bounded actual final-frontier and stratified historical-label reductions from the saved round18 provisional program. Complete cold source proof replay precedes terminal binding. This diagnostic neither establishes final closure nor changes production candidate collection. No source discovery or period values.","rules":count,"source_identity":record.measure,"timings":timings,"points":points,"limits":{"max_rule_applications":work.max_rule_applications,"max_pending_integrals":work.max_pending_integrals},"timing_scope":"Single pass shared host; CPU tick granularity retained; conversion measured separately from native reduction; no whole-workflow speed claim."});
    std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();progress(&out,"complete",&timings);
}
