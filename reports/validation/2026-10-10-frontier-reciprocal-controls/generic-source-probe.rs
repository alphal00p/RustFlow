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

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let bytes = std::fs::read(&args[1]).unwrap();
    let limits = BinaryIoLimits::default();
    let envelope = inspect_program(&bytes, limits).unwrap();
    let (record, used): (Rec, usize) = bincode::decode_from_slice(
        envelope.section(SectionTag::PROGRAM).unwrap(),
        bincode::config::standard(),
    )
    .unwrap();
    assert_eq!(used, envelope.section(SectionTag::PROGRAM).unwrap().len());
    assert!(
        matches!(
            record.schema.as_str(),
            "rustred.guarded-source-program.v1" | "rustred.guarded-source-program.v2"
        ),
        "only explicitly known native source schemas may be imported"
    );
    let table = DecodedCoefficientTable::import_generated_normalized(
        envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),
        envelope.section(SectionTag::COEFFICIENTS).unwrap(),
        limits,
    )
    .unwrap();
    let polynomial = |id| {
        let c = table
            .coefficient(CoefficientId::try_from_index(id).unwrap())
            .unwrap()
            .clone();
        assert!(c.denominator.is_one());
        c.numerator
    };
    let roles: [IndexRole; 12] = record
        .roles
        .iter()
        .map(|r| match r {
            0 => IndexRole::Ordinary,
            1 => IndexRole::RequiredCut,
            2 => IndexRole::Occupation,
            _ => panic!(),
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let indices: [usize; 12] = record.indices.clone().try_into().unwrap();
    let mode = &args[4];
    assert_eq!(
        mode, "original",
        "portfolio control must retain full original corpus"
    );
    let order: Vec<usize> = (0..record.sources.len()).collect();
    let rows = order
        .iter()
        .map(|&i| {
            let s = &record.sources[i];
            GuardedSource::new(
                s.id.clone(),
                s.terms
                    .iter()
                    .map(|t| Term {
                        integral: label(&t.integral),
                        coefficient: polynomial(t.coefficient),
                    })
                    .collect(),
                domain(&s.domain),
            )
            .with_nonzero_conditions(s.conditions.iter().map(|&id| polynomial(id)).collect())
        })
        .collect();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            format!(
                "{};singleton-source-order-diagnostic={mode}",
                record.measure
            ),
            roles,
            indices,
            rows,
        )
        .unwrap()
        .with_zero_domains(record.zero_domains.iter().map(domain).collect())
        .unwrap(),
    );
    if let Some(other) = args.get(5) {
        let encoded = std::fs::read(other).unwrap();
        assert!(
            GuardedProgram::decode_generated(&encoded, sources, limits).is_err(),
            "cross-experiment program accepted"
        );
        println!("{{\"cross_schema_rejected\":true}}");
        return;
    }
    let proof_dir = std::path::Path::new(&args[3]).with_extension("proofs");
    std::fs::create_dir_all(&proof_dir).unwrap();
    let raw: Vec<Vec<i64>> = serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();
    let points: Vec<[i64; 12]> = raw.into_iter().map(|p| p.try_into().unwrap()).collect();
    let mut probes = Vec::new();
    for (point_ordinal, target) in points.into_iter().enumerate() {
        let started = std::time::Instant::now();
        let found = sources
            .solve_domains(
                vec![IndexDomain::new(target.map(IndexBounds::fixed)).unwrap()],
                SearchOptions {
                    max_depth: Some(3),
                    sample_seed: 0,
                    ..Default::default()
                },
                1,
            )
            .unwrap();
        let gaps = found
            .unresolved
            .iter()
            .map(|g| serde_json::json!({"reason":format!("{:?}",g.reason),"detail":g.detail}))
            .collect::<Vec<_>>();
        let program = GuardedProgram::new(sources.clone(), found.rules, []).unwrap();
        let encoded = program.encode_native(limits).unwrap();
        let proof_path = proof_dir.join(format!("point-{point_ordinal:02}.bin"));
        std::fs::write(&proof_path, &encoded).unwrap();
        let decoded = GuardedProgram::decode_generated(&encoded, sources.clone(), limits).unwrap();
        let applied = program.apply(&target).unwrap();
        let replayed = decoded.apply(&target).unwrap();
        assert_eq!(
            format!("{:?}", applied.status),
            format!("{:?}", replayed.status)
        );
        assert_eq!(applied.terms, replayed.terms);
        assert_eq!(applied.nonzero_conditions, replayed.nonzero_conditions);
        probes.push(serde_json::json!({"target":target.to_vec(),"proof_path":proof_path,"persistence_roundtrip_replay":true,"search_stats":program.rules().iter().map(|r|serde_json::json!({"rows":r.candidate().stats.rows,"seeds":r.candidate().stats.seeds,"independent_rows":r.candidate().stats.independent_rows,"exact_trace_rows":r.candidate().stats.exact_trace_rows,"direct_hit":r.candidate().stats.direct_hit,"portfolio":experiment_stats(r)})).collect::<Vec<_>>(),"mode":mode,"seconds":started.elapsed().as_secs_f64(),"status":format!("{:?}",applied.status),"rhs":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"discovery_gaps":gaps,"rule_sources":program.rules().iter().flat_map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":record.sources[order[s.basis_row]].id,"original_ordinal":order[s.basis_row],"seed":(0..12).map(|i|s.seed.integral[i].value()).collect::<Vec<_>>()}))).collect::<Vec<_>>() }));
    }
    let sources_json=order.iter().map(|&i|{let s=&record.sources[i];serde_json::json!({"original_ordinal":i,"source_id":s.id,"domain":s.domain,"conditions":s.conditions.iter().map(|&id|polynomial(id).to_expression().to_string()).collect::<Vec<_>>(),"terms":s.terms.iter().map(|t|serde_json::json!({"label":t.integral,"coefficient":polynomial(t.coefficient).to_expression().to_string()})).collect::<Vec<_>>()})}).collect::<Vec<_>>();
    let output = serde_json::json!({"scope":"Isolated point-only original-source portfolio: complete original source corpus retained, fresh native discovery, exact replay and generated-program roundtrip. No physical values or full closure.", "experiment_enabled":cfg!(source_portfolio), "historical_rules_imported":false, "mode":mode,"source_program":args[1],"imported_source_schema":record.schema,"original_measure_id":record.measure,"roles":record.roles,"index_variable_positions":record.indices,"zero_domains":record.zero_domains,"point_fixture":args[2],"points":probes.len(),"max_domains_per_point":1,"max_depth":3,"full_corpus_preserved":mode=="original"||mode=="virtual-first","original_source_count":record.sources.len(),"source_count":order.len(),"source_term_count":order.iter().map(|&i|record.sources[i].terms.len()).sum::<usize>(),"probes":probes,"original_source_rows":sources_json});
    std::fs::write(&args[3], serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    println!(
        "{}",
        serde_json::json!({"mode":mode,"points":output["points"],"source_count":output["source_count"],"applied":output["probes"].as_array().unwrap().iter().filter(|p|p["status"].as_str().unwrap().starts_with("Applied")).count(),"output":args[3]})
    );
}
