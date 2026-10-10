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
    let roles: [IndexRole; 16] = record
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
    let indices: [usize; 16] = record.indices.clone().try_into().unwrap();
    let mode = &args[4];
    assert!(["original", "global-euler"].contains(&mode.as_str()));
    let order: Vec<usize> = (0..record.sources.len()).collect();
    let rows: Vec<GuardedSource<16>> = order
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
        GuardedSourceSystem::new(record.measure.clone(), roles, indices, rows)
            .unwrap()
            .with_zero_domains(record.zero_domains.iter().map(domain).collect())
            .unwrap(),
    );
    let depth: u32 = args[5].parse().unwrap();
    assert!([3, 4, 5].contains(&depth));
    let selected: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();
    let selected = selected.as_array().unwrap();
    assert!((1..=16).contains(&selected.len()));
    let physical = if roles[12] == IndexRole::Occupation {
        13
    } else {
        11
    };
    let mut bounds0 = *IndexDomain::for_roles(&roles).bounds();
    for b in &mut bounds0[5..9] {
        *b = IndexBounds::new(None, Some(0)).unwrap()
    }
    for b in &mut bounds0[physical..] {
        *b = IndexBounds::fixed(0)
    }
    let admitted = IndexDomain::new(bounds0).unwrap();
    use rustred::solver::guarded::{
        GuardedApplicationFailure as Failure, GuardedApplicationStatus as Status,
    };
    let historical = GuardedProgram::decode_generated(&bytes, sources.clone(), limits).unwrap();
    assert!(historical.terminals().is_empty());
    let mut observations = Vec::new();
    for point in selected {
        let target: [i64; 16] = serde_json::from_value(point["point"].clone()).unwrap();
        let applied = historical.apply(&target).unwrap();
        if point["kind"] == "actual-frontier" {
            assert_eq!(
                applied.status,
                Status::Unresolved(Failure::NoApplicableRule)
            )
        } else {
            assert!(matches!(applied.status, Status::Applied { .. }))
        }
        observations.push(serde_json::json!({"point":target,"kind":point["kind"],"status":format!("{:?}",applied.status)}));
    }
    drop(historical);
    assert_eq!(record.sources.len(), 232);
    let source_epoch = record.measure.clone();
    let variables = sources.native_sources().coefficient_variables();
    let symbols = variables
        .iter()
        .map(|v| {
            let PolyVariable::Symbol(s) = v else { panic!() };
            *s
        })
        .collect::<Vec<_>>();
    let index_symbols: [Symbol; 16] = record
        .indices
        .iter()
        .map(|&i| symbols[i])
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let parameters = symbols
        .iter()
        .enumerate()
        .filter(|(i, _)| !record.indices.contains(i))
        .map(|(_, s)| *s)
        .collect::<Vec<_>>();
    let eta = symbol!("rustflow_occupied::eta");
    let epsilon = symbol!("rustflow_occupied::epsilon");
    let input: symbolica_amflow::finite_density::DensityInput =
        serde_json::from_slice(&std::fs::read(&args[6]).unwrap()).unwrap();
    let input = input.prepare().unwrap();
    assert!(record.measure.contains(input.identity()));
    let family = input
        .occupied_cut(&[0, 3], 4096)
        .unwrap()
        .at_physical_masses();
    let measure = family.deformed_measure::<16>(eta, &[1, 2, 4]).unwrap();
    assert_eq!(measure.roles(), &roles);
    assert_eq!(family.loops(), 3);
    assert_eq!(family.input_slots(), 9);
    assert_eq!(physical, 13);
    let coordinates = family.coordinates();
    let dimension = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let direction = coordinates
        .iter()
        .enumerate()
        .map(|(i, q)| Atom::num(if i < 6 { 2 } else { 1 }) * q)
        .collect::<Vec<_>>();
    let mut added = measure
        .ibp(
            "global-euler/3",
            &index_symbols,
            &direction,
            &(Atom::num(3) * &dimension),
            16,
        )
        .unwrap();
    assert_eq!(added.len(), 16);
    for row in &mut added {
        row.domain = row.domain.intersection(&admitted).unwrap();
    }
    let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
    let locals = (0..3)
        .map(|leg| {
            let mut v = pairs
                .iter()
                .enumerate()
                .map(|(k, &(i, j))| {
                    Atom::num(i32::from(i == leg) + i32::from(j == leg)) * &coordinates[k]
                })
                .collect::<Vec<_>>();
            v.extend((0..3).map(|i| Atom::num(i32::from(i == leg)) * &coordinates[6 + i]));
            measure
                .ibp(
                    &format!("audit-local/{leg}"),
                    &index_symbols,
                    &v,
                    &dimension,
                    16,
                )
                .unwrap()
        })
        .collect::<Vec<_>>();
    for global in &added {
        let witness = std::array::from_fn(|i| {
            let b = global.domain.bounds()[i];
            b.lower().unwrap_or(0).min(b.upper().unwrap_or(i64::MAX))
        });
        assert!(global.domain.contains(&witness));
        let mut difference: BTreeMap<[i16; 16], Atom> = BTreeMap::new();
        for term in &global.terms {
            *difference.entry(term.shift).or_default() += &term.coefficient;
        }
        for local in &locals {
            let rows = local
                .iter()
                .filter(|row| row.domain.contains(&witness))
                .collect::<Vec<_>>();
            assert_eq!(rows.len(), 1);
            for term in &rows[0].terms {
                *difference.entry(term.shift).or_default() -= &term.coefficient;
            }
        }
        assert!(
            difference
                .values()
                .all(|c| c.expand().together().cancel().is_zero())
        );
    }
    let added_audit=added.iter().map(|s|serde_json::json!({"id":s.id,"domain":bounds(&s.domain),"conditions":s.nonzero_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),"terms":s.terms.iter().map(|t|serde_json::json!({"shift":t.shift,"coefficient":t.coefficient.to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
    let sources = if mode == "original" {
        sources
    } else {
        use symbolica_amflow::finite_density::guarded::{GuardedContext, GuardedMeasureIdentity};
        let ctx = GuardedContext::new_with_physical_arity(
            GuardedMeasureIdentity {
                measure: "validation global Euler exact sum of original local dilations".into(),
                support: "original polynomial and distribution domains".into(),
                orientation: "unchanged occupied cut0,3".into(),
                normalization: "homogeneous source only".into(),
                branch: "unchanged joint origin and zeros".into(),
                deformation: "original shifted1,2,4".into(),
            },
            roles,
            index_symbols,
            parameters,
            added,
            13,
        )
        .unwrap();
        assert_eq!(
            ctx.sources().native_sources().coefficient_variables(),
            variables
        );
        let mut rows = sources
            .sources()
            .iter()
            .zip(sources.native_sources().rows())
            .map(|(m, r)| {
                GuardedSource::new(m.id.clone(), r.clone(), m.domain.clone())
                    .with_nonzero_conditions(m.nonzero_conditions.clone())
            })
            .collect::<Vec<_>>();
        rows.extend(
            ctx.sources()
                .sources()
                .iter()
                .zip(ctx.sources().native_sources().rows())
                .map(|(m, r)| {
                    GuardedSource::new(m.id.clone(), r.clone(), m.domain.clone())
                        .with_nonzero_conditions(m.nonzero_conditions.clone())
                }),
        );
        Arc::new(
            GuardedSourceSystem::new(
                format!("{};global-euler-presentation-v1", source_epoch),
                roles,
                indices,
                rows,
            )
            .unwrap()
            .with_zero_domains(record.zero_domains.iter().map(domain).collect())
            .unwrap(),
        )
    };
    let out = std::path::PathBuf::from(&args[3]);
    std::fs::create_dir_all(&out).unwrap();
    let describe = |applied: &rustred::solver::guarded::GuardedApplication<16>| serde_json::json!({"status":format!("{:?}",applied.status),"rhs":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()});
    let fallback = |s: &Status| {
        matches!(
            s,
            Status::Unresolved(Failure::NoApplicableRule | Failure::ConditionVanished { .. })
        )
    };
    let mut probes = Vec::new();
    for (ordinal, selection) in selected.iter().enumerate() {
        let target: [i64; 16] = serde_json::from_value(selection["point"].clone()).unwrap();
        assert!(admitted.contains(&target));
        let search = |domain| {
            sources
                .solve_domains(
                    vec![domain],
                    SearchOptions {
                        max_depth: Some(depth),
                        sample_seed: 0,
                        ..Default::default()
                    },
                    1,
                )
                .unwrap()
        };
        let ray = IndexDomain::new(std::array::from_fn(|i| {
            if target[i] == 0 {
                IndexBounds::fixed(0)
            } else if target[i] > 0 {
                IndexBounds::new(Some(target[i]), None).unwrap()
            } else {
                IndexBounds::new(None, Some(target[i])).unwrap()
            }
        }))
        .unwrap()
        .intersection(&admitted)
        .unwrap();
        let ray = search(ray);
        let mut gaps=ray.unresolved.iter().map(|g|serde_json::json!({"phase":"ray","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
        let mut program = GuardedProgram::new(sources.clone(), ray.rules, []).unwrap();
        let ray_result = program.apply(&target).unwrap();
        let needs_point = fallback(&ray_result.status);
        assert!(needs_point || matches!(ray_result.status, Status::Applied { .. } | Status::Zero));
        let mut point_result = serde_json::Value::Null;
        if needs_point {
            let exact = search(IndexDomain::new(target.map(IndexBounds::fixed)).unwrap());
            gaps.extend(exact.unresolved.iter().map(|g|serde_json::json!({"phase":"point","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})));
            let exact = GuardedProgram::new(sources.clone(), exact.rules, []).unwrap();
            point_result = describe(&exact.apply(&target).unwrap());
            program = program.union_replayed(exact, [], 65536).unwrap();
        }
        let encoded = program.encode_native(limits).unwrap();
        std::fs::write(out.join(format!("point-{ordinal:02}.bin")), &encoded).unwrap();
        let replay = GuardedProgram::decode_generated(&encoded, sources.clone(), limits).unwrap();
        let result = program.apply(&target).unwrap();
        let check = replay.apply(&target).unwrap();
        assert_eq!(result.status, check.status);
        assert_eq!(result.terms, check.terms);
        assert_eq!(result.nonzero_conditions, check.nonzero_conditions);
        assert!(result.terms.keys().all(|p| admitted.contains(p)));
        probes.push(serde_json::json!({"selection":selection,"ray":describe(&ray_result),"point_fallback":needs_point,"point":point_result,"application":describe(&result),"roundtrip":true,"rules":program.rules().len(),"gaps":gaps,"proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>()})).collect::<Vec<_>>()).collect::<Vec<_>>() }));
    }
    let report = serde_json::json!({"scope":"Bounded depth control on exact historical NoApplicableRule leaves, with separately verified successful root control. Historical program is replayed ONLY to verify selection status, then discarded before fresh discovery; no old rule union. Original232 sources retained in order; optional16 global Euler branches generated by exact summed local dilations. No periods/closure inference.","depth":depth,"ray_domains":1,"point_domains":1,"mode":mode,"source_count":sources.sources().len(),"source_identity":sources.measure_id(),"original_source_identity":record.measure,"added_source_rows":added_audit,"global_is_exact_sum_of_local_dilations":true,"same_original_domains_and_zeros":true,"source_schema":record.schema,"selection_observations":observations,"all_rhs_labels_admitted":true,"probes":probes});
    std::fs::write(
        out.join("result.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!(
        "depth={depth} points={} source={} selection and fresh replay passed",
        selected.len(),
        record.sources.len()
    );
}
