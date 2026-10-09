//! Explicit source/closure pilot for runtime graph inputs. This test does not
//! admit a physical contour, construct a boundary, or evaluate an amplitude.
use std::path::PathBuf;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::guarded::{
    GuardedDiscoveryOptions, GuardedMeasureIdentity, GuardedUnresolved, IndexDomain, IndexRole,
};
use symbolica_amflow::finite_density::preparation::{WeightedSourceOptions, WeightedSourcePolicy};
use symbolica_amflow::finite_density::reduction::{
    GuardRefinementOptions, WeightedClosureOptions, WeightedClosureOutcome, prepare_weighted_system,
};
use symbolica_amflow::{Result, RunContext};

fn budget(name: &str, default: usize) -> usize {
    std::env::var(name)
        .map(|value| value.parse().expect("invalid diagnostic budget"))
        .unwrap_or(default)
}
fn enabled(name: &str, default: bool) -> bool {
    match std::env::var(name).as_deref() {
        Ok("1") | Ok("true") => true,
        Ok("0") | Ok("false") => false,
        Err(_) => default,
        Ok(value) => panic!("invalid boolean {name}={value}"),
    }
}
fn bounds<const N: usize>(domain: &IndexDomain<N>) -> Vec<[Option<i64>; 2]> {
    domain
        .bounds()
        .iter()
        .map(|bound| [bound.lower(), bound.upper()])
        .collect()
}
fn gaps<const N: usize>(unresolved: &[GuardedUnresolved<N>]) -> Vec<serde_json::Value> {
    unresolved.iter().map(|gap| serde_json::json!({"domain":bounds(&gap.domain),"reason":format!("{:?}",gap.reason),"detail":gap.detail})).collect()
}
fn save(report: &std::path::Path, name: &str, value: &serde_json::Value) {
    std::fs::write(report.join(name), serde_json::to_vec_pretty(value).unwrap()).unwrap();
}
fn pilot<const N: usize>(input: DensityInput, cuts: &[usize], root: PathBuf) -> Result<()> {
    let started = Instant::now();
    let family = input
        .prepare()?
        .occupied_cut(cuts, 65536)?
        .at_physical_masses();
    let physical_arity = family.factors().len();
    let shifted = (0..family.physical_slots())
        .filter(|&slot| family.roles()[slot] == IndexRole::Ordinary)
        .collect::<Vec<_>>();
    let policy = std::env::var("RUSTFLOW_WEIGHTED_SOURCE_POLICY")
        .map(|value| value.parse::<WeightedSourcePolicy>().unwrap())
        .unwrap_or_default();
    let source_options = WeightedSourceOptions {
        policy,
        positive_compact_energy_powers: enabled("RUSTFLOW_WEIGHTED_POSITIVE_ENERGY_POWERS", false),
        free_virtual_zero_sectors: false,
    };
    let report = root.join(format!(
        "physical-{physical_arity}-capacity-{N}-cuts-{}",
        cuts.iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_")
    ));
    std::fs::create_dir_all(&report)?;
    let options = WeightedClosureOptions {
        max_rounds: budget("RUSTFLOW_WEIGHTED_ROUNDS", 2),
        max_frontier: budget("RUSTFLOW_WEIGHTED_FRONTIER", 128),
        max_requested: budget("RUSTFLOW_WEIGHTED_REQUESTED", 512),
        discovery: GuardedDiscoveryOptions {
            max_depth: budget("RUSTFLOW_WEIGHTED_DEPTH", 2).try_into().unwrap(),
            max_domains: budget("RUSTFLOW_WEIGHTED_DOMAINS", 256),
            sample_seed: 0,
        },
        guard_refinement: GuardRefinementOptions {
            max_passes: budget("RUSTFLOW_WEIGHTED_GUARD_PASSES", 0),
            max_added_domains: budget("RUSTFLOW_WEIGHTED_GUARD_DOMAINS", 64),
            max_interval_width: budget("RUSTFLOW_WEIGHTED_GUARD_WIDTH", 2)
                .try_into()
                .unwrap(),
        },
        search_frontier_sectors: enabled("RUSTFLOW_WEIGHTED_FRONTIER_SECTORS", true),
        split_ordinary_zero_faces: enabled("RUSTFLOW_WEIGHTED_ZERO_FACES", false),
        checkpoints: Some(report.clone()),
        ..Default::default()
    };
    save(
        &report,
        "input.json",
        &serde_json::json!({"input":input,"cut_slots":cuts,"shifted_slots":shifted,"physical_arity":physical_arity,"storage_capacity":N,"source_policy":policy.as_str(),"positive_compact_energy_powers":source_options.positive_compact_energy_powers,"free_virtual_zero_sectors":source_options.free_virtual_zero_sectors,"source_domain_budget":budget("RUSTFLOW_WEIGHTED_SOURCE_DOMAIN_BUDGET",65536),"max_rounds":options.max_rounds,"max_frontier":options.max_frontier,"max_requested":options.max_requested,"max_depth":options.discovery.max_depth,"max_domains":options.discovery.max_domains,"sample_seed":options.discovery.sample_seed,"guard_refinement":options.guard_refinement,"search_frontier_sectors":options.search_frontier_sectors,"split_ordinary_zero_faces":options.split_ordinary_zero_faces,"numerical_prediction":false,"contour_admission":false}),
    );
    let prepared=family.guarded_sources_with_options::<N>(symbol!("generic_closure::epsilon"),4,symbol!("generic_closure::eta"),&shifted,budget("RUSTFLOW_WEIGHTED_SOURCE_DOMAIN_BUDGET",65536),vec![],GuardedMeasureIdentity {
        measure:format!("runtime graph input={}; assigned factors={:?}; coordinates={:?}",input.prepare()?.identity(),family.factors(),family.coordinates()),
        support:format!("formal meromorphic weighted source identities on occupied shells {:?}; no endpoint-contour certificate",family.shells()),
        orientation:format!("future occupied momenta; inverse routing={:?}; determinant={}",family.inverse_routing(),family.routing_determinant()),
        normalization:"unscaled source identities; original target Wick phases retained; no numerical normalization or oracle inputs".into(),
        branch:"formal exact rational eta,epsilon coefficients with all nonzero conditions retained; no endpoint or dimensional sample admission".into(),
        deformation:format!("native D-eta on physical slots {shifted:?}; cuts and occupation weights fixed"),
    },source_options)?;
    let sources = prepared.context.sources();
    save(
        &report,
        "sources.json",
        &serde_json::json!({"measure_id":sources.measure_id(),"physical_arity":prepared.context.physical_arity(),"storage_capacity":N,"roles":sources.roles().iter().map(|role|format!("{role:?}")).collect::<Vec<_>>(),"coefficient_variables":format!("{:?}",sources.native_sources().coefficient_variables()),"source_count":sources.sources().len(),"sources":sources.sources().iter().zip(sources.native_sources().rows()).map(|(info,row)|serde_json::json!({"id":info.id,"domain":bounds(&info.domain),"nonzero_conditions":info.nonzero_conditions.iter().map(|c|c.to_expression().to_canonical_string()).collect::<Vec<_>>(),"terms":row.iter().map(|term|serde_json::json!({"powers":term.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"coefficient":term.coefficient.to_expression().to_canonical_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>(),"zero_domains":sources.zero_domains().iter().map(bounds).collect::<Vec<_>>(),"targets":prepared.targets.iter().map(|target|target.iter().map(|(indices,coefficient)|serde_json::json!({"indices":indices.to_vec(),"coefficient":coefficient.to_canonical_string()})).collect::<Vec<_>>()).collect::<Vec<_>>()}),
    );
    let preparation_seconds = started.elapsed().as_secs_f64();
    eprintln!(
        "generic physical={physical_arity} capacity={N} sources={}: bounded native closure starting",
        sources.sources().len()
    );
    let run = RunContext {
        progress: Some(std::sync::Arc::new(move |event| {
            eprintln!("generic capacity={N}: {event:?}")
        })),
        ..Default::default()
    };
    let closure_started = Instant::now();
    let result = prepare_weighted_system(
        &prepared.context,
        &prepared.targets,
        &prepared.deformation,
        options,
        &run,
    );
    let closure_seconds = closure_started.elapsed().as_secs_f64();
    let outcome = match result {
        Ok(outcome) => outcome,
        Err(error) => {
            save(
                &report,
                "result.json",
                &serde_json::json!({"status":"error","error":error.to_string(),"preparation_seconds":preparation_seconds,"closure_seconds":closure_seconds,"physical_arity":physical_arity,"storage_capacity":N,"numerical_prediction":false}),
            );
            return Err(error);
        }
    };
    let (status, diagnostics, detail) = match &outcome {
        WeightedClosureOutcome::Closed(closed) => {
            if let Some(system) = closed.differential_system() {
                system.validate()?;
            }
            (
                "closed",
                &closed.diagnostics,
                serde_json::json!({"basis":closed.reduced.basis.iter().map(|i|&i.0).collect::<Vec<_>>(),"matrix":closed.reduced.matrix.iter().map(|row|row.iter().map(|a|a.to_canonical_string()).collect::<Vec<_>>()).collect::<Vec<_>>(),"target_weights":closed.reduced.targets.iter().map(|target|target.iter().map(|(i,c)|serde_json::json!({"indices":i.0,"coefficient":c.to_canonical_string()})).collect::<Vec<_>>()).collect::<Vec<_>>(),"nonzero_conditions":closed.reduced.nonzero_conditions.iter().map(|a|a.to_canonical_string()).collect::<Vec<_>>(),"discovery_unresolved":gaps(&closed.discovery_unresolved)}),
            )
        }
        WeightedClosureOutcome::Unresolved(unresolved) => (
            "unresolved",
            &unresolved.diagnostics,
            serde_json::json!({"reason":unresolved.reason,"frontier":unresolved.provisional_frontier.iter().map(|i|i.to_vec()).collect::<Vec<_>>(),"application_unresolved":unresolved.unresolved.iter().map(|term|serde_json::json!({"indices":term.integral.to_vec(),"coefficient":term.coefficient.to_canonical_string(),"reason":format!("{:?}",term.reason)})).collect::<Vec<_>>(),"discovery_unresolved":gaps(&unresolved.discovery_unresolved),"nonzero_conditions":unresolved.nonzero_conditions.iter().map(|a|a.to_canonical_string()).collect::<Vec<_>>()}),
        ),
    };
    save(
        &report,
        "result.json",
        &serde_json::json!({"status":status,"preparation_seconds":preparation_seconds,"closure_seconds":closure_seconds,"physical_arity":physical_arity,"storage_capacity":N,"rounds":diagnostics.rounds,"requested":diagnostics.requested,"native_frontier_requests":diagnostics.native_frontier_requests,"provisional_sizes":diagnostics.provisional_sizes,"native_rules":diagnostics.native_rules,"native_rule_applications":diagnostics.native_rule_applications,"guard_refinement_passes":diagnostics.guard_refinement_passes,"guard_refinement_added_domains":diagnostics.guard_refinement_added_domains,"guard_refinement_budget_exhausted":diagnostics.guard_refinement_budget_exhausted,"guard_refinements":diagnostics.guard_refinements,"detail":detail,"numerical_prediction":false,"contour_admission":false}),
    );
    eprintln!(
        "generic physical={physical_arity} capacity={N}: {status}; closure {closure_seconds:.3}s; report={}",
        report.display()
    );
    Ok(())
}

#[test]
#[ignore = "bounded runtime graph source pilot; explicitly inspect closed/unresolved report; no numerical acceptance"]
fn runtime_graph_native_closure() {
    let input_path = std::env::var_os("RUSTFLOW_WEIGHTED_INPUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from("examples/finite_density/chain_of_three_parallel_pairs.json")
        });
    let input: DensityInput = serde_json::from_slice(&std::fs::read(input_path).unwrap()).unwrap();
    let cuts = std::env::var("RUSTFLOW_WEIGHTED_CUTS")
        .unwrap_or_else(|_| "0".into())
        .split(',')
        .map(|value| value.trim().parse::<usize>().unwrap())
        .collect::<Vec<_>>();
    let physical = input
        .prepare()
        .unwrap()
        .occupied_cut(&cuts, 65536)
        .unwrap()
        .factors()
        .len();
    let capacity = budget(
        "RUSTFLOW_WEIGHTED_CAPACITY",
        if physical <= 16 {
            16
        } else if physical <= 20 {
            20
        } else {
            24
        },
    );
    let root = std::env::var_os("RUSTFLOW_WEIGHTED_CLOSURE_REPORT")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-generic-weighted-closure"));
    match capacity {
        16 => pilot::<16>(input, &cuts, root),
        20 => pilot::<20>(input, &cuts, root),
        24 => pilot::<24>(input, &cuts, root),
        other => panic!("diagnostic compiled capacities are 16,20,24; requested {other}"),
    }
    .unwrap();
}
