//! Standalone matched partial-placement source control, no flow admission.
use rustred::persistence::BinaryIoLimits;
use rustred::solver::SearchOptions;
use rustred::solver::guarded::{
    GuardedApplicationFailure as Failure, GuardedApplicationStatus as Status, GuardedProgram,
    GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain,
};
use std::collections::BTreeSet;
use std::sync::Arc;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::guarded::GuardedMeasureIdentity;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::preparation::{WeightedSourceOptions, WeightedSourcePolicy};

const N: usize = 16;
fn bounds(d: &IndexDomain<N>) -> Vec<(Option<i64>, Option<i64>)> {
    d.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let out = std::path::PathBuf::from(&args[2]);
    assert!(!out.exists());
    std::fs::create_dir_all(&out).unwrap();
    let mode = &args[3];
    let (cuts, shifted, all): (Vec<usize>, Vec<usize>, Vec<usize>) = match mode.as_str() {
        "single-all" => (vec![3], vec![0, 1, 2, 4], vec![0, 1, 2, 4]),
        "single-slot0" => (vec![3], vec![0], vec![0, 1, 2, 4]),
        "double-all" => (vec![0, 3], vec![1, 2, 4], vec![1, 2, 4]),
        "double-slots12" => (vec![0, 3], vec![1, 2], vec![1, 2, 4]),
        "double-slots24" => (vec![0, 3], vec![2, 4], vec![1, 2, 4]),
        _ => panic!("unknown mode"),
    };
    let proof_binding = &args[4];
    let input: DensityInput = serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();
    let input = input.prepare().unwrap();
    let family = input
        .occupied_cut(&cuts, 4096)
        .unwrap()
        .at_physical_masses();
    assert_eq!(family.loops(), 3);
    assert_eq!(family.physical_slots(), 5);
    assert_eq!(family.input_slots(), 9);
    assert!(
        family
            .shells()
            .iter()
            .all(|s| s.mass_squared.is_zero() && s.chemical_potential > Rational::zero())
    );
    let eta = symbol!("partial_placement_probe::eta");
    let epsilon = symbol!("partial_placement_probe::epsilon");
    let identity = GuardedMeasureIdentity {
        measure: format!("validation-partial-placement-v1; input={}; cuts={cuts:?}; mask={shifted:?}",input.identity()),
        support: "polynomial completion<=0; full-positive physical support lower-origin zeros only; no free-virtual zeros".into(),
        orientation: "future occupied shells from exact input geometry".into(),
        normalization: "native homogeneous weighted identities, no values".into(),
        branch: format!("full-positive-support resolved germ; line regulators at fixed eta/T first; joint high-D finite jets; proof SHA256 {proof_binding}"),
        deformation: format!("fixed shells/occupations; selected uncut native factors minus eta: {shifted:?}"),
    };
    let prepared = family
        .guarded_sources_with_options::<N>(
            epsilon,
            4,
            eta,
            &shifted,
            4096,
            vec![],
            identity,
            WeightedSourceOptions {
                policy: WeightedSourcePolicy::PolynomialClosure,
                positive_compact_energy_powers: false,
                free_virtual_zero_sectors: false,
            },
        )
        .unwrap();
    let admitted = prepared.deformation.admitted_domain().clone();
    let original = prepared.context.sources();
    assert!(original.zero_domains().is_empty());
    assert!(!original.sources().iter().any(|s| s.id.contains("raw-ward")));
    let mut zeros = Vec::new();
    for shell in family.shells() {
        let mut b = *admitted.bounds();
        for bound in &mut b[..family.physical_slots()] {
            *bound = IndexBounds::new(Some(1), None).unwrap()
        }
        b[shell.lower_slot] = IndexBounds::new(Some(1), None).unwrap();
        zeros.push(IndexDomain::new(b).unwrap());
    }
    let rows = original
        .sources()
        .iter()
        .zip(original.native_sources().rows())
        .map(|(s, r)| {
            GuardedSource::new(s.id.clone(), r.clone(), s.domain.clone())
                .with_nonzero_conditions(s.nonzero_conditions.clone())
        })
        .collect();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            original.measure_id(),
            *original.roles(),
            *original.native_sources().index_variables(),
            rows,
        )
        .unwrap()
        .with_zero_domains(zeros.clone())
        .unwrap(),
    );
    // All matched masks query the same initial roots and all-uncut raw derivatives.
    let roots = prepared
        .targets
        .iter()
        .flat_map(|t| t.keys().copied())
        .collect::<BTreeSet<_>>();
    assert!(!roots.is_empty() && roots.len() <= 16);
    let mut derivatives = BTreeSet::new();
    for root in &roots {
        for &slot in &all {
            if root[slot] != 0 {
                let mut p = *root;
                p[slot] += 1;
                derivatives.insert(p);
            }
        }
    }
    let mut selected = roots.clone();
    for p in &derivatives {
        if selected.len() >= 16 {
            break;
        }
        selected.insert(*p);
    }
    let selection = selected.iter().copied().collect::<Vec<_>>();
    let limits = BinaryIoLimits::default();
    let describe = |a: &rustred::solver::guarded::GuardedApplication<N>| {
        serde_json::json!({
        "status":format!("{:?}",a.status),"rhs":a.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),
        "conditions":a.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()})
    };
    let mut probes = Vec::new();
    for (ordinal, target) in selection.iter().enumerate() {
        assert!(admitted.contains(target));
        let search = |domain| {
            sources
                .solve_domains(
                    vec![domain],
                    SearchOptions {
                        max_depth: Some(3),
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
        let found = search(ray);
        let mut gaps=found.unresolved.iter().map(|g|serde_json::json!({"phase":"ray","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
        let mut program = GuardedProgram::new(sources.clone(), found.rules, []).unwrap();
        let ray_result = program.apply(target).unwrap();
        let fallback = matches!(
            ray_result.status,
            Status::Unresolved(Failure::NoApplicableRule | Failure::ConditionVanished { .. })
        );
        assert!(fallback || matches!(ray_result.status, Status::Applied { .. } | Status::Zero));
        let mut point_result = serde_json::Value::Null;
        if fallback {
            let found = search(IndexDomain::new(target.map(IndexBounds::fixed)).unwrap());
            gaps.extend(found.unresolved.iter().map(|g|serde_json::json!({"phase":"point","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})));
            let exact = GuardedProgram::new(sources.clone(), found.rules, []).unwrap();
            point_result = describe(&exact.apply(target).unwrap());
            program = program.union_replayed(exact, [], 65536).unwrap();
        }
        let encoded = program.encode_native(limits).unwrap();
        std::fs::write(out.join(format!("point-{ordinal:02}.bin")), &encoded).unwrap();
        let replay = GuardedProgram::decode_generated(&encoded, sources.clone(), limits).unwrap();
        let result = program.apply(target).unwrap();
        let check = replay.apply(target).unwrap();
        assert_eq!(result.status, check.status);
        assert_eq!(result.terms, check.terms);
        assert_eq!(result.nonzero_conditions, check.nonzero_conditions);
        assert!(result.terms.keys().all(|p| admitted.contains(p)));
        probes.push(serde_json::json!({"point":target.to_vec(),"root":roots.contains(target),"ray":describe(&ray_result),"point_fallback":fallback,"point_result":point_result,"application":describe(&result),"roundtrip":true,"rule_count":program.rules().len(),"gaps":gaps,
            "proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|sources.sources()[s.basis_row].id.clone()).collect::<Vec<_>>()).collect::<Vec<_>>() }));
    }
    let measure = family.deformed_measure::<N>(eta, &shifted).unwrap();
    let report = serde_json::json!({"scope":"Fresh matched restricted source contexts only; no numerical/closure/admission claim. Public PolynomialClosure without sealed origin omits raw Ward for ALL arms. No old rules or zero boxes imported.",
        "mode":mode,"cuts":cuts,"shifted":shifted,"all_uncut":all,"input_identity":input.identity(),"source_identity":sources.measure_id(),"source_count":sources.sources().len(),
        "source_ids":sources.sources().iter().map(|s|s.id.clone()).collect::<Vec<_>>(),"zero_domains":zeros.iter().map(bounds).collect::<Vec<_>>(),"admitted_domain":bounds(&admitted),
        "proof_binding_sha256":proof_binding,"physical_arity":prepared.deformation.physical_arity(),"factor_basis":measure.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"coordinates":family.coordinates().iter().map(ToString::to_string).collect::<Vec<_>>(),
        "inverse_routing":family.inverse_routing().iter().map(|row|row.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "roots":roots.iter().map(|p|p.to_vec()).collect::<Vec<_>>(),"all_initial_derivatives":derivatives.iter().map(|p|p.to_vec()).collect::<Vec<_>>(),"selected":selection.iter().map(|p|p.to_vec()).collect::<Vec<_>>(),
        "actual_mask_derivatives":roots.iter().map(|p|serde_json::json!({"root":p.to_vec(),"derivative":prepared.deformation.derivative(*p).unwrap().iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>(),
        "limits":{"depth":3,"ray_domains":1,"point_domains":1,"max_selected_points":16},"probes":probes});
    std::fs::write(
        out.join("result.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!(
        "mode={mode} sources={} points={} fresh replay PASS",
        sources.sources().len(),
        selection.len()
    );
}
