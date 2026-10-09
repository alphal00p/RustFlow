//! Independent exact verification of optional externally proposed virtual IBP
//! relations. No AMF evaluation, compact boundary, oracle, or numerical answer
//! enters this test. Nonzero residuals mean native replay is inconclusive.
use std::collections::BTreeSet;
use std::path::PathBuf;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::family::LinearCombination;
use symbolica_amflow::finite_density::InversePropagatorBasis;
use symbolica_amflow::{
    Integral, IntegralFamily, Propagator, ReductionBackend, RunContext, RustRedBackend,
};

fn p(text: &str) -> Atom {
    Atom::parse(text, "prism_reference", Default::default()).unwrap()
}
fn candidate_path() -> PathBuf {
    std::env::var_os("RUSTFLOW_PRISM_VIRTUAL_CANDIDATE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("fixtures/finite_density/prism_virtual_candidate_relations.json")
        })
}
fn candidate() -> serde_json::Value {
    serde_json::from_slice(&std::fs::read(candidate_path()).unwrap()).unwrap()
}
fn indices(value: &serde_json::Value) -> Integral {
    Integral(
        value
            .as_array()
            .unwrap()
            .iter()
            .map(|x| i16::try_from(x.as_i64().unwrap()).unwrap())
            .collect(),
    )
}
fn family() -> IntegralFamily {
    let gram = vec![vec![p("-1"), p("-1/2")], vec![p("-1/2"), p("0")]];
    let routings = [
        ([1, 0], [0, 0]),
        ([0, 1], [0, 0]),
        ([1, 0], [-1, 0]),
        ([0, 1], [-1, 0]),
        ([1, 0], [0, -1]),
        ([1, -1], [0, 0]),
        ([0, 1], [0, -1]),
    ];
    IntegralFamily {
        name: "independent_prism_virtual_ibp_certificate".into(),
        loops: vec!["K".into(), "L".into()],
        external: vec!["p".into(), "q".into()],
        propagators: routings
            .iter()
            .map(|(loops, external)| Propagator::quadratic(loops, external, p("0"), &gram).unwrap())
            .collect(),
        external_gram: gram,
        physical_propagators: 6,
        epsilon: match p("eps").as_view() {
            AtomView::Var(v) => v.get_symbol(),
            _ => unreachable!(),
        },
        dimension: 4,
    }
}
fn originals() -> Vec<LinearCombination> {
    let d = (0..7).map(|i| p(&format!("d{i}"))).collect::<Vec<_>>();
    let z = (0..7).map(|i| p(&format!("z{i}"))).collect::<Vec<_>>();
    // This is solely exact polynomial lowering. Native IntegralFamily separately
    // declares the seventh slot to be a numerator-only completion.
    let basis = InversePropagatorBasis::new(d.clone(), d, z).unwrap();
    let na = p("((d0-d4)^2-d0-d1+d5)/4");
    let na_prime = p("(d0-d4)*(d4-d2-1)/2-(d0+d1-d5)/4");
    let nb = p("(1+(d1-1-d3)*(d0-d4))/4");
    let nb_prime = p("1/2+(d1-1-d3)*(d4-d2-1)/4");
    let nc = p("((d0+d4)^2+(d0-d1+d5)*(d4+d2))/4");
    let c1a = p("d4") * na_prime + &na * p("d4-d2");
    let c1b = p("d4") * nb_prime + &nb * p("d4-d2");
    [
        (na, [1, 1, 1, 1, 1, 1, 0]),
        (c1a, [1, 1, 1, 1, 2, 1, 0]),
        (nb, [1, 1, 1, 1, 1, 1, 0]),
        (c1b, [1, 1, 1, 1, 2, 1, 0]),
        (nc, [1, 1, 1, 1, 2, 1, 0]),
    ]
    .into_iter()
    .map(|(n, powers)| basis.convert(&n, &powers).unwrap())
    .collect()
}
fn exact_clean(terms: &mut LinearCombination) {
    for coefficient in terms.values_mut() {
        *coefficient = coefficient.together().cancel();
    }
    terms.retain(|_, c| !c.is_zero());
}
fn row_json(terms: &LinearCombination) -> Vec<serde_json::Value> {
    terms
        .iter()
        .map(|(i, c)| serde_json::json!({"indices":i.0,"coefficient":c.to_canonical_string()}))
        .collect()
}
fn save(path: &std::path::Path, value: &serde_json::Value) {
    std::fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
}

#[test]
fn original_virtual_targets_and_candidate_metadata_match_exactly() {
    let c = candidate();
    assert_eq!(c["dimension"], "4-2*eps");
    assert_eq!(c["complete_prism_reference"], false);
    assert_eq!(c["native_predictions_read"], 0);
    assert_eq!(c["oracle_records_read"], 0);
    let f = family();
    f.validate().unwrap();
    let master_indices = c["master_indices"]
        .as_array()
        .unwrap()
        .iter()
        .map(indices)
        .collect::<Vec<_>>();
    assert_eq!(
        master_indices,
        vec![
            Integral(vec![0, 1, 0, 1, 1, 1, 0]),
            Integral(vec![0, 1, 1, 0, 0, 1, 0]),
            Integral(vec![1, 1, 1, 1, 0, 0, 0])
        ]
    );
    let expected_labels = ["C0A", "C1A", "C0B", "C1B", "CC"];
    let originals = originals();
    for ((row, original), label) in c["relations"]
        .as_array()
        .unwrap()
        .iter()
        .zip(&originals)
        .zip(expected_labels)
    {
        assert_eq!(row["label"], label);
        assert_eq!(row["candidate_coefficients"].as_array().unwrap().len(), 3);
        let external_original = row["original_terms"]
            .as_array()
            .unwrap()
            .iter()
            .map(|term| {
                (
                    indices(&term["indices"]),
                    p(term["coefficient"].as_str().unwrap()),
                )
            })
            .collect::<LinearCombination>();
        assert_eq!(
            &external_original, original,
            "original numerator/jet changed for {label}"
        );
        for integral in original.keys().chain(&master_indices) {
            f.validate_integral(integral).unwrap();
            assert_eq!(integral.0[6], 0);
        }
    }
}

#[test]
fn optional_prism_candidates_are_native_exact_source_consequences() {
    original_virtual_targets_and_candidate_metadata_match_exactly();
    let started = Instant::now();
    let candidate = candidate();
    let originals = originals();
    let masters = candidate["master_indices"]
        .as_array()
        .unwrap()
        .iter()
        .map(indices)
        .collect::<Vec<_>>();
    let targets = originals
        .iter()
        .flat_map(|t| t.keys().cloned())
        .chain(masters.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let requested_report = std::env::var_os("RUSTFLOW_PRISM_VIRTUAL_REPLAY_REPORT");
    let keep_report = requested_report.is_some();
    let report = requested_report.map(PathBuf::from).unwrap_or_else(|| {
        std::env::temp_dir().join(format!(
            "rustflow-prism-native-replay-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    });
    std::fs::create_dir_all(&report).unwrap();
    let depth = std::env::var("RUSTFLOW_PRISM_VIRTUAL_REPLAY_DEPTH")
        .map(|x| x.parse().unwrap())
        .unwrap_or(3);
    let backend = RustRedBackend {
        max_depth: depth,
        bubble_subloops: false,
        symmetry_rules: true,
        checkpoints: Some(report.join("native-checkpoints")),
        ..Default::default()
    };
    let definition = serde_json::json!({"scope":"validation-only generic-epsilon exact source replay; no AMF or numerical values","candidate_path":candidate_path(),"family":{"loops":["K","L"],"external":["p","q"],"external_gram":[["-1","-1/2"],["-1/2","0"]],"denominators":["K^2","L^2","(K-p)^2","(L-p)^2","(K-q)^2","(K-L)^2","(L-q)^2"],"physical_propagators":6,"dimension":"4-2*eps"},"max_depth":depth,"bubble_subloops":false,"symmetry_rules":true,"targets":targets.iter().map(|i|i.0.clone()).collect::<Vec<_>>(),"originals":originals.iter().map(row_json).collect::<Vec<_>>(),"candidate_raw_denominators":["8*eps^3","8*eps^3*(1+eps)","8*eps^3","8*eps^3*(1+eps)*(1+2*eps)","8*eps^3*(1+eps)*(1+2*eps)"],"candidate_nonzero_conditions":["eps","1+eps","1+2*eps"],"backend_identity":backend.identity()});
    save(&report.join("input.json"), &definition);
    let context = RunContext {
        progress: Some(std::sync::Arc::new(|event| {
            eprintln!("prism native replay: {event:?}")
        })),
        ..Default::default()
    };
    let reduced = match backend.reduce(&family(), &targets, &context) {
        Ok(reduced) => reduced,
        Err(error) => {
            save(
                &report.join("result.json"),
                &serde_json::json!({"status":"native_replay_failed","error":error.to_string(),"elapsed_seconds":started.elapsed().as_secs_f64(),"complete_prism_reference":false,"native_predictions_read":0,"oracle_records_read":0}),
            );
            panic!("native prism replay did not complete: {error}");
        }
    };
    save(
        &report.join("native-reduction.json"),
        &serde_json::json!({"rules":reduced.rules.iter().map(|(i,terms)|serde_json::json!({"indices":i.0,"terms":row_json(terms)})).collect::<Vec<_>>(),"residuals":reduced.residuals.iter().map(|i|i.0.clone()).collect::<Vec<_>>(),"nonzero_conditions":reduced.nonzero_conditions.iter().map(|c|c.to_canonical_string()).collect::<Vec<_>>() }),
    );
    let expanded = reduced.expand_many(&targets).unwrap();
    let mut results = Vec::new();
    let mut all_zero = true;
    for (original, row) in originals
        .iter()
        .zip(candidate["relations"].as_array().unwrap())
    {
        let mut difference = LinearCombination::new();
        for (target, coefficient) in original {
            for (residual, factor) in &expanded[target] {
                *difference.entry(residual.clone()).or_default() += coefficient * factor;
            }
        }
        for (master, coefficient) in masters
            .iter()
            .zip(row["candidate_coefficients"].as_array().unwrap())
        {
            let coefficient = p(coefficient.as_str().unwrap());
            for (residual, factor) in &expanded[master] {
                *difference.entry(residual.clone()).or_default() -= &coefficient * factor;
            }
        }
        exact_clean(&mut difference);
        all_zero &= difference.is_empty();
        results.push(serde_json::json!({"label":row["label"],"exact_zero":difference.is_empty(),"residual_difference":row_json(&difference)}));
    }
    save(
        &report.join("result.json"),
        &serde_json::json!({"status":if all_zero {"all_five_candidate_relations_verified_by_native_exact_source_replay"} else {"native_residuals_inconclusive"},"relations":results,"nonzero_conditions":reduced.nonzero_conditions.iter().map(|c|c.to_canonical_string()).collect::<Vec<_>>(),"candidate_nonzero_conditions":["eps","1+eps","1+2*eps"],"elapsed_seconds":started.elapsed().as_secs_f64(),"complete_prism_reference":false,"common_thermal_continuation_proved":false,"native_predictions_read":0,"oracle_records_read":0}),
    );
    assert!(
        all_zero,
        "bounded native replay left a nonzero common residual; candidate proof remains inconclusive"
    );
    if !keep_report {
        std::fs::remove_dir_all(report).unwrap();
    }
}
