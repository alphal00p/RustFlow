//! Explicit expensive full-amplitude validation; no oracle values are loaded.
use std::sync::Arc;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput,
    flow::PreparedOccupiedFlow,
    guarded::GuardedDiscoveryOptions,
    preparation::WeightedSourceOptions,
    reduction::{GuardRefinementOptions, WeightedClosureOptions},
    vacuum::PreparedDensityVacuum,
};
use symbolica_amflow::*;

#[test]
#[ignore = "explicit vacuum and single-cut transport, independent of double-cut closure"]
fn massive_vacuum_and_single_cut_transport() {
    let input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    let prepared = input.prepare().unwrap();
    let context = RunContext {
        progress: Some(Arc::new(|event| eprintln!("{event:?}"))),
        ..Default::default()
    };
    let options = FlowOptions {
        digits: 18,
        guard_digits: 20,
        series_order: 60,
        mass_mode: MassMode::All,
        ..Default::default()
    };
    let report = std::env::var_os("RUSTFLOW_DENSITY_FLOW_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-density-single-sunset"));
    std::fs::create_dir_all(&report).unwrap();
    std::fs::write(
        report.join("input.json"),
        serde_json::to_vec_pretty(&input).unwrap(),
    )
    .unwrap();
    let vacuum = PreparedDensityVacuum::prepare(&prepared, &options, &context).unwrap();
    let single = PreparedOccupiedFlow::<7>::prepare(
        &prepared,
        &[0],
        &options,
        WeightedClosureOptions {
            max_rounds: 12,
            discovery: GuardedDiscoveryOptions {
                max_depth: 3,
                max_domains: 8192,
                sample_seed: 0,
            },
            checkpoints: Some(report.join("cut-0")),
            ..Default::default()
        },
        &context,
    )
    .unwrap();
    let epsilon = Rational::from((4, 5));
    for (digits, order, start) in [(18, 60, 8), (28, 60, 8), (28, 80, 8), (28, 80, 12)] {
        let refined = FlowOptions {
            digits,
            series_order: order,
            ..options.clone()
        };
        for sector in ["vacuum", "cut-0"] {
            let values = if sector == "vacuum" {
                vacuum.evaluate(&epsilon, &refined, &context)
            } else {
                single.evaluate_with_start_scale(&epsilon, &refined, &context, start)
            }
            .unwrap();
            std::fs::write(report.join(format!("prediction-{sector}-{digits}-{order}-{start}.json")),
                serde_json::to_vec_pretty(&serde_json::json!({
                    "sector":sector, "epsilon":"4/5", "normalization":"unscaled Euclidean amplitude",
                    "digits":digits,"series_order":order,"occupied_start_scale":start,
                    "values":values.iter().map(ToString::to_string).collect::<Vec<_>>(),
                    "independent_reference_comparisons":0,"full_amplitude":false
                })).unwrap()).unwrap();
        }
    }
}

#[test]
#[ignore = "explicit native reduction and full massive amplitude refinement"]
fn complete_massive_sunset_assembles_vacuum_and_all_occupied_sectors() {
    let input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    let prepared = input.prepare().unwrap();
    let context = RunContext {
        progress: Some(Arc::new(|event| eprintln!("{event:?}"))),
        ..Default::default()
    };
    let options = FlowOptions {
        digits: 18,
        guard_digits: 20,
        series_order: 80,
        mass_mode: MassMode::All,
        ..Default::default()
    };
    let report = std::env::var_os("RUSTFLOW_DENSITY_FLOW_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-density-full-sunset"));
    std::fs::create_dir_all(&report).unwrap();
    std::fs::write(
        report.join("input.json"),
        serde_json::to_vec_pretty(&input).unwrap(),
    )
    .unwrap();
    let budget = |name: &str, default: usize| {
        std::env::var(name)
            .ok()
            .map(|s| s.parse().unwrap())
            .unwrap_or(default)
    };
    let guard_refinement = GuardRefinementOptions {
        max_passes: budget("RUSTFLOW_WEIGHTED_GUARD_PASSES", 0),
        max_added_domains: budget("RUSTFLOW_WEIGHTED_GUARD_DOMAINS", 256),
        max_interval_width: budget("RUSTFLOW_WEIGHTED_GUARD_WIDTH", 2)
            .try_into()
            .unwrap(),
    };
    let closure = |name: &str| WeightedClosureOptions {
        max_rounds: budget("RUSTFLOW_WEIGHTED_ROUNDS", 12),
        discovery: GuardedDiscoveryOptions {
            max_depth: budget("RUSTFLOW_WEIGHTED_DEPTH", 3).try_into().unwrap(),
            max_domains: budget("RUSTFLOW_WEIGHTED_DOMAINS", 8192),
            sample_seed: 0,
        },
        guard_refinement,
        max_frontier: budget("RUSTFLOW_WEIGHTED_FRONTIER", 256),
        max_requested: budget("RUSTFLOW_WEIGHTED_REQUESTED", 4096),
        checkpoints: Some(report.join(name)),
        ..Default::default()
    };
    let vacuum = PreparedDensityVacuum::prepare(&prepared, &options, &context).unwrap();
    let first =
        PreparedOccupiedFlow::<7>::prepare(&prepared, &[0], &options, closure("cut-0"), &context)
            .unwrap();
    let second =
        PreparedOccupiedFlow::<7>::prepare(&prepared, &[1], &options, closure("cut-1"), &context)
            .unwrap();
    let source_options = WeightedSourceOptions {
        policy: std::env::var("RUSTFLOW_WEIGHTED_SOURCE_POLICY")
            .ok()
            .map(|s| s.parse().unwrap())
            .unwrap_or_default(),
        positive_compact_energy_powers: std::env::var("RUSTFLOW_WEIGHTED_POSITIVE_ENERGY_POWERS")
            .is_ok_and(|value| value == "1"),
    };
    let both = PreparedOccupiedFlow::<9>::prepare_with_source_options(
        &prepared,
        &[0, 1],
        &options,
        closure("cut-01"),
        &context,
        source_options,
    )
    .unwrap();
    let epsilon_text =
        std::env::var("RUSTFLOW_DENSITY_FLOW_EPSILON").unwrap_or_else(|_| "1/7".into());
    let epsilon = Rational::try_from(
        Atom::parse(&epsilon_text, "density_flow_test", Default::default())
            .unwrap()
            .as_view(),
    )
    .unwrap();
    let mut previous: Option<Vec<ComplexFloat>> = None;
    for (digits, order, start_scale) in [(18, 60, 8), (28, 60, 8), (28, 80, 8), (28, 80, 12)] {
        let options = FlowOptions {
            digits,
            series_order: order,
            ..options.clone()
        };
        let p = Precision::decimal(digits + options.guard_digits).unwrap();
        let parts = [
            vacuum.evaluate(&epsilon, &options, &context).unwrap(),
            first
                .evaluate_with_start_scale(&epsilon, &options, &context, start_scale)
                .unwrap(),
            second
                .evaluate_with_start_scale(&epsilon, &options, &context, start_scale)
                .unwrap(),
            both.evaluate_with_start_scale(&epsilon, &options, &context, start_scale)
                .unwrap(),
        ];
        assert!(
            p.norm(&parts[0][0]) > p.tolerance(10),
            "massive vacuum was discarded"
        );
        let assembled = (0..input.targets.len())
            .map(|i| parts.iter().fold(p.zero(), |v, part| p.add(&v, &part[i])))
            .collect::<Vec<_>>();
        // Predictions are saved before any independent comparison. This run
        // establishes refinement stability only, not reference agreement.
        std::fs::write(report.join(format!("prediction-{digits}-{order}-{start_scale}.json")),serde_json::to_vec_pretty(&serde_json::json!({
            "normalization":"unscaled Euclidean amplitude", "epsilon":epsilon_text,
            "digits":digits,"series_order":order,"occupied_start_scale":start_scale,
            "guard_refinement":guard_refinement,
            "double_cut_closure_diagnostics":format!("{:?}",both.closure_diagnostics()),
            "double_cut_source_policy":source_options.policy.as_str(),
            "double_cut_positive_compact_energy_powers":source_options.positive_compact_energy_powers,
            "stability_criterion":"relative 1e-12 for resolved nonzero values; absolute 1e-25 Euclidean units for values below 1e-20",
            "contributions":parts.iter().map(|v|v.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
            "assembled":assembled.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "independent_reference_comparisons":0
        })).unwrap()).unwrap();
        if let Some(previous) = previous {
            for (a, b) in previous.iter().zip(&assembled) {
                let an = p.norm(a);
                let bn = p.norm(b);
                let scale = if an > bn { an } else { bn };
                let tolerance = if scale < p.tolerance(20) {
                    p.tolerance(25)
                } else {
                    p.tolerance(12) * scale
                };
                assert!(
                    p.norm(&p.sub(a, b)) <= tolerance,
                    "full massive amplitude failed independent precision/order/start refinement"
                );
            }
        }
        previous = Some(assembled);
    }
}

#[test]
#[ignore = "explicit complete massive native AMF Laurent refinement"]
fn complete_massive_sunset_laurent_refinement() {
    use symbolica_amflow::finite_density::{
        assembly::PreparedDensityFlow, reduction::GuardRefinementOptions,
    };
    let input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    let context = RunContext {
        progress: Some(Arc::new(|event| eprintln!("{event:?}"))),
        ..Default::default()
    };
    let options = FlowOptions {
        digits: 18,
        guard_digits: 24,
        series_order: 60,
        mass_mode: MassMode::All,
        ..Default::default()
    };
    let report = std::env::var_os("RUSTFLOW_DENSITY_FLOW_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-density-full-sunset-laurent"));
    std::fs::create_dir_all(&report).unwrap();
    std::fs::write(
        report.join("input.json"),
        serde_json::to_vec_pretty(&input).unwrap(),
    )
    .unwrap();
    let budget = |name: &str, default: usize| {
        std::env::var(name)
            .ok()
            .map(|value| value.parse().unwrap())
            .unwrap_or(default)
    };
    let closure = WeightedClosureOptions {
        max_rounds: budget("RUSTFLOW_WEIGHTED_ROUNDS", 12),
        discovery: GuardedDiscoveryOptions {
            max_depth: budget("RUSTFLOW_WEIGHTED_DEPTH", 3).try_into().unwrap(),
            max_domains: budget("RUSTFLOW_WEIGHTED_DOMAINS", 8192),
            sample_seed: 0,
        },
        guard_refinement: GuardRefinementOptions {
            max_passes: budget("RUSTFLOW_WEIGHTED_GUARD_PASSES", 3),
            max_added_domains: budget("RUSTFLOW_WEIGHTED_GUARD_DOMAINS", 256),
            max_interval_width: budget("RUSTFLOW_WEIGHTED_GUARD_WIDTH", 2)
                .try_into()
                .unwrap(),
        },
        checkpoints: Some(report.join("native-closure")),
        ..Default::default()
    };
    let guard_refinement = closure.guard_refinement;
    let flow = PreparedDensityFlow::prepare(&input, &options, closure, &context).unwrap();
    let mut previous: Option<Vec<symbolica_amflow::LaurentExpansion>> = None;
    for (digits, order, start) in [(18, 60, 8), (28, 60, 8), (28, 80, 8), (28, 80, 12)] {
        let refined = FlowOptions {
            digits,
            series_order: order,
            ..options.clone()
        };
        let values = flow
            .solve_with_start_scale(&refined, &context, start)
            .unwrap();
        // Save predictions before any independent-reference comparison.
        std::fs::write(report.join(format!("prediction-{digits}-{order}-{start}.json")),
            serde_json::to_vec_pretty(&serde_json::json!({
                "full_amplitude":true,"normalization":"unscaled Euclidean amplitude",
                "digits":digits,"series_order":order,"occupied_start_scale":start,
                "independent_reference_comparisons":0, "guard_refinement":guard_refinement,
                "expansions":values.iter().map(|value|serde_json::json!({
                    "verified_digits":value.verified_digits,"working_bits":value.working_bits,
                    "samples":value.samples,"validation_samples":value.validation_samples,
                    "refinements":value.refinements,
                    "coefficients":value.coefficients.iter().map(|(power,value)|
                        (power.to_string(),value.to_string())).collect::<std::collections::BTreeMap<_,_>>()
                })).collect::<Vec<_>>()
            })).unwrap()).unwrap();
        let p = Precision::decimal(digits + options.guard_digits).unwrap();
        assert!(values.iter().all(|value| {
            value
                .verified_digits
                .is_some_and(|verified| verified >= digits)
        }));
        if let Some(previous) = previous {
            for (a, b) in previous.iter().zip(&values) {
                for (power, coefficient) in &b.coefficients {
                    let scale = p.norm(coefficient);
                    let tolerance = if scale < p.tolerance(18) {
                        p.tolerance(20)
                    } else {
                        scale * p.tolerance(12)
                    };
                    assert!(
                        p.norm(&p.sub(&a.coefficients[power], coefficient)) <= tolerance,
                        "full massive Laurent coefficient {power} failed independent refinement"
                    );
                }
            }
        }
        previous = Some(values);
    }
}
