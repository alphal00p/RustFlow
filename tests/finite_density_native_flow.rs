//! Explicit expensive full-amplitude validation; no oracle values are loaded.
use std::sync::Arc;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput, flow::PreparedOccupiedFlow, guarded::GuardedDiscoveryOptions,
    reduction::WeightedClosureOptions, vacuum::PreparedDensityVacuum,
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
    let closure = |name: &str| WeightedClosureOptions {
        max_rounds: budget("RUSTFLOW_WEIGHTED_ROUNDS", 12),
        discovery: GuardedDiscoveryOptions {
            max_depth: budget("RUSTFLOW_WEIGHTED_DEPTH", 3).try_into().unwrap(),
            max_domains: budget("RUSTFLOW_WEIGHTED_DOMAINS", 8192),
            sample_seed: 0,
        },
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
    let both = PreparedOccupiedFlow::<9>::prepare(
        &prepared,
        &[0, 1],
        &options,
        closure("cut-01"),
        &context,
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
