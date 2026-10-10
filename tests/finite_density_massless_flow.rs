//! Native predictions only. Independent references are generated and compared
//! separately, after these files have been saved.
use std::{path::PathBuf, sync::Arc};
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput, assembly::PreparedDensityFlow, reduction::WeightedClosureOptions,
};
use symbolica_amflow::*;

fn profiles() -> Vec<(u32, usize, u32)> {
    std::env::var("RUSTFLOW_DENSITY_FLOW_PROFILES")
        .map(|text| {
            text.split(',')
                .map(|entry| {
                    let fields = entry.split(':').collect::<Vec<_>>();
                    assert_eq!(fields.len(), 3);
                    (
                        fields[0].parse().unwrap(),
                        fields[1].parse().unwrap(),
                        fields[2].parse().unwrap(),
                    )
                })
                .collect()
        })
        .unwrap_or_else(|_| vec![(18, 60, 8), (28, 60, 8), (28, 80, 8), (28, 80, 12)])
}

fn prepare(suffix: &str) -> (PreparedDensityFlow, FlowOptions, RunContext, PathBuf) {
    let input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massless_two_loop_sunset.json"
    ))
    .unwrap();
    let report = std::env::var_os("RUSTFLOW_DENSITY_FLOW_REPORT")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join(format!("rustflow-massless-sunset-{suffix}")));
    std::fs::create_dir_all(&report).unwrap();
    std::fs::write(
        report.join("input.json"),
        serde_json::to_vec_pretty(&input).unwrap(),
    )
    .unwrap();
    let options = FlowOptions {
        digits: 18,
        guard_digits: 40,
        series_order: 60,
        mass_mode: MassMode::All,
        ..Default::default()
    };
    let context = RunContext {
        progress: Some(Arc::new(|event| {
            if !matches!(event, Progress::Step { .. }) {
                eprintln!("{event:?}");
            }
        })),
        ..Default::default()
    };
    let flow = PreparedDensityFlow::prepare(
        &input,
        &options,
        WeightedClosureOptions {
            checkpoints: Some(report.join("native-closure")),
            ..Default::default()
        },
        &context,
    )
    .unwrap();
    assert_eq!(
        flow.closure_diagnostics().len(),
        1,
        "the nonzero double cut must use native AMF"
    );
    assert_eq!(flow.physical_zero_certificates().len(), 2);
    std::fs::write(report.join("physical-zero-certificates.json"),
        serde_json::to_vec_pretty(&serde_json::json!(flow.physical_zero_certificates()
            .iter().map(|(cuts, proof)| serde_json::json!({"cut_slots":cuts,"certificate":proof.report()}))
            .collect::<Vec<_>>())).unwrap()).unwrap();
    std::fs::write(
        report.join("closure.json"),
        serde_json::to_vec_pretty(&serde_json::json!(
            flow.closure_diagnostics()
                .iter()
                .map(|(cuts, details)| serde_json::json!({
                    "cut_slots":cuts, "diagnostics":format!("{details:?}")
                }))
                .collect::<Vec<_>>()
        ))
        .unwrap(),
    )
    .unwrap();
    (flow, options, context, report)
}

fn stable(previous: &[ComplexFloat], current: &[ComplexFloat], p: Precision) {
    assert_eq!(previous.len(), current.len());
    for (a, b) in previous.iter().zip(current) {
        let an = p.norm(a);
        let bn = p.norm(b);
        let scale = if an > bn { an } else { bn };
        let tolerance = if scale < p.tolerance(20) {
            p.tolerance(25)
        } else {
            scale * p.tolerance(12)
        };
        assert!(
            p.norm(&p.sub(a, b)) <= tolerance,
            "independent native refinement changed {a} to {b}"
        );
    }
}

#[test]
#[ignore = "bounded complete massless native AMF, with independent numerical refinements"]
fn complete_massless_sunset_fixed_dimension() {
    let (flow, base, context, report) = prepare("fixed");
    let epsilon_text =
        std::env::var("RUSTFLOW_DENSITY_FLOW_EPSILON").unwrap_or_else(|_| "-5/4".into());
    let epsilon = Rational::try_from(
        Atom::parse(&epsilon_text, "massless_flow_test", Default::default())
            .unwrap()
            .as_view(),
    )
    .unwrap();
    let mut previous: Option<Vec<ComplexFloat>> = None;
    for (digits, order, start) in profiles() {
        let options = FlowOptions {
            digits,
            series_order: order,
            ..base.clone()
        };
        let result = flow.evaluate(&epsilon, &options, &context, start).unwrap();
        assert!(result.empty_support.is_empty());
        assert_eq!(
            result
                .contributions
                .iter()
                .map(|(cuts, _)| cuts.clone())
                .collect::<Vec<_>>(),
            vec![vec![], vec![0], vec![1], vec![0, 1]]
        );
        assert!(
            result
                .occupied_reports
                .iter()
                .all(|(cuts, r)| r.construction
                    == if cuts.len() == 1 {
                        "certified_massless_singleton_physical_zero"
                    } else {
                        "weighted_amf"
                    })
        );
        std::fs::write(report.join(format!("prediction-{digits}-{order}-{start}.json")),
            serde_json::to_vec_pretty(&serde_json::json!({
                "full_amplitude":true,"epsilon":epsilon_text,"normalization":"unscaled Euclidean amplitude",
                "digits":digits,"guard_digits":options.guard_digits,"series_order":order,"occupied_start_scale":start,
                "contributions":result.contributions.iter().map(|(cuts,values)| serde_json::json!({
                    "cut_slots":cuts,"values":values.iter().map(ToString::to_string).collect::<Vec<_>>()
                })).collect::<Vec<_>>(),
                "values":result.values.iter().map(ToString::to_string).collect::<Vec<_>>(),
                "occupied_reports":result.occupied_reports.iter().map(|(cuts,r)| serde_json::json!({
                    "cut_slots":cuts,"construction":r.construction,"basis_size":r.basis_size,
                    "physical_arity":r.physical_arity,"native_storage_capacity":r.native_storage_capacity,
                    "contour_admission":r.contour_admission,"boundary":format!("{:?}",r.boundary),
                    "massless_endpoint":r.massless_endpoint,
                    "nonzero_conditions":r.nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>()
                })).collect::<Vec<_>>(),
                "independent_reference_comparisons":0
            })).unwrap()).unwrap();
        let values = result
            .contributions
            .into_iter()
            .flat_map(|(_, v)| v)
            .chain(result.values)
            .collect::<Vec<_>>();
        if let Some(previous) = &previous {
            stable(
                previous,
                &values,
                Precision::decimal(digits + options.guard_digits).unwrap(),
            );
        }
        previous = Some(values);
    }
}

#[test]
#[ignore = "bounded complete massless Laurent fit with independent grid refinement"]
fn complete_massless_sunset_laurent() {
    let (flow, base, context, report) = prepare("laurent");
    let mut profiles = profiles()
        .into_iter()
        .map(|(d, o, s)| (d, o, s, 1000))
        .collect::<Vec<_>>();
    let &(digits, order, start, _) = profiles.last().unwrap();
    profiles.push((digits, order, start, 2000));
    let mut previous: Option<Vec<ComplexFloat>> = None;
    for (digits, order, start, grid) in profiles {
        let options = FlowOptions {
            digits,
            series_order: order,
            ..base.clone()
        };
        let expansions = flow
            .solve_with_sampling_grid(&options, &context, start, grid)
            .unwrap();
        std::fs::write(report.join(format!("prediction-{digits}-{order}-{start}-grid-{grid}.json")),
            serde_json::to_vec_pretty(&serde_json::json!({
                "full_amplitude":true,"normalization":"unscaled Euclidean amplitude",
                "digits":digits,"guard_digits":options.guard_digits,"series_order":order,
                "occupied_start_scale":start,"epsilon_grid_denominator":grid,
                "expansions":expansions.iter().map(|e| serde_json::json!({
                    "verified_digits":e.verified_digits,"coefficients":e.coefficients.iter()
                        .map(|(order,value)| (order.to_string(),value.to_string())).collect::<std::collections::BTreeMap<_,_>>()
                })).collect::<Vec<_>>(),"independent_reference_comparisons":0
            })).unwrap()).unwrap();
        assert!(
            expansions
                .iter()
                .all(|e| e.verified_digits.is_some_and(|verified| verified >= digits))
        );
        let values = expansions
            .iter()
            .flat_map(|e| e.coefficients.values().cloned())
            .collect::<Vec<_>>();
        if let Some(previous) = &previous {
            stable(
                previous,
                &values,
                Precision::decimal(digits + options.guard_digits).unwrap(),
            );
        }
        previous = Some(values);
    }
}
