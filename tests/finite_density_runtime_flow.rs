//! Runtime graph predictions through the native public flow owners. No reference
//! values are loaded. A failed or timed-out sector never counts as an amplitude.
use serde_json::{Value, json};
use std::{path::PathBuf, sync::Arc};
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput,
    assembly::PreparedDensityFlow,
    flow::{OccupiedFlowEvaluation, PreparedOccupiedFlow},
    preparation::WeightedSourceOptions,
    reduction::{GuardRefinementOptions, WeightedClosureOptions},
};
use symbolica_amflow::*;

fn setting<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name)
        .map(|v| v.parse().unwrap_or_else(|_| panic!("invalid {name}={v}")))
        .unwrap_or(default)
}

struct Run {
    input: DensityInput,
    report: PathBuf,
    options: FlowOptions,
    closure: WeightedClosureOptions,
    context: RunContext,
}

impl Run {
    fn new() -> Self {
        let input_path = std::env::var("RUSTFLOW_WEIGHTED_INPUT").unwrap_or_else(|_| {
            "examples/finite_density/chain_of_three_parallel_pairs.json".into()
        });
        let input: DensityInput =
            serde_json::from_slice(&std::fs::read(&input_path).unwrap()).unwrap();
        let report = std::env::var_os("RUSTFLOW_DENSITY_FLOW_REPORT")
            .map(PathBuf::from)
            .expect("set a distinct RUSTFLOW_DENSITY_FLOW_REPORT for each bounded run");
        std::fs::create_dir_all(&report).unwrap();
        let options = FlowOptions {
            digits: input.digits,
            guard_digits: 40,
            series_order: 60,
            mass_mode: MassMode::All,
            ..Default::default()
        };
        let closure = WeightedClosureOptions {
            max_rounds: setting("RUSTFLOW_WEIGHTED_ROUNDS", 12),
            max_frontier: setting("RUSTFLOW_WEIGHTED_FRONTIER", 256),
            max_requested: setting("RUSTFLOW_WEIGHTED_REQUESTED", 4096),
            discovery: symbolica_amflow::finite_density::guarded::GuardedDiscoveryOptions {
                max_depth: setting("RUSTFLOW_WEIGHTED_DEPTH", 3),
                max_domains: setting("RUSTFLOW_WEIGHTED_DOMAINS", 8192),
                sample_seed: 0,
            },
            guard_refinement: GuardRefinementOptions {
                max_passes: setting("RUSTFLOW_WEIGHTED_GUARD_PASSES", 3),
                max_added_domains: setting("RUSTFLOW_WEIGHTED_GUARD_DOMAINS", 256),
                max_interval_width: setting("RUSTFLOW_WEIGHTED_GUARD_WIDTH", 2),
            },
            split_ordinary_zero_faces: setting("RUSTFLOW_WEIGHTED_ZERO_FACES", false),
            checkpoints: Some(report.join("native-closure")),
            ..Default::default()
        };
        let run = Self {
            input,
            report,
            options,
            closure,
            context: RunContext {
                progress: Some(Arc::new(|event| {
                    if !matches!(event, Progress::Step { .. }) {
                        eprintln!("{event:?}");
                    }
                })),
                ..Default::default()
            },
        };
        run.save("input.json", json!(run.input));
        run.save("configuration.json", json!({
            "input_path":input_path,
            "native_closure":format!("{:?}",run.closure),
            "guard_digits":run.options.guard_digits,
            "profiles":run.profiles(),
            "source_options":run.source_options(),
            "initial_epsilon_grid":setting("RUSTFLOW_DENSITY_FLOW_GRID",1000_i64),
            "epsilon":std::env::var("RUSTFLOW_DENSITY_FLOW_EPSILON").unwrap_or_else(|_|"1/8".into()),
            "independent_reference_comparisons":0,
        }));
        run
    }

    fn save(&self, file: &str, data: Value) {
        std::fs::write(
            self.report.join(file),
            serde_json::to_vec_pretty(&data).unwrap(),
        )
        .unwrap();
    }

    fn checked<T>(&self, stage: &str, result: Result<T>) -> T {
        result.unwrap_or_else(|error| {
            self.save(
                "failure.json",
                json!({
                    "stage":stage,"error":error.to_string(),"numerical_acceptance":false,
                }),
            );
            panic!("{stage}: {error}")
        })
    }

    fn epsilon(&self) -> Rational {
        symbolica_amflow::finite_density::interface::parse_epsilon(
            &std::env::var("RUSTFLOW_DENSITY_FLOW_EPSILON").unwrap_or_else(|_| "1/8".into()),
        )
        .unwrap()
    }

    fn profiles(&self) -> Vec<(u32, usize, u32)> {
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
            .unwrap_or_else(|_| vec![(18, 60, 8)])
    }

    fn source_options(&self) -> WeightedSourceOptions {
        WeightedSourceOptions {
            policy: std::env::var("RUSTFLOW_WEIGHTED_SOURCE_POLICY")
                .map(|value| value.parse().unwrap())
                .unwrap_or_default(),
            free_virtual_zero_sectors: setting("RUSTFLOW_WEIGHTED_FREE_VIRTUAL_ZEROS", false),
            ..Default::default()
        }
    }

    fn full(&self) -> PreparedDensityFlow {
        self.checked(
            "complete preparation",
            PreparedDensityFlow::prepare_with_source_options(
                &self.input,
                &self.options,
                self.closure.clone(),
                &self.context,
                self.source_options(),
            ),
        )
    }
}

fn report(value: &OccupiedFlowEvaluation) -> Value {
    json!({
        "construction":value.construction,"values":value.values.iter().map(ToString::to_string).collect::<Vec<_>>(),
        "basis_size":value.basis_size,"physical_arity":value.physical_arity,
        "native_storage_capacity":value.native_storage_capacity,"source_options":value.source_options,
        "contour_admission":value.contour_admission,"shifted_slots":value.shifted_slots,
        "massless_endpoint":value.massless_endpoint,"boundary":format!("{:?}",value.boundary),
        "nonzero_conditions":value.nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),
    })
}

fn stable(
    previous: &mut Option<Vec<ComplexFloat>>,
    values: Vec<ComplexFloat>,
    options: &FlowOptions,
) {
    if let Some(old) = previous.as_ref() {
        assert_eq!(old.len(), values.len());
        let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
        for (a, b) in old.iter().zip(&values) {
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
                "native refinement changed {a} to {b}"
            );
        }
    }
    *previous = Some(values);
}

fn occupied<const N: usize>(run: &Run, cuts: &[usize]) {
    let input = run.input.prepare().unwrap();
    let source_options = run.source_options();
    run.save("sector.json",json!({"cut_slots":cuts,"storage_capacity":N,"source_options":source_options,"full_amplitude":false}));
    let flow = run.checked(
        "occupied preparation",
        PreparedOccupiedFlow::<N>::prepare_with_source_options(
            &input,
            cuts,
            &run.options,
            run.closure.clone(),
            &run.context,
            source_options,
        ),
    );
    run.save("closure.json",json!({
        "diagnostics":format!("{:?}",flow.closure_diagnostics()),
        "basis":flow.reduced().basis.iter().map(|i|&i.0).collect::<Vec<_>>(),
        "nonzero_conditions":flow.reduced().nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),
    }));
    let mut previous = None;
    for (digits, order, start) in run.profiles() {
        let options = FlowOptions {
            digits,
            series_order: order,
            ..run.options.clone()
        };
        let result = run.checked(
            "occupied evaluation",
            flow.evaluate_report(&run.epsilon(), &options, &run.context, start),
        );
        run.save(&format!("prediction-{digits}-{order}-{start}.json"),json!({
            "full_amplitude":false,"epsilon":run.epsilon().to_string(),
            "normalization":"unscaled Euclidean cut contribution", "cut_slots":cuts,
            "digits":digits,"guard_digits":options.guard_digits,"series_order":order,"occupied_start_scale":start,
            "report":report(&result),"independent_reference_comparisons":0,
        }));
        stable(&mut previous, result.values, &options);
    }
}

#[test]
#[ignore = "bounded runtime graph occupied AMF; explicit partial sector predictions only"]
fn runtime_graph_occupied_flow() {
    let run = Run::new();
    let cuts = std::env::var("RUSTFLOW_WEIGHTED_CUTS")
        .unwrap_or_else(|_| "0".into())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect::<Vec<usize>>();
    let arity = run
        .input
        .prepare()
        .unwrap()
        .occupied_cut(&cuts, 65536)
        .unwrap()
        .factors()
        .len();
    match arity {
        1..=16 => occupied::<16>(&run, &cuts),
        17..=20 => occupied::<20>(&run, &cuts),
        21..=24 => occupied::<24>(&run, &cuts),
        _ => panic!("runtime diagnostic supports physical arity up to24"),
    }
}

#[test]
#[ignore = "bounded complete runtime graph native amplitude with saved predictions"]
fn runtime_graph_full_amplitude() {
    let run = Run::new();
    let flow = run.full();
    let mut previous = None;
    for (digits, order, start) in run.profiles() {
        let options = FlowOptions {
            digits,
            series_order: order,
            ..run.options.clone()
        };
        let result = run.checked(
            "complete evaluation",
            flow.evaluate(&run.epsilon(), &options, &run.context, start),
        );
        run.save(&format!("prediction-{digits}-{order}-{start}.json"),json!({
            "full_amplitude":true,"epsilon":run.epsilon().to_string(),"normalization":"unscaled Euclidean amplitude",
            "digits":digits,"guard_digits":options.guard_digits,"series_order":order,"occupied_start_scale":start,
            "values":result.values.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "contributions":result.contributions.iter().map(|(cuts,values)|json!({
                "cut_slots":cuts,"values":values.iter().map(ToString::to_string).collect::<Vec<_>>()
            })).collect::<Vec<_>>(),
            "occupied_reports":result.occupied_reports.iter().map(|(cuts,r)|json!({"cut_slots":cuts,"report":report(r)})).collect::<Vec<_>>(),
            "vacuum_zero_certificates":flow.vacuum_zero_certificates().iter().map(|c|c.report()).collect::<Vec<_>>(),
            "empty_support":result.empty_support,"independent_reference_comparisons":0,
        }));
        let values = result
            .contributions
            .into_iter()
            .flat_map(|(_, v)| v)
            .chain(result.values)
            .collect();
        stable(&mut previous, values, &options);
    }
}

#[test]
#[ignore = "bounded complete runtime graph Laurent reconstruction, no reference constraints"]
fn runtime_graph_full_laurent() {
    let run = Run::new();
    let flow = run.full();
    let grid = setting("RUSTFLOW_DENSITY_FLOW_GRID", 1000_i64);
    let mut profiles = run
        .profiles()
        .into_iter()
        .map(|(digits, order, start)| (digits, order, start, grid))
        .collect::<Vec<_>>();
    let &(digits, order, start, _) = profiles.last().expect("at least one numerical profile");
    profiles.push((
        digits,
        order,
        start,
        grid.checked_mul(2).expect("grid overflow"),
    ));
    let mut previous = None;
    for (digits, order, start, grid) in profiles {
        let options = FlowOptions {
            digits,
            series_order: order,
            ..run.options.clone()
        };
        let result = run.checked(
            "complete Laurent reconstruction",
            flow.solve_with_sampling_grid(&options, &run.context, start, grid),
        );
        run.save(&format!("prediction-{digits}-{order}-{start}-grid-{grid}.json"),json!({
            "full_amplitude":true,"normalization":"unscaled Euclidean amplitude",
            "digits":digits,"guard_digits":options.guard_digits,"series_order":order,"occupied_start_scale":start,"epsilon_grid_denominator":grid,
            "expansions":result.iter().map(|e|json!({"verified_digits":e.verified_digits,
                "coefficients":e.coefficients.iter().map(|(order,value)|(order.to_string(),value.to_string())).collect::<std::collections::BTreeMap<_,_>>()
            })).collect::<Vec<_>>(),"independent_reference_comparisons":0,
        }));
        assert!(
            result
                .iter()
                .all(|e| e.verified_digits.is_some_and(|d| d >= digits))
        );
        stable(
            &mut previous,
            result
                .iter()
                .flat_map(|e| e.coefficients.values().cloned())
                .collect(),
            &options,
        );
    }
}
