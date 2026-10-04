//! Check a subsector of the physical two-loop benchmark at two precisions.
//! Optional argument: nine comma-separated signed propagator powers.
use symbolica::prelude::*;
use symbolica_amflow::*;

fn main() -> Result<()> {
    let args = std::env::args().collect::<Vec<_>>();
    let depth = args
        .iter()
        .position(|s| s == "--depth")
        .map(|i| {
            args.get(i + 1)
                .ok_or_else(|| Error::InvalidInput("missing depth".into()))?
                .parse::<u32>()
                .map_err(|e| Error::InvalidInput(e.to_string()))
        })
        .transpose()?
        .unwrap_or(2);
    let max_exact_frontier = args
        .iter()
        .position(|s| s == "--max-exact-frontier")
        .map(|i| {
            args.get(i + 1)
                .ok_or_else(|| Error::InvalidInput("missing frontier limit".into()))?
                .parse::<usize>()
                .map_err(|e| Error::InvalidInput(e.to_string()))
        })
        .transpose()?
        .unwrap_or(512);
    let max_targets = args
        .iter()
        .position(|s| s == "--max-targets")
        .map(|i| {
            args.get(i + 1)
                .ok_or_else(|| Error::InvalidInput("missing target limit".into()))?
                .parse::<usize>()
                .map_err(|e| Error::InvalidInput(e.to_string()))
        })
        .transpose()?
        .unwrap_or(32768);
    let max_backward_frontier = args
        .iter()
        .position(|s| s == "--max-backward-frontier")
        .map(|i| {
            args.get(i + 1)
                .ok_or_else(|| Error::InvalidInput("missing backward frontier limit".into()))?
                .parse::<usize>()
                .map_err(|e| Error::InvalidInput(e.to_string()))
        })
        .transpose()?
        .unwrap_or(RustRedBackend::default().max_backward_frontier);
    let powers = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "0,1,1,0,1,0,1,0,0".into())
        .split(',')
        .map(|v| {
            v.parse::<i16>()
                .map_err(|e| Error::InvalidInput(e.to_string()))
        })
        .collect::<Result<Vec<_>>>()?;
    let (family, _) = benchmarks::paper_two_loop()?;
    let target = Integral(powers);
    family.validate_integral(&target)?;
    let backend = cache::CachedBackend {
        backend: RustRedBackend {
            max_depth: depth,
            parametric_rules: args.iter().any(|s| s == "--parametric-rules"),
            symmetry_rules: args.iter().any(|s| s == "--symmetry-rules"),
            max_targets,
            max_sector_batch: 32,
            max_exact_frontier,
            max_backward_frontier,
            checkpoints: Some(".amflow-cache/subsectors".into()),
            ..Default::default()
        },
        directory: ".amflow-cache/subsectors".into(),
    };
    let context = RunContext {
        progress: Some(std::sync::Arc::new(|event| {
            if !matches!(event, Progress::Step { .. }) {
                eprintln!("{event:?}");
            }
        })),
        ..Default::default()
    };
    let options = FlowOptions {
        cache_directory: Some(".amflow-cache/subsectors".into()),
        ..Default::default()
    };
    let epsilon = Rational::from((1, 2700));
    let prepared = PreparedFlow::new_at_epsilon(
        &family,
        std::slice::from_ref(&target),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
        &epsilon,
    )?;
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    if std::env::args().any(|s| s == "--indicial") {
        let p = Precision::decimal(60)?;
        let values = ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(&epsilon))]);
        for (name, system) in [
            ("zero", prepared.system.clone()),
            (
                "infinity",
                prepared.system.invert_variable(symbol!("subsector_z")),
            ),
        ] {
            let basis = system.frobenius(p, &values, 8)?;
            eprintln!(
                "{name} indicial roots: {:?}",
                basis
                    .columns
                    .iter()
                    .map(|c| &c.exponent)
                    .collect::<Vec<_>>()
            );
        }
    }
    let first = prepared.evaluate(&epsilon, &options, &boundary, &context)?;
    let refined = FlowOptions {
        guard_digits: options.guard_digits + 20,
        series_order: options.series_order + 32,
        ..options.clone()
    };
    let second = prepared.evaluate(&epsilon, &refined, &boundary, &context)?;
    let p = Precision::decimal(refined.digits + refined.guard_digits)?;
    if !p.close(&first[0], &second[0], options.digits) {
        return Err(Error::Numerical(
            "subsector precision refinement failed".into(),
        ));
    }
    println!("{:?} at epsilon={epsilon}: {}", target.0, second[0]);
    println!(
        "stable to {} digits at increased precision/order",
        options.digits
    );
    Ok(())
}
