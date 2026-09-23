//! Required two-loop acceptance calculation. Any unsupported stage is an error.
use symbolica_amflow::*;

fn main() -> Result<()> {
    let mut native = RustRedBackend {
        checkpoints: Some(".amflow-cache".into()),
        ..Default::default()
    };
    let mut options = FlowOptions {
        workers: 4,
        cache_directory: Some(".amflow-cache".into()),
        ..Default::default()
    };
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--max-targets" => native.max_targets = integer_argument(&mut arguments, &argument)?,
            "--max-exact-frontier" => {
                native.max_exact_frontier = integer_argument(&mut arguments, &argument)?
            }
            "--depth" => {
                native.max_depth = u32::try_from(integer_argument(&mut arguments, &argument)?)
                    .map_err(|_| Error::InvalidInput("depth exceeds u32".into()))?
            }
            "--workers" => options.workers = integer_argument(&mut arguments, &argument)?,
            "--factorized" => native.factorized = true,
            "--plain" => native.factorized = false,
            "--symbolic-epsilon" => options.sampled_reduction = false,
            "--skip-reduction" => options.skip_reduction = true,
            "--help" => {
                println!(
                    "two_loop_acceptance [--max-targets N] [--max-exact-frontier N] [--depth N] [--workers N] [--symbolic-epsilon] [--skip-reduction] [--factorized | --plain]"
                );
                return Ok(());
            }
            _ => return Err(Error::InvalidInput(format!("unknown option {argument}"))),
        }
    }
    let (family, targets) = benchmarks::paper_two_loop()?;
    let backend = cache::CachedBackend {
        backend: native,
        directory: ".amflow-cache".into(),
    };
    let context = RunContext {
        progress: Some(std::sync::Arc::new(|event| {
            if !matches!(event, Progress::Step { .. }) {
                eprintln!("{event:?}");
            }
        })),
        ..Default::default()
    };
    let result = solve_integrals(
        &family,
        &targets,
        &KinematicPoint::default(),
        0,
        &options,
        &backend,
        &context,
    )?;
    benchmarks::validate_paper_two_loop(&targets, &result)?;
    for (target, value) in targets.iter().zip(result) {
        println!("{:?}: {:?}", target.0, value);
    }
    Ok(())
}

fn integer_argument(arguments: &mut impl Iterator<Item = String>, option: &str) -> Result<usize> {
    arguments
        .next()
        .ok_or_else(|| Error::InvalidInput(format!("missing value for {option}")))?
        .parse()
        .map_err(|e| Error::InvalidInput(format!("{option}: {e}")))
}
