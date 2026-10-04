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
    let mut cancel_file: Option<std::path::PathBuf> = None;
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--checkpoint-seconds" => {
                native.checkpoint_interval = std::time::Duration::from_secs(
                    u64::try_from(integer_argument(&mut arguments, &argument)?).map_err(|_| {
                        Error::InvalidInput("checkpoint interval exceeds u64".into())
                    })?,
                );
            }
            "--cancel-file" => {
                cancel_file = Some(
                    arguments
                        .next()
                        .ok_or_else(|| Error::InvalidInput("missing cancellation file".into()))?
                        .into(),
                )
            }
            "--max-targets" => native.max_targets = integer_argument(&mut arguments, &argument)?,
            "--native-workers" => {
                native.native_workers = integer_argument(&mut arguments, &argument)?
            }
            "--case-batch" => {
                native.max_sector_batch = integer_argument(&mut arguments, &argument)?
            }
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
            "--no-bubble-subloops" => native.bubble_subloops = false,
            "--parametric-rules" => native.parametric_rules = true,
            "--symbolic-epsilon" => options.sampled_reduction = false,
            "--skip-reduction" => options.skip_reduction = true,
            "--help" => {
                println!(
                    "two_loop_acceptance [--max-targets N] [--case-batch N] [--max-exact-frontier N] [--depth N] [--workers N] [--native-workers N] [--checkpoint-seconds N] [--cancel-file PATH] [--symbolic-epsilon] [--skip-reduction] [--factorized | --plain] [--no-bubble-subloops] [--parametric-rules]"
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
    if cancel_file.as_ref().is_some_and(|p| p.exists()) {
        context.cancellation.cancel();
    }
    let (stop_monitor, stopped) = std::sync::mpsc::channel::<()>();
    let monitor = cancel_file
        .map(|path| {
            let token = context.cancellation.clone();
            std::thread::Builder::new()
                .name("amflow-cancel".into())
                .spawn(move || {
                    loop {
                        if path.exists() {
                            token.cancel();
                            return;
                        }
                        if !matches!(
                            stopped.recv_timeout(std::time::Duration::from_millis(250)),
                            Err(std::sync::mpsc::RecvTimeoutError::Timeout)
                        ) {
                            return;
                        }
                    }
                })
        })
        .transpose()?;
    let result = solve_integrals(
        &family,
        &targets,
        &KinematicPoint::default(),
        0,
        &options,
        &backend,
        &context,
    );
    drop(stop_monitor);
    if let Some(monitor) = monitor {
        monitor
            .join()
            .map_err(|_| std::io::Error::other("cancellation monitor panicked"))?;
    }
    let result = result?;
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
