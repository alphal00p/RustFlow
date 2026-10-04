//! Required two-loop acceptance calculation. Any unsupported stage is an error.
use symbolica::prelude::*;
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
    let mut sample = None;
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--sample" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| Error::InvalidInput("missing epsilon fraction".into()))?;
                let (numerator, denominator) = value.split_once('/').ok_or_else(|| {
                    Error::InvalidInput("--sample expects an exact NUM/DEN fraction".into())
                })?;
                let numerator = numerator
                    .parse::<i64>()
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                let denominator = denominator
                    .parse::<i64>()
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                if numerator == 0 || denominator == 0 {
                    return Err(Error::InvalidInput(
                        "sample numerator and denominator must be nonzero".into(),
                    ));
                }
                sample = Some(Rational::from((numerator, denominator)));
            }
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
            "--max-backward-frontier" => {
                native.max_backward_frontier = integer_argument(&mut arguments, &argument)?
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
            "--symmetry-rules" => native.symmetry_rules = true,
            "--symbolic-epsilon" => options.sampled_reduction = false,
            "--skip-reduction" => options.skip_reduction = true,
            "--help" => {
                println!(
                    "two_loop_acceptance [--sample NUM/DEN] [--max-targets N] [--case-batch N] [--max-exact-frontier N] [--max-backward-frontier N] [--depth N] [--workers N] [--native-workers N] [--checkpoint-seconds N] [--cancel-file PATH] [--symbolic-epsilon] [--skip-reduction] [--factorized | --plain] [--no-bubble-subloops] [--parametric-rules] [--symmetry-rules]"
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
    let result = if let Some(epsilon) = sample {
        evaluate_sample(&family, &targets, &backend, &options, &context, &epsilon).map(|()| None)
    } else {
        solve_integrals(
            &family,
            &targets,
            &KinematicPoint::default(),
            0,
            &options,
            &backend,
            &context,
        )
        .map(Some)
    };
    drop(stop_monitor);
    if let Some(monitor) = monitor {
        monitor
            .join()
            .map_err(|_| std::io::Error::other("cancellation monitor panicked"))?;
    }
    let Some(result) = result? else {
        return Ok(());
    };
    benchmarks::validate_paper_two_loop(&targets, &result)?;
    for (target, value) in targets.iter().zip(result) {
        println!("{:?}: {:?}", target.0, value);
    }
    Ok(())
}

/// A bounded prerequisite for the complete Laurent acceptance calculation.
fn evaluate_sample(
    family: &IntegralFamily,
    targets: &[Integral],
    backend: &dyn ReductionBackend,
    options: &FlowOptions,
    context: &RunContext,
    epsilon: &Rational,
) -> Result<()> {
    let prepared = if options.sampled_reduction {
        PreparedFlow::new_at_epsilon(
            family,
            targets,
            &KinematicPoint::default(),
            backend,
            options,
            context,
            epsilon,
        )?
    } else {
        PreparedFlow::new(
            family,
            targets,
            &KinematicPoint::default(),
            backend,
            options,
            context,
        )?
    };
    let boundary = recursive::RecursiveBoundary::new(backend, options, context);
    let first = prepared.evaluate(epsilon, options, &boundary, context)?;
    let refined = FlowOptions {
        guard_digits: options.guard_digits + 20,
        series_order: options.series_order + 32,
        ..options.clone()
    };
    let second = prepared.evaluate(epsilon, &refined, &boundary, context)?;
    let precision = Precision::decimal(refined.digits + refined.guard_digits)?;
    for ((target, initial), value) in targets.iter().zip(first).zip(second) {
        if !precision.close(&initial, &value, options.digits) {
            return Err(Error::Accuracy(format!(
                "sample precision refinement failed for {:?}",
                target.0
            )));
        }
        println!("{:?} at epsilon={epsilon}: {value}", target.0);
        println!(
            "change under increased precision/order: {}",
            precision.norm(&precision.sub(&initial, &value))
        );
    }
    println!(
        "all four targets at this epsilon stable to {} digits; Laurent acceptance requires the full run",
        options.digits
    );
    Ok(())
}

fn integer_argument(arguments: &mut impl Iterator<Item = String>, option: &str) -> Result<usize> {
    arguments
        .next()
        .ok_or_else(|| Error::InvalidInput(format!("missing value for {option}")))?
        .parse()
        .map_err(|e| Error::InvalidInput(format!("{option}: {e}")))
}
