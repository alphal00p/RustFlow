//! Inspect native reduction closure for one exact epsilon sample.
use symbolica::prelude::*;
use symbolica_amflow::*;
fn main() -> Result<()> {
    let depth = std::env::args()
        .nth(1)
        .map(|s| s.parse::<u32>())
        .transpose()
        .map_err(|e| Error::InvalidInput(e.to_string()))?
        .unwrap_or(1);
    let max_targets = std::env::args()
        .nth(2)
        .map(|s| s.parse::<usize>())
        .transpose()
        .map_err(|e| Error::InvalidInput(e.to_string()))?
        .unwrap_or(4096);
    let (family, targets) = benchmarks::paper_two_loop()?;
    let backend = cache::CachedBackend {
        backend: RustRedBackend {
            max_depth: depth,
            max_targets,
            factorized: std::env::args().any(|a| a == "factorized"),
            checkpoints: Some(".amflow-cache".into()),
            ..Default::default()
        },
        directory: ".amflow-cache".into(),
    };
    let context = RunContext {
        progress: Some(std::sync::Arc::new(|event| eprintln!("{event:?}"))),
        ..Default::default()
    };
    let options = FlowOptions {
        cache_directory: Some(".amflow-cache".into()),
        ..Default::default()
    };
    let flow = if std::env::args().nth(3).as_deref() == Some("symbolic") {
        PreparedFlow::new(
            &family,
            &targets,
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )?
    } else {
        PreparedFlow::new_at_epsilon(
            &family,
            &targets,
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
            &Rational::from((1, 2700)),
        )?
    };
    println!(
        "closed basis: {} integrals; blocks: {:?}",
        flow.reduced.basis.len(),
        flow.system
            .blocks()?
            .iter()
            .map(Vec::len)
            .collect::<Vec<_>>()
    );
    Ok(())
}
