//! Load test for RustRed harness solutions: read family + solution, build the
//! supplied fixed-d0 eta connection, prepare it with
//! `PreparedFlow::from_supplied_at_epsilon`, compile it numerically and print
//! dimensions, singular points and peak memory. Computes no boundary and no transport.
//!
//! usage: rustred_harness_load FAMILY.json SOLUTION.json [--targets NAMES.json] [--digits N]
//!   NAMES.json: class config with "names": {name: index vector} (default: none -> first master)
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::rustred_harness::{self, HarnessFamily, HarnessSolution};
use symbolica_amflow::*;

fn rss_gb() -> (f64, f64) {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let field = |key: &str| {
        status
            .lines()
            .find(|l| l.starts_with(key))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0)
            / 1e6
    };
    (field("VmRSS:"), field("VmHWM:"))
}

type AnyResult<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

/// Exact symbolic preparation of a several-hundred-master system recurses
/// deeply; run it on a thread with a large stack.
fn main() -> AnyResult<()> {
    std::thread::Builder::new()
        .name("harness-load".into())
        .stack_size(1 << 30)
        .spawn(run)?
        .join()
        .map_err(|_| "load thread panicked")?
}

fn run() -> AnyResult<()> {
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() < 3 {
        eprintln!(
            "usage: rustred_harness_load FAMILY.json SOLUTION.json [--targets CONFIG.json] [--digits N] [--homogeneous-eblock]"
        );
        std::process::exit(2);
    }
    let opt = |key: &str| {
        args.iter()
            .position(|a| a == key)
            .map(|i| args[i + 1].clone())
    };
    let digits: u32 = opt("--digits")
        .map(|d| d.parse())
        .transpose()?
        .unwrap_or(30);
    let t0 = Instant::now();
    let harness: HarnessFamily = serde_json::from_slice(&std::fs::read(&args[1])?)?;
    let bytes = std::fs::read(&args[2])?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let solution = HarnessSolution::from_slice(&bytes)?;
    drop(bytes);
    println!(
        "read: {} masters, {} rules, fixed {:?}, blake3 {hash}, {:.1} s, rss/peak {:?} GB",
        solution.masters.len(),
        solution.rules.len(),
        solution.fixed,
        t0.elapsed().as_secs_f64(),
        rss_gb()
    );
    let mut names = vec![];
    let targets = match opt("--targets") {
        Some(path) => {
            let config: serde_json::Value = serde_json::from_slice(&std::fs::read(path)?)?;
            let map = config["names"].as_object().ok_or("config has no names")?;
            map.iter()
                .map(|(name, v)| {
                    names.push(name.clone());
                    Ok(Integral(serde_json::from_value::<Vec<i16>>(v.clone())?))
                })
                .collect::<std::result::Result<Vec<_>, serde_json::Error>>()?
        }
        None => {
            names.push("master0".into());
            vec![Integral(solution.masters[0].clone())]
        }
    };
    let epsilon = symbol!("rflift::eps");
    let options = FlowOptions::default();
    let sample = solution.epsilon_sample(options.dimension)?;
    let t1 = Instant::now();
    println!("preparing (exact connection + supplied contract checks) ...");
    let flow = rustred_harness::prepare(
        &harness,
        &solution,
        &targets,
        epsilon,
        &options,
        &RunContext::default(),
        &format!("solution blake3 {hash}"),
        args.iter().any(|a| a == "--homogeneous-eblock"),
    )?;
    let mask = flow.deformation_mask().to_vec();
    let a0 = mask.iter().position(|&b| b).unwrap();
    let free = flow.basis().iter().filter(|m| m.0[a0] <= 0).count();
    let nonzero = flow
        .system()
        .matrix
        .iter()
        .flatten()
        .filter(|a| !a.is_zero())
        .count();
    let blocks = flow.system().blocks()?;
    println!(
        "prepared: eps sample {sample}, basis {} ({} eta-independent), nonzero matrix entries {nonzero}, {} blocks (largest {}), mask {:?}, {:.1} s, rss/peak {:?} GB",
        flow.basis().len(),
        free,
        blocks.len(),
        blocks.iter().map(Vec::len).max().unwrap_or(0),
        mask,
        t1.elapsed().as_secs_f64(),
        rss_gb()
    );
    for (name, t) in names.iter().zip(&targets) {
        println!("target {name} {:?}", t.0);
    }
    let t2 = Instant::now();
    let p = Precision::decimal(digits)?;
    let compiled = flow.system().compile(p, &Default::default())?;
    let mut poles = compiled
        .poles
        .iter()
        .map(|c| (c.re.to_f64(), c.im.to_f64()))
        .collect::<Vec<_>>();
    poles.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!(
        "compiled at {digits} digits: {} finite singular points in eta, {:.1} s, rss/peak {:?} GB",
        poles.len(),
        t2.elapsed().as_secs_f64(),
        rss_gb()
    );
    for (re, im) in &poles {
        println!(
            "  pole eta = {re:.12e} {im:+.3e}i   (z = 1 + eta = {:.12e})",
            1.0 + re
        );
    }
    println!(
        "total {:.1} s, rss/peak {:?} GB",
        t0.elapsed().as_secs_f64(),
        rss_gb()
    );
    Ok(())
}
