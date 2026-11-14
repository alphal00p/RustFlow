//! Benchmark/cross-check: the SAME path segment through RustFlow's existing
//! adaptive certified transport (`DifferentialSystem::compile` +
//! `CompiledSystem::transport`, Symbolica floats rounded per operation), for
//! comparison with the fixed-precision patch kernel.
//!
//! usage: compiled_segment FAMILY.json SOLUTION.json PATH.json --segment A:B --digits D
//!        [--order N] [--target-digits T] [--max-seconds S] --out RESULT.json
//! The start vector is the patch_transport 'ramp'; waypoints are the patch end
//! points of the segment (path variable z = 1 + eta). RESULT.json has the
//! patch-kernel state format, so `patch_transport compare` reads it.
use rug::{Complex, Rational};
use std::sync::Arc;
use std::time::Instant;
use symbolica_amflow::options::CancellationToken;
use symbolica_amflow::patch_kernel::checkpoint::Fingerprint;
use symbolica_amflow::patch_kernel::{PlannedPath, TransportState};
use symbolica_amflow::rustred_harness::{self, HarnessFamily, HarnessSolution};
use symbolica_amflow::*;

type R<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn rss_gb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .unwrap_or_default()
        .lines()
        .find(|l| l.starts_with("VmHWM:"))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.0)
        / 1e6
}

fn main() -> R<()> {
    let args: Vec<String> = std::env::args().collect();
    let opt = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .map(|i| args[i + 1].clone())
    };
    let t0 = Instant::now();
    let harness: HarnessFamily = serde_json::from_slice(&std::fs::read(&args[1])?)?;
    let solution = HarnessSolution::from_slice(&std::fs::read(&args[2])?)?;
    let path_json: serde_json::Value = serde_json::from_slice(&std::fs::read(&args[3])?)?;
    let mut path = PlannedPath::from_json(path_json, &Rational::from(1))?;
    let seg = opt("--segment").ok_or("--segment")?;
    let (a, b) = seg.split_once(':').ok_or("A:B")?;
    path = path.segment(a.parse()?..b.parse()?)?;
    let digits: u32 = opt("--digits").ok_or("--digits")?.parse()?;
    let (_, mask) =
        rustred_harness::family_from_harness(&harness, symbolica::symbol!("rflift::eps"), 4)?;
    let reduced = rustred_harness::reduced_system_with(&solution, &mask, &[], true)?;
    eprintln!(
        "reduced e-block {} ({:.1} s, peak {:.2} GB)",
        reduced.basis.len(),
        t0.elapsed().as_secs_f64(),
        rss_gb()
    );
    let system = DifferentialSystem {
        variable: rustred_harness::eta_symbol(),
        matrix: reduced.matrix,
    };
    let p = Precision::decimal(digits)?;
    let t1 = Instant::now();
    let compiled = system.compile(p, &Default::default())?;
    let t_compile = t1.elapsed().as_secs_f64();
    eprintln!(
        "compiled at {digits} digits ({} bits): {:.1} s, peak {:.2} GB",
        p.bits,
        t_compile,
        rss_gb()
    );
    let dim = compiled.dimension();
    let sf = |x: &Rational| symbolica::domains::float::Float::from(rug::Float::with_val(p.bits, x));
    let to_c = |re: &Rational, im: &Rational| ComplexFloat::new(sf(re), sf(im));
    let one = Rational::from(1);
    let start = &path.patches[0].centre;
    let boundary = BoundaryData {
        point: to_c(&Rational::from(&start.0 - &one), &start.1),
        values: (0..dim)
            .map(|j| {
                to_c(
                    &Rational::from((8 + (5 * j as i64) % 13, 8)),
                    &Rational::new(),
                )
            })
            .collect(),
    };
    let waypoints: Vec<ComplexFloat> = path
        .patches
        .iter()
        .map(|s| {
            let (re, im) = s.end();
            to_c(&Rational::from(&re - &one), &im)
        })
        .collect();
    let mut options = FlowOptions::default();
    options.digits = opt("--target-digits")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(digits.saturating_sub(30).max(10));
    options.guard_digits = 30;
    options.series_order = opt("--order")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(300);
    options.max_steps = 100_000;
    let cancel = CancellationToken::default();
    if let Some(s) = opt("--max-seconds") {
        let secs: f64 = s.parse()?;
        let c = cancel.clone();
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_secs_f64(secs));
            c.cancel();
        });
    }
    let t_start = Instant::now();
    let progress: Arc<dyn Fn(Progress) + Send + Sync> = Arc::new(move |e| {
        if let Progress::Step { index } = e {
            eprintln!(
                "step {index} at {:.1} s, peak {:.2} GB",
                t_start.elapsed().as_secs_f64(),
                rss_gb()
            );
        }
    });
    let context = RunContext {
        cancellation: cancel,
        progress: Some(progress),
    };
    let t2 = Instant::now();
    let result = compiled.transport(&boundary, &waypoints, &options, &context);
    let t_transport = t2.elapsed().as_secs_f64();
    match &result {
        Ok(r) => eprintln!(
            "transport done: {:.1} s, steps {} rejected {} order {} bits {}, peak {:.2} GB",
            t_transport,
            r.diagnostics.steps,
            r.diagnostics.rejected_steps,
            r.diagnostics.expansion_order,
            r.diagnostics.working_bits,
            rss_gb()
        ),
        Err(e) => eprintln!("transport stopped after {t_transport:.1} s: {e}"),
    }
    let r = result?;
    let y: Vec<Complex> = r
        .values
        .iter()
        .map(|c| Complex::with_val(p.bits, (c.re.as_raw(), c.im.as_raw())))
        .collect();
    let state = TransportState {
        next_patch: path.patches.len(),
        y,
        integral: None,
        probe: None,
        records: vec![],
        overlap_log10: vec![-1e9; dim],
        tail_log10: vec![-1e9; dim],
    };
    let fp = Fingerprint {
        kernel_version: "compiled-system-transport".into(),
        system: String::new(),
        path: path.fingerprint(),
        bits: p.bits,
        mode: "ode.rs".into(),
        cauchy_pole: None,
        project_real_axis: false,
        overlap_check: false,
        initial: String::new(),
        binary: String::new(),
    };
    let mut json = state.to_json(&fp);
    json["summary"] = serde_json::json!({"seconds_compile": t_compile, "seconds_transport": t_transport,
        "steps": r.diagnostics.steps, "rejected": r.diagnostics.rejected_steps, "peak_rss_gb": rss_gb(),
        "series_order": options.series_order, "target_digits": options.digits});
    std::fs::write(opt("--out").ok_or("--out")?, serde_json::to_string(&json)?)?;
    Ok(())
}
