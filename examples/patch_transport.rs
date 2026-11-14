//! Driver of RustFlow's fixed-precision patch-transport kernel (`patch_kernel`).
//!
//! usage:
//!   patch_transport system-harness FAMILY.json SOLUTION.json --out SYS.json [--full]
//!       RustRed harness solution -> exact system (homogeneous e-block unless --full);
//!       also writes SYS.json.basis.json (master index vectors in basis order)
//!   patch_transport system-terms INP.json --out SYS.json
//!       campaign in-house dispersive input (terms, masters) -> exact system (+ basis file)
//!   patch_transport system-diff A.json B.json
//!       exact entry-by-entry comparison (basis files give the permutation)
//!   patch_transport plan PLAN.json --digits W --out PATH.json [--offset 1] [--check ZPATH.json] [--system SYS.json]
//!       port of zplan.py (PLAN.json: s0, s_end, singular_points, axis_singular, ratio); --system checks that
//!       every singular point of the system (denominator roots, z = x + offset) is in PLAN.json
//!   patch_transport run --config dispersive|generic --system SYS.json --path PATH.json [--segment A:B]
//!       [--orders-digits W --plan PLAN.json] (--digits W | --bits B) [--threads T] [--mode fixed|tracked]
//!       [--initial ramp|unit|FILE] [--no-overlap] [--no-pin]
//!       generic only: [--cauchy-pole A] [--project-real]   (dispersive = real projection + pole 1)
//!       [--checkpoint FILE [--every-patches N] [--every-seconds S] [--fresh] [--require-same-binary]]
//!       [--stop-after K] --out RESULT.json
//!       An existing checkpoint is RESUMED (any thread count); --fresh moves it aside first.
//!       SIGTERM/SIGINT -> checkpoint flushed -> exit 75.
//!   patch_transport compare RUN.json REF.json [--out REPORT.json]
//!       per-entry digits (two working precisions; refuses runs that are not comparable)
use rug::{Complex, Rational};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::Instant;
use symbolica_amflow::options::CancellationToken;
use symbolica_amflow::patch_kernel::{
    self, ArithmeticMode, DigitsReport, PlanInput, PlannedPath, RationalSystem, Transport,
    TransportConfig, TransportOutcome, TransportState, checkpoint,
};
use symbolica_amflow::rustred_harness::{self, HarnessFamily, HarnessSolution};

type R<T> = std::result::Result<T, Box<dyn std::error::Error>>;

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

struct Args(Vec<String>);
impl Args {
    fn opt(&self, key: &str) -> Option<String> {
        self.0
            .iter()
            .position(|a| a == key)
            .and_then(|i| self.0.get(i + 1).cloned())
    }
    fn flag(&self, key: &str) -> bool {
        self.0.iter().any(|a| a == key)
    }
    fn req(&self, key: &str) -> R<String> {
        self.opt(key).ok_or_else(|| format!("missing {key}").into())
    }
    fn pos(&self, k: usize) -> R<String> {
        self.0
            .get(k)
            .cloned()
            .ok_or_else(|| format!("missing positional argument {k}").into())
    }
}

fn read_json(path: &str) -> R<serde_json::Value> {
    Ok(serde_json::from_slice(&std::fs::read(path)?)?)
}

fn write_json(path: &str, v: &serde_json::Value) -> R<()> {
    checkpoint::write_atomic(
        std::path::Path::new(path),
        serde_json::to_string(v)?.as_bytes(),
    )?;
    Ok(())
}

fn save_system(sys: &RationalSystem, basis: &[Vec<i64>], out: &str, t0: Instant) -> R<()> {
    let fp = sys.fingerprint();
    write_json(out, &sys.to_json())?;
    write_json(
        &format!("{out}.basis.json"),
        &serde_json::json!({"basis": basis, "fingerprint": fp}),
    )?;
    println!(
        "system: dim {} entries {} distinct polys {} products/order {} fingerprint {fp} ({:.1} s, rss/peak {:?} GB)",
        sys.dim,
        sys.entries.len(),
        sys.polys.len(),
        sys.products_per_order(),
        t0.elapsed().as_secs_f64(),
        rss_gb()
    );
    Ok(())
}

fn system_harness(a: &Args) -> R<()> {
    let t0 = Instant::now();
    let harness: HarnessFamily = serde_json::from_slice(&std::fs::read(a.pos(2)?)?)?;
    let bytes = std::fs::read(a.pos(3)?)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let solution = HarnessSolution::from_slice(&bytes)?;
    drop(bytes);
    println!(
        "read solution blake3 {hash}: {} masters, fixed {:?} ({:.1} s)",
        solution.masters.len(),
        solution.fixed,
        t0.elapsed().as_secs_f64()
    );
    let epsilon = symbolica::symbol!("rflift::eps");
    let (_, mask) = rustred_harness::family_from_harness(&harness, epsilon, 4)?;
    let homogeneous = !a.flag("--full");
    let reduced = rustred_harness::reduced_system_with(&solution, &mask, &[], homogeneous)?;
    println!(
        "reduced (homogeneous e-block {homogeneous}): basis {} ({:.1} s, rss/peak {:?} GB)",
        reduced.basis.len(),
        t0.elapsed().as_secs_f64(),
        rss_gb()
    );
    let sys = RationalSystem::from_matrix_atoms(&reduced.matrix, rustred_harness::eta_symbol())?;
    let basis: Vec<Vec<i64>> = reduced
        .basis
        .iter()
        .map(|i| i.0.iter().map(|&v| v as i64).collect())
        .collect();
    save_system(&sys, &basis, &a.req("--out")?, t0)
}

fn system_terms(a: &Args) -> R<()> {
    let t0 = Instant::now();
    let inp = read_json(&a.pos(2)?)?;
    let dim = inp["E"].as_u64().ok_or("no E")? as usize;
    let sys = RationalSystem::from_terms_json(dim, &inp["terms"])?;
    let basis: Vec<Vec<i64>> = serde_json::from_value(inp["masters"].clone())?;
    save_system(&sys, &basis, &a.req("--out")?, t0)
}

fn system_diff(a: &Args) -> R<()> {
    let (pa, pb) = (a.pos(2)?, a.pos(3)?);
    let sa = RationalSystem::from_json(read_json(&pa)?)?;
    let sb = RationalSystem::from_json(read_json(&pb)?)?;
    let ba: Vec<Vec<i64>> =
        serde_json::from_value(read_json(&format!("{pa}.basis.json"))?["basis"].clone())?;
    let bb: Vec<Vec<i64>> =
        serde_json::from_value(read_json(&format!("{pb}.basis.json"))?["basis"].clone())?;
    let index: BTreeMap<&Vec<i64>, usize> = bb.iter().enumerate().map(|(i, m)| (m, i)).collect();
    let perm = ba
        .iter()
        .map(|m| index.get(m).copied().ok_or("basis differs"))
        .collect::<Result<Vec<_>, _>>()?;
    let diff = sa.differences(&sb, &perm);
    let identity = perm.iter().enumerate().all(|(i, &p)| i == p);
    println!(
        "system-diff: dims {} / {}, entries {} / {}, basis order identical {identity}, differing entries {}",
        sa.dim,
        sb.dim,
        sa.entries.len(),
        sb.entries.len(),
        diff.len()
    );
    for d in diff.iter().take(20) {
        println!("  differs at {d:?}");
    }
    if !diff.is_empty() {
        std::process::exit(1);
    }
    Ok(())
}

fn plan(a: &Args) -> R<()> {
    let input: PlanInput = serde_json::from_value(read_json(&a.pos(2)?)?)?;
    let digits: f64 = a.req("--digits")?.parse()?;
    let offset: Rational = a.opt("--offset").unwrap_or("1".into()).parse()?;
    if let Some(sys) = a.opt("--system") {
        let system = RationalSystem::from_json(read_json(&sys)?)?;
        let t = Instant::now();
        let off = patch_kernel::path::to_f64(&offset);
        let derived: Vec<(f64, f64)> = system
            .singular_points()?
            .into_iter()
            .map(|(re, im)| (re + off, im))
            .collect();
        patch_kernel::path::verify_singular_points(&input.singular_points, &derived, 1e-9)?;
        println!(
            "plan: the {} singular points of the system are all in the planner input ({:.1} s)",
            derived.len(),
            t.elapsed().as_secs_f64()
        );
    }
    let path = patch_kernel::path::plan_axis_detours(&input, digits, offset.clone())?;
    let json = path.to_json();
    println!(
        "plan: {} patches, {} complex, sum of orders {}",
        json["n"], json["n_complex"], json["sum_nord"]
    );
    if let Some(check) = a.opt("--check") {
        let reference = read_json(&check)?;
        let same = reference["patches"] == json["patches"];
        println!("check against {check}: patches identical = {same}");
        if !same {
            let (x, y) = (
                reference["patches"].as_array().unwrap(),
                json["patches"].as_array().unwrap(),
            );
            if let Some(k) = (0..x.len().min(y.len())).find(|&k| x[k] != y[k]) {
                println!(
                    "  first difference at patch {k}: reference {} ours {}",
                    x[k], y[k]
                );
            }
            println!("  lengths {} / {}", x.len(), y.len());
        }
    }
    if let Some(out) = a.opt("--out") {
        write_json(&out, &json)?;
    }
    Ok(())
}

fn initial_vector(spec: &str, dim: usize, bits: u32) -> R<Vec<Complex>> {
    Ok(match spec {
        // deterministic, well-conditioned: y_j = (8 + (5 j mod 13)) / 8
        "ramp" => (0..dim)
            .map(|j| Complex::with_val(bits, (Rational::from((8 + (5 * j as i64) % 13, 8)), 0)))
            .collect(),
        "unit" => (0..dim).map(|_| Complex::with_val(bits, (1, 0))).collect(),
        file => {
            let v: Vec<(String, String)> = serde_json::from_value(read_json(file)?)?;
            v.iter()
                .map(|(re, im)| {
                    Ok(Complex::with_val(
                        bits,
                        (re.parse::<Rational>()?, im.parse::<Rational>()?),
                    ))
                })
                .collect::<R<Vec<_>>>()?
        }
    })
}

/// SIGTERM and SIGINT set `flag`: the run flushes its checkpoint and exits 75.
fn cancel_on_signals(flag: Arc<AtomicBool>) -> R<()> {
    for sig in [signal_hook::consts::SIGTERM, signal_hook::consts::SIGINT] {
        signal_hook::flag::register(sig, flag.clone())?;
    }
    Ok(())
}

fn run(a: &Args) -> R<()> {
    let t0 = Instant::now();
    let flag = Arc::new(AtomicBool::new(false));
    cancel_on_signals(flag.clone())?;
    let cancel = CancellationToken::from_flag(flag);
    let sys_path = a.req("--system")?;
    let system = RationalSystem::from_json(read_json(&sys_path)?)?;
    let system_fp = system.fingerprint();
    let t_load = t0.elapsed().as_secs_f64();
    let mut path = PlannedPath::from_json(read_json(&a.req("--path")?)?, &Rational::from(1))?;
    if let Some(seg) = a.opt("--segment") {
        let (x, y) = seg.split_once(':').ok_or("--segment A:B")?;
        path = path.segment(x.parse()?..y.parse()?)?;
    }
    if let Some(w) = a.opt("--orders-digits") {
        let input: PlanInput = serde_json::from_value(read_json(&a.req("--plan")?)?)?;
        path = path.with_orders(&input.singular_points, w.parse()?)?;
    }
    let bits = match (a.opt("--bits"), a.opt("--digits")) {
        (Some(b), _) => b.parse()?,
        (None, Some(d)) => patch_kernel::bits_for_digits(d.parse()?),
        _ => return Err("--digits or --bits".into()),
    };
    let mut config = match a.req("--config")?.as_str() {
        "dispersive" => {
            if a.opt("--cauchy-pole").is_some() || a.flag("--project-real") {
                return Err("--cauchy-pole/--project-real are fixed by --config dispersive".into());
            }
            TransportConfig::dispersive(bits)
        }
        "generic" => {
            let mut c = TransportConfig::generic(bits);
            c.cauchy_pole = a
                .opt("--cauchy-pole")
                .map(|s| s.parse::<Rational>())
                .transpose()?;
            c.project_real_axis = a.flag("--project-real");
            c
        }
        other => return Err(format!("--config {other}: dispersive|generic").into()),
    };
    config.threads = a.opt("--threads").unwrap_or("1".into()).parse()?;
    config.mode = a
        .opt("--mode")
        .unwrap_or("fixed".into())
        .parse::<ArithmeticMode>()?;
    config.overlap_check = !a.flag("--no-overlap");
    config.pin_threads = !a.flag("--no-pin");
    let initial = initial_vector(
        &a.opt("--initial").unwrap_or("ramp".into()),
        system.dim,
        bits,
    )?;
    let transport =
        Transport::new(&system, &path, config.clone())?.with_system_fingerprint(system_fp.clone());
    let policy = a.opt("--checkpoint").map(|f| {
        let mut p = checkpoint::CheckpointPolicy::new(f);
        p.every_patches = a
            .opt("--every-patches")
            .and_then(|v| v.parse().ok())
            .unwrap_or(1);
        p.every_seconds = a
            .opt("--every-seconds")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0);
        p
    });
    if let Some(p) = &policy {
        if a.flag("--fresh") && p.file.exists() {
            let secs = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0);
            let aside = format!("{}.aside-{secs}", p.file.display());
            std::fs::rename(&p.file, &aside)?;
            eprintln!("--fresh: moved the old checkpoint to {aside}");
        }
    }
    let mut state = match &policy {
        Some(p) if p.file.exists() => {
            let s = transport
                .resume(&p.file, &initial, a.flag("--require-same-binary"))
                .map_err(|e| {
                    format!("{e} (use --fresh to start over; the old file is kept aside)")
                })?;
            eprintln!(
                "resumed at patch {} from {}",
                s.next_patch,
                p.file.display()
            );
            s
        }
        _ => transport.initial_state(&initial)?,
    };
    let until = a
        .opt("--stop-after")
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(usize::MAX);
    eprintln!(
        "run: dim {} entries {} patches {} (complex {}) bits {bits} threads {} mode {:?}; load {t_load:.1} s, rss/peak {:?} GB",
        system.dim,
        system.entries.len(),
        path.patches.len(),
        path.patches.iter().filter(|p| !p.is_real()).count(),
        config.threads,
        config.mode,
        rss_gb()
    );
    let t_run = Instant::now();
    let outcome = transport.run(&mut state, until, &initial, &cancel, policy.as_ref(), |s| {
        let r = s.records.last().unwrap();
        eprintln!(
            "patch {} order {} complex {} setup {:.2} s rec {:.2} s (busy {:.2}..{:.2}) fin {:.2} s overlap {:.1} tail {:.1} imag {:.1} rss/peak {:?} GB",
            r.index, r.order, r.complex, r.seconds_setup, r.seconds_recurrence, r.worker_busy_min, r.worker_busy_max, r.seconds_finish,
            r.overlap_log10, r.tail_log10, r.dropped_imag_log10, rss_gb()
        );
    })?;
    let seconds = t_run.elapsed().as_secs_f64();
    let fp = transport.fingerprint(&initial);
    let recs = &state.records;
    let sum = |f: &dyn Fn(&patch_kernel::PatchRecord) -> f64| recs.iter().map(f).sum::<f64>();
    let summary = serde_json::json!({
        "outcome": format!("{outcome:?}"),
        "patches_done": state.next_patch,
        "patches_total": path.patches.len(),
        "bits": bits,
        "threads": config.threads,
        "mode": format!("{:?}", config.mode),
        "seconds_this_run": seconds,
        "seconds_load": t_load,
        "seconds_setup": sum(&|r| r.seconds_setup),
        "seconds_recurrence": sum(&|r| r.seconds_recurrence),
        "seconds_finish": sum(&|r| r.seconds_finish),
        "peak_rss_gb": rss_gb().1,
        "identity_hash": state.identity_hash(),
        "system_fingerprint": system_fp,
        "max_dropped_imag_log10": state.max_dropped_imag_log10(),
        "overlap_skipped": state.records.iter().filter(|r| r.overlap_skipped).count(),
    });
    eprintln!("summary: {summary}");
    if let Some(out) = a.opt("--out") {
        let mut json = state.to_json(&fp);
        json["summary"] = summary;
        write_json(&out, &json)?;
    }
    if outcome == TransportOutcome::Cancelled {
        eprintln!(
            "stopped by signal after patch {}; checkpoint flushed; exit {}",
            state.next_patch,
            checkpoint::EXIT_RESUMABLE
        );
        std::process::exit(checkpoint::EXIT_RESUMABLE);
    }
    Ok(())
}

fn compare(a: &Args) -> R<()> {
    let load = |p: &str| -> R<(TransportState, checkpoint::Fingerprint)> {
        Ok(TransportState::from_json(read_json(p)?)?)
    };
    let (run, f1) = load(&a.pos(2)?)?;
    let (reference, f2) = load(&a.pos(3)?)?;
    let report = DigitsReport::new(&run, &f1, &reference, &f2)?;
    println!(
        "compare: run W~{:.0} vs reference W~{:.0}: agreement min {:.2} median {:.2}; integral min {:?}; overlap digits min {:.2} ({} skipped); dropped imag 10^{:.1}",
        report.digits,
        report.reference_digits,
        report.min_agreement,
        report.median_agreement,
        report.min_integral_agreement,
        report.min_overlap_digits,
        report.overlap_skipped,
        report.dropped_imag_log10
    );
    if let Some(out) = a.opt("--out") {
        write_json(&out, &serde_json::to_value(&report)?)?;
    }
    Ok(())
}

fn main() {
    let args = Args(std::env::args().collect());
    let result = match args.0.get(1).map(String::as_str) {
        Some("system-harness") => system_harness(&args),
        Some("system-terms") => system_terms(&args),
        Some("system-diff") => system_diff(&args),
        Some("plan") => plan(&args),
        Some("run") => run(&args),
        Some("compare") => compare(&args),
        _ => Err("usage: patch_transport system-harness|system-terms|system-diff|plan|run|compare ... (see source header)".into()),
    };
    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
