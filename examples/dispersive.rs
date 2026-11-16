//! RustFlow DISPERSIVE MODE, one sample. Usage and options: src/dispersive/cli.rs. Example (class 16359, node a1 idx 022):
//!   cargo run --release --example dispersive -- --de classes/16359/de/solution-299999_75000.json --class-config classes/16359/config.json \
//!     --m0 24 --minf 50 --w 165 --wf 205 --threads 16 --workdir RUN --out RUN/out.json --two-precision -15,-15 \
//!     --v native:SOFTPLAN.json,classes/16359/soft/solution-299999_75000.json,v9a/pr/samples/attempt-1/sample-022.json [--pr-values LANE/out.json]
//!
//! `dispersive --system-info DE.json CLASS_CONFIG.json` prints the e-block size and the kernel system fingerprint.
//! Propagator configuration: a class config with "kind": "propagator" (graph, e, w, wf, m, moments, ...; src/pflow/cli.rs) runs
//! route (B) for massless p-integrals with the same command (`dispersive --de DE.json --class-config P.json --out OUT.json`).
//! `dispersive --progress WORKDIR` prints one JSON progress line (cli::progress) and exits.
//!
//! Everything class-specific comes from the DE file, the class config and the command line. Restarting the same command resumes
//! from the stage files and the transport checkpoint in the workdir (any thread count); a workdir of another input/configuration is
//! refused. SIGTERM/SIGINT: inside the transport the kernel finishes its check, flushes the checkpoint and the run exits with 75;
//! outside it the process exits with 75 at once (stage files are written atomically; the interrupted stage is recomputed).
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use symbolica_amflow::dispersive::{cli, harness, run};
use symbolica_amflow::options::CancellationToken;
use symbolica_amflow::rustred_harness::HarnessSolution;

const EXIT_RESUMABLE: i32 = 75;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.len() == 2 && argv[0] == "--progress" {
        // progress from the checkpoint files of a workdir (drivers; never ahead of what a restart resumes)
        println!("{}", cli::progress(&argv[1]));
        return;
    }
    if argv.len() == 3 && argv[0] == "--system-info" {
        // DE + class config -> e-block through the harness reader, kernel system fingerprint (checks after reader changes)
        let t = std::time::Instant::now();
        let cfg: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&argv[2]).expect("class config"))
                .expect("class config JSON");
        let e = cfg["e"].as_u64().expect("e slot") as usize;
        let sol = HarnessSolution::from_slice(&std::fs::read(&argv[1]).expect("DE"))
            .unwrap_or_else(|err| panic!("reader: {err}"));
        let sys = harness::eblock(&sol, e).unwrap_or_else(|err| panic!("e-block: {err}"));
        let k =
            symbolica_amflow::dispersive::kernel::PatchKernel::new(&sys, 1).expect("kernel system");
        println!(
            "{}",
            serde_json::json!({"de": argv[1], "d0": sys.d0.to_string(), "e_slot": e, "dim": sys.dim(), "terms": sys.terms.len(),
            "kernel_system_fingerprint": k.system_fingerprint, "products_per_order": k.system.products_per_order(), "seconds": t.elapsed().as_secs_f64()})
        );
        return;
    }
    // propagator configuration (class config "kind": "propagator"): massless p-integrals by route (B), src/pflow/cli.rs
    if let Some(i) = argv.iter().position(|x| x == "--class-config") {
        if argv.get(i + 1).is_some_and(|c| symbolica_amflow::pflow::cli::is_propagator_config(c)) {
            match symbolica_amflow::pflow::cli::run(&argv) {
                Ok(v) => {
                    let s = serde_json::to_string_pretty(&v).unwrap();
                    match argv.iter().position(|x| x == "--out") {
                        Some(j) => std::fs::write(&argv[j + 1], &s).expect("--out"),
                        None => println!("{s}"),
                    }
                    let code = symbolica_amflow::pflow::cli::exit_code(&v);
                    eprintln!("dispersive (propagator): {}", symbolica_amflow::pflow::cli::summary_line(&v));
                    if code != 0 {
                        eprintln!("dispersive (propagator): NOT ACCEPTED (output written; failures in diag.acceptance)");
                        std::process::exit(code)
                    }
                    return;
                }
                Err(e) => {
                    eprintln!("dispersive (propagator): {e}");
                    std::process::exit(1)
                }
            }
        }
    }
    let args = cli::Args::parse(&argv).unwrap_or_else(|e| {
        eprintln!("dispersive: {e}");
        std::process::exit(2)
    });
    let flag = Arc::new(AtomicBool::new(false));
    for sig in [signal_hook::consts::SIGTERM, signal_hook::consts::SIGINT] {
        signal_hook::flag::register(sig, flag.clone()).expect("signal handler");
    }
    let cancel = CancellationToken::from_flag(flag.clone());
    let transport_active = Arc::new(AtomicBool::new(false));
    {
        let (flag, active) = (flag.clone(), transport_active.clone());
        std::thread::spawn(move || {
            loop {
                std::thread::sleep(std::time::Duration::from_millis(200));
                if flag.load(Ordering::SeqCst) && !active.load(Ordering::SeqCst) {
                    // a signal outside the checkpointed transport: give a just-starting transport a moment to own it
                    std::thread::sleep(std::time::Duration::from_millis(500));
                    if !active.load(Ordering::SeqCst) {
                        eprintln!(
                            "dispersive: signal outside the transport: exit {EXIT_RESUMABLE} (stage files are atomic; rerun the same command to resume)"
                        );
                        std::process::exit(EXIT_RESUMABLE);
                    }
                }
            }
        });
    }
    let run = || -> Result<serde_json::Value, String> {
        let (e_slot, _, _) = args.class_settings()?;
        let bytes = std::fs::read(&args.de).map_err(|e| format!("{}: {e}", args.de))?;
        let sol = HarnessSolution::from_slice(&bytes).map_err(|e| format!("{}: {e}", args.de))?;
        let sys = harness::eblock(&sol, e_slot).map_err(|e| e.to_string())?;
        drop(sol);
        // optional lane stage: regenerate the 4L sample's PR rows from its P14 rows first (--lane)
        let load = |de: &str,
                    e: usize|
         -> Result<symbolica_amflow::dispersive::system::EBlockSystem, String> {
            let b = std::fs::read(de).map_err(|e| format!("{de}: {e}"))?;
            let s = HarnessSolution::from_slice(&b).map_err(|e| format!("{de}: {e}"))?;
            harness::eblock(&s, e).map_err(|e| e.to_string())
        };
        let (args5, _) = cli::lane_stage(&args, &sys.d0, &load, &cancel, &transport_active)?;
        if args.lane_only {
            // hidden diagnostic flag --lane-only 1: stop after the lane stage
            return Ok(serde_json::json!({"lane_stage": args5.lane_info, "v": args5.v}));
        }
        cli::run_with(&args5, &sys, &cancel, &transport_active)
    };
    match run() {
        Ok(out) if out.get("lane_stage").is_some() => println!("{out}"),
        Ok(out) if out.get("family_check").is_some() => println!("{out}"),
        Ok(out) => println!(
            "{}",
            serde_json::json!({"boundary_fit": out["boundary_fit"], "matching": out["matching"], "transport": out["transport"],
            "residuals": out["residuals"], "digits_summary": out.get("digits_summary"), "named": out.get("named"),
            "oracle_comparison_min_median": out["oracle_comparison"].as_object().map(|m| m.iter().map(|(k, v)| (k.clone(), serde_json::json!([v["min"], v["median"], v["named"]]))).collect::<serde_json::Map<_, _>>()),
            "cost": out["cost"], "wall_seconds": out["wall_seconds"]})
        ),
        Err(e) if e.starts_with(run::CANCELLED) => {
            eprintln!("dispersive: {e}");
            std::process::exit(EXIT_RESUMABLE);
        }
        Err(e) => {
            eprintln!("dispersive: {e}");
            std::process::exit(1);
        }
    }
}
