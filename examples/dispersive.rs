//! RustFlow DISPERSIVE MODE, one sample. Usage and options: src/dispersive/cli.rs. Example (class 16359, node a1 idx 022):
//!   cargo run --release --example dispersive -- --de classes/16359/de/solution-299999_75000.json --class-config classes/16359/config.json \
//!     --m0 24 --minf 50 --w 165 --wf 205 --threads 16 --workdir RUN --out RUN/out.json --two-precision -15,-15 \
//!     --v native:SOFTPLAN.json,classes/16359/soft/solution-299999_75000.json,v9a/pr/samples/attempt-1/sample-022.json [--pr-values LANE/out.json]
//!
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
                        eprintln!("dispersive: signal outside the transport: exit {EXIT_RESUMABLE} (stage files are atomic; rerun the same command to resume)");
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
        cli::run_with(&args, &sys, &cancel, &transport_active)
    };
    match run() {
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
