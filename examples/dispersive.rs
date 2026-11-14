//! RustFlow DISPERSIVE MODE, one sample (W3a). Usage and options: src/dispersive/cli.rs. Example (class 16359, node a1 idx 022):
//!   cargo run --release --example dispersive -- --de classes/16359/de/solution-299999_75000.json --class-config classes/16359/config.json \
//!     --m0 24 --minf 50 --w 60 --wf 200 --threads 4 --workdir RUN --out RUN/out.json \
//!     --v native:SOFTPLAN.json,classes/16359/soft/solution-299999_75000.json,v9a/pr/samples/attempt-1/sample-022.json [--pr-values LANE/out.json]
use symbolica_amflow::dispersive::{cli, harness};
use symbolica_amflow::rustred_harness::HarnessSolution;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = cli::Args::parse(&argv).unwrap_or_else(|e| panic!("{e}"));
    let (e_slot, _, _) = args.class_settings().unwrap_or_else(|e| panic!("{e}"));
    let sol = HarnessSolution::from_slice(&std::fs::read(&args.de).unwrap()).unwrap();
    let sys = harness::eblock(&sol, e_slot).unwrap();
    match cli::run(&args, &sys) {
        Ok(out) => println!("{}", serde_json::json!({"boundary_fit": out["boundary_fit"], "matching": out["matching"], "transport": out["transport"],
            "named": out.get("named"), "oracle_comparison": out["oracle_comparison"], "wall_seconds": out["wall_seconds"]})),
        Err(e) => {
            eprintln!("dispersive: {e}");
            std::process::exit(1);
        }
    }
}
