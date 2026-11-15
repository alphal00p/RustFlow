//! PFLOW route (B), one d0: massless propagator masters from a one-mass RustRed DE. Thin wrapper of `pflow::cli::run`; the same run
//! is the propagator configuration of `examples/dispersive.rs` (class config with "kind": "propagator").
//!   cargo run --release --example pflow -- --de solution-<d0>.json --graph family.json --e 0 --loops 3 --w 80 --wf 100 --m 120 \
//!       [--moments 1,2] [--s0 1/4] [--S 4] [--r 1/4] [--threads 1] [--two-precision 15] [--provided LOWER.json] --out out.json
//! family.json = one family of RustRed's p catalogue schema (lines[].endpoints, external.in/out, loops) or {"catalogue": PATH,
//! "name": NAME}. Output: for every e-block master and moment m, FP int rho s^-m ds in the Euclidean E-normalisation (m = 1:
//! the massless value), with two-precision digits when requested.
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    match symbolica_amflow::pflow::cli::run(&a) {
        Ok(v) => {
            let s = serde_json::to_string_pretty(&v).unwrap();
            match a.iter().position(|x| x == "--out") {
                Some(i) => std::fs::write(&a[i + 1], &s).unwrap(),
                None => println!("{s}"),
            }
        }
        Err(e) => {
            eprintln!("pflow: {e}");
            std::process::exit(1)
        }
    }
}
