//! RustFlow DISPERSIVE MODE, one sample (W3a): native plan (singular points, mod-p family structure, path), soft boundary
//! (V solved natively from the 4L inputs, or the prepared V file), transport (stand-in kernel behind the W2 interface),
//! matching at infinity, dispersion integral, export in the Store schema (out.json I vector) + J = -(sin(pi d0/2)/pi) I.
//! usage: cargo run --release --example dispersive -- DE.json V.json|native:PLAN.json,SOFTSOL.json,SAMPLE.json W WF NTHREADS WORKDIR OUT.json [REF.json ...]
use rug::Float;
use std::time::Instant;
use symbolica_amflow::dispersive::{boundary, harness, run, transport};
use symbolica_amflow::rustred_harness::HarnessSolution;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (w, wf, nt): (u32, u32, usize) = (a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap());
    let t0 = Instant::now();
    let sol = HarnessSolution::from_slice(&std::fs::read(&a[1]).unwrap()).unwrap();
    let es: usize = std::env::var("DISPERSIVE_E_SLOT").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
    let sys = harness::eblock(&sol, es).unwrap();
    let cfg = run::Config { w, wf, nthreads: nt, workdir: a[6].clone(), ..Default::default() };
    std::fs::create_dir_all(&cfg.workdir).unwrap();
    let plan = run::plan(&sys, &cfg).unwrap();
    let pf = run::bits(wf);
    let (v, vsrc) = if let Some(spec) = a[2].strip_prefix("native:") {
        let f: Vec<&str> = spec.split(',').collect();
        let rd = |p: &str| -> serde_json::Value { serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap() };
        let smp = rd(f[2]);
        let ep: Vec<(String, Vec<i64>, Float)> = smp["endpoint"].as_array().unwrap().iter().map(|x| {
            let t: Vec<i64> = x["target"].as_array().unwrap().iter().map(|v| v.as_i64().unwrap()).collect();
            let key = format!("{}{}", x.get("family").and_then(|f| f.as_str()).unwrap_or(""), serde_json::to_string(&t).unwrap().replace(',', ", "));
            (key, t, Float::with_val(pf + 64, Float::parse(x["value"]["re"].as_str().unwrap()).unwrap()))
        }).collect();
        let (v, d) = boundary::solve_v(&rd(f[0]), &rd(f[1]), &ep, 0, pf + 64).unwrap();
        (v, serde_json::json!({"V": "native (soft_gen solve port)", "plan": f[0], "soft_solution": f[1], "sample": f[2], "diag": d}))
    } else {
        let (v, m) = boundary::read_v_file(&a[2]).unwrap();
        (v, serde_json::json!({"V": "read", "file": a[2], "meta": m}))
    };
    let bnd = boundary::rows_from_v(&sys, &v, 5, pf).unwrap();
    let kernel = transport::TaylorRecurrence::new(&sys, nt);
    let prov = serde_json::json!({"de": a[1], "boundary": vsrc, "transport": "TaylorRecurrence stand-in (W3a)", "patches": plan.patches.len()});
    let out = run::run_sample(&sys, &cfg, &plan, &bnd, &kernel, prov).unwrap();
    std::fs::write(&a[7], serde_json::to_string_pretty(&out).unwrap()).unwrap();
    println!("wall {:.1} s; boundary_fit {} matching {} transport {}", t0.elapsed().as_secs_f64(), out["boundary_fit"], out["matching"], out["transport"]);
    for r in &a[8..] {
        let o: serde_json::Value = serde_json::from_slice(&std::fs::read(r).unwrap()).unwrap();
        let dig = |i: usize| -> f64 {
            let x = Float::with_val(2000, Float::parse(out["I"][i]["I"].as_str().unwrap()).unwrap());
            let y = Float::with_val(2000, Float::parse(o["I"][i]["I"].as_str().unwrap()).unwrap());
            let d = Float::with_val(2000, (x - &y) / y).abs();
            if d.is_zero() { 99.0 } else { -d.log10().to_f64() }
        };
        let all: Vec<f64> = (0..sys.dim()).map(dig).collect();
        let mut s = all.clone();
        s.sort_by(|x, y| x.partial_cmp(y).unwrap());
        println!("vs {r}: idx 564..568 digits {:?}; all entries min {:.1} median {:.1}", (564..569).map(|i| (all[i] * 10.0).round() / 10.0).collect::<Vec<_>>(), s[0], s[s.len() / 2]);
    }
}
