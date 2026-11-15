//! PFLOW route (B), one d0: massless propagator masters from a one-mass RustRed DE.
//!   cargo run --release --example pflow -- --de solution-<d0>.json --graph family.json --e 0 --loops 3 --w 80 --wf 100 --m 120 \
//!       [--moments 1,2] [--s0 1/4] [--S 4] [--r 1/4] [--threads 1] [--two-precision 15] [--provided LOWER.json] --out out.json
//! family.json = one family of RustRed's p catalogue schema (lines[].endpoints, external.in/out, loops) or {"catalogue": PATH,
//! "name": NAME}. Output: for every e-block master and moment m, FP int rho s^-m ds in the Euclidean E-normalisation (m = 1:
//! the massless value), with two-precision digits when requested.
use rug::Rational;
use serde_json::{Value, json};
use symbolica_amflow::dispersive::system::{EBlockSystem, RuleSet, parse_coefficient};
use symbolica_amflow::pflow::engine::{PConfig, run_b};
use symbolica_amflow::pflow::graph::Graph;

fn arg(a: &[String], k: &str) -> Option<String> {
    a.iter().position(|x| x == k).map(|i| a[i + 1].clone())
}
fn load_graph(p: &str) -> Result<Graph, String> {
    let v: Value = serde_json::from_str(&std::fs::read_to_string(p).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    let fam = if let Some(c) = v.get("catalogue") {
        let cat: Value = serde_json::from_str(&std::fs::read_to_string(c.as_str().unwrap()).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        cat["families"].as_array().unwrap().iter().find(|f| f["name"] == v["name"]).cloned().ok_or("family not in catalogue")?
    } else {
        v
    };
    let edges = fam["lines"].as_array().ok_or("no lines")?.iter().map(|l| (l["endpoints"][0].as_u64().unwrap() as usize, l["endpoints"][1].as_u64().unwrap() as usize)).collect();
    Ok(Graph { edges, ext: (fam["external"]["in"].as_u64().unwrap() as usize, fam["external"]["out"].as_u64().unwrap() as usize), loops: fam["loops"].as_u64().unwrap() as usize })
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let r = (|| -> Result<Value, String> {
        let de: Value = serde_json::from_str(&std::fs::read_to_string(arg(&a, "--de").ok_or("--de")?).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let var = de["provenance"]["free"].as_str().unwrap_or("z").to_string();
        let rs = RuleSet::from_json(&de)?;
        let e: usize = arg(&a, "--e").ok_or("--e")?.parse().map_err(|_| "--e")?;
        let sys = EBlockSystem::from_rules_with(&rs, e, &|c: &str| parse_coefficient(c, &var))?;
        let g = load_graph(&arg(&a, "--graph").ok_or("--graph")?)?;
        let num = |k: &str, d: &str| arg(&a, k).unwrap_or(d.into());
        let mut cfg = PConfig::new(num("--w", "60").parse().unwrap(), num("--wf", "80").parse().unwrap(), num("--m", "80").parse().unwrap(), num("--loops", &g.loops.to_string()).parse().unwrap());
        cfg.s0 = num("--s0", "1/4").parse::<Rational>().map_err(|e| e.to_string())?;
        cfg.s_end = num("--S", "4").parse::<Rational>().map_err(|e| e.to_string())?;
        cfg.r = num("--r", "1/4").parse::<Rational>().map_err(|e| e.to_string())?;
        cfg.threads = num("--threads", "1").parse().unwrap();
        cfg.moments = num("--moments", "1").split(',').map(|x| x.parse().unwrap()).collect();
        cfg.nlo = num("--nlo", "-10").parse().unwrap();
        if let Some(f) = arg(&a, "--provided") {
            // lower-loop massless values (E-normalisation, P^2 = 1) for graphs the evaluator cannot reduce
            let v: Value = serde_json::from_str(&std::fs::read_to_string(&f).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
            cfg.provided = symbolica_amflow::pflow::graph::Provided::from_json(&v, symbolica_amflow::dispersive::run::bits(cfg.wf))?;
        }
        let res = run_b(&sys, &g, e, &cfg)?;
        let tp: i64 = num("--two-precision", "0").parse().unwrap();
        let res2 = if tp > 0 {
            let mut c2 = cfg.clone();
            c2.w -= tp as u32;
            c2.wf -= tp as u32;
            c2.m0 -= tp;
            c2.minf -= tp;
            c2.mthr -= tp;
            Some(run_b(&sys, &g, e, &c2)?)
        } else {
            None
        };
        let mut out = Vec::new();
        for (k, (m, v)) in res.moments.iter().enumerate() {
            for (j, mast) in sys.masters.iter().enumerate() {
                let sgn = if mast.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
                let val = rug::Float::with_val(v[j].prec(), &v[j] * sgn);
                let dig = res2.as_ref().map(|r2| {
                    let w = &r2.moments[k].1[j];
                    let d = rug::Float::with_val(v[j].prec(), &v[j] - w);
                    if d.is_zero() { 999.0 } else { -(d.abs() / rug::Float::with_val(v[j].prec(), v[j].abs_ref())).log10().to_f64() }
                });
                out.push(json!({"m": m, "master": mast, "value": val.to_string_radix(10, Some((cfg.wf as usize).max(20))), "digits_two_precision": dig}));
            }
        }
        // non-master targets: --extra "a1,a2,...;b1,..." (pinched-master values from --efree FILE {"[..]": "value"} if needed)
        let mut extras = Vec::new();
        if let Some(ex) = arg(&a, "--extra") {
            let efree: std::collections::HashMap<Vec<i64>, rug::Float> = match arg(&a, "--efree") {
                Some(f) => {
                    let v: Value = serde_json::from_str(&std::fs::read_to_string(f).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
                    v.as_object().unwrap().iter().map(|(k, x)| (serde_json::from_str::<Vec<i64>>(k).unwrap(), rug::Float::with_val(symbolica_amflow::dispersive::run::bits(cfg.wf), rug::Float::parse(x.as_str().unwrap()).unwrap()))).collect()
                }
                None => Default::default(),
            };
            for t in ex.split(';') {
                let xv: Vec<i64> = t.split(',').map(|y| y.trim().parse().unwrap()).collect();
                let prec = symbolica_amflow::dispersive::run::bits(cfg.wf);
                let r1 = symbolica_amflow::pflow::engine::extra_value(&rs.rules, &sys, &xv, &res.moments, &efree, &var, prec);
                let dig = match (&r1, &res2) {
                    (Ok((v, _)), Some(r2)) => symbolica_amflow::pflow::engine::extra_value(&rs.rules, &sys, &xv, &r2.moments, &efree, &var, prec).ok().map(|(w, _)| {
                        let d = rug::Float::with_val(prec, v - &w);
                        if d.is_zero() { 999.0 } else { -(d.abs() / rug::Float::with_val(prec, v.abs_ref())).log10().to_f64() }
                    }),
                    _ => None,
                };
                extras.push(match r1 {
                    Ok((v, n)) => json!({"target": xv, "value": v.to_string_radix(10, Some((cfg.wf as usize).max(20))), "max_moment_n": n, "digits_two_precision": dig}),
                    Err(e) => json!({"target": xv, "error": e}),
                });
            }
        }
        Ok(json!({"schema": "pflow-route-b.v1", "extras": extras, "d0": rs.d0.to_string(), "e_slot": e, "config": {"w": cfg.w, "wf": cfg.wf, "m": cfg.m0, "s0": cfg.s0.to_string(), "S": cfg.s_end.to_string(), "r": cfg.r.to_string(), "moments": cfg.moments},
            "values": out, "diag": res.diag, "diag_second": res2.map(|r| r.diag)}))
    })();
    match r {
        Ok(v) => {
            let s = serde_json::to_string_pretty(&v).unwrap();
            if let Some(o) = arg(&a, "--out") { std::fs::write(o, &s).unwrap() } else { println!("{s}") }
        }
        Err(e) => {
            eprintln!("pflow: {e}");
            std::process::exit(1)
        }
    }
}
