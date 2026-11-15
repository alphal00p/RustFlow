//! Propagator configuration of the RustFlow dispersive driver (route B, PFLOW-S1*.md). One entry point for both front ends:
//! `examples/dispersive.rs` (class config with "kind": "propagator") and `examples/pflow.rs` (flags only).
//!
//! Options (flag, or the same key without dashes in the class config; flags win):
//!   --de solution-<d0>.json  --graph family.json (or config "graph": path relative to the config, or an inline family)
//!   --e SLOT  --w W --wf WF --m M  [--loops L] [--moments 1,2] [--s0 1/4] [--S 4] [--r 1/4] [--threads 1] [--nlo -10]
//!   [--two-precision dW] [--provided LOWER.json] [--extra "a1,..;b1,.."] [--efree EFREE.json] [--workdir DIR] [--out OUT.json]
//! --de may list several per-cut DE BLOCKS on one common basis (comma list, or a manifest {"blocks": [...]}): see `load_blocks`.
//! --workdir: checkpoints (Frobenius families, every transport after every patch); rerunning the same command resumes.
//! family.json = one family of RustRed's p catalogue schema (lines[].endpoints, external.in/out, loops) or {"catalogue": PATH,
//! "name": NAME}. Output: for every e-block master and moment m, FP int rho s^-m ds in the Euclidean E-normalisation (m = 1:
//! the massless value), with two-precision digits when requested.
use super::engine::{PConfig, run_b};
use super::graph::{Graph, Provided};
use crate::dispersive::system::{EBlockSystem, RuleSet, parse_coefficient};
use rug::Rational;
use serde_json::{Value, json};
use std::path::Path;

/// option lookup: command line first, then the class config (key without the dashes); relative paths in the config are taken
/// relative to the config file
struct Opts<'a> {
    argv: &'a [String],
    cfg: Value,
    dir: std::path::PathBuf,
}

impl Opts<'_> {
    fn get(&self, k: &str) -> Option<String> {
        if let Some(i) = self.argv.iter().position(|x| x == k) {
            return self.argv.get(i + 1).cloned();
        }
        match self.cfg.get(k.trim_start_matches('-')) {
            Some(Value::String(s)) => Some(s.clone()),
            Some(Value::Number(n)) => Some(n.to_string()),
            Some(Value::Array(a)) => Some(a.iter().map(|x| x.to_string().trim_matches('"').to_string()).collect::<Vec<_>>().join(",")),
            _ => None,
        }
    }
    fn path(&self, k: &str) -> Option<String> {
        if self.argv.iter().any(|x| x == k) {
            return self.get(k);
        }
        self.get(k).map(|p| if Path::new(&p).is_absolute() { p } else { self.dir.join(p).to_string_lossy().into_owned() })
    }
}

fn read_json(p: &str) -> Result<Value, String> {
    serde_json::from_str(&std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}"))?).map_err(|e| format!("{p}: {e}"))
}

pub fn graph_from(v: &Value) -> Result<Graph, String> {
    let fam = if let Some(c) = v.get("catalogue") {
        let cat = read_json(c.as_str().ok_or("catalogue: a path")?)?;
        cat["families"].as_array().ok_or("catalogue: families")?.iter().find(|f| f["name"] == v["name"]).cloned().ok_or("family not in catalogue")?
    } else {
        v.clone()
    };
    let edges = fam["lines"].as_array().ok_or("no lines")?.iter().map(|l| (l["endpoints"][0].as_u64().unwrap() as usize, l["endpoints"][1].as_u64().unwrap() as usize)).collect();
    Ok(Graph { edges, ext: (fam["external"]["in"].as_u64().ok_or("external.in")? as usize, fam["external"]["out"].as_u64().ok_or("external.out")? as usize), loops: fam["loops"].as_u64().ok_or("loops")? as usize })
}

/// true when `class_config` selects this configuration ("kind": "propagator")
pub fn is_propagator_config(class_config: &str) -> bool {
    read_json(class_config).map(|v| v["kind"] == "propagator").unwrap_or(false)
}

/// the DE of one d0: one RustRed harness solution, or several BLOCKS (per-cut reductions on one common master basis) given as a
/// comma list or as a manifest JSON {"blocks": [paths relative to it]}. Blocks are united: masters and rules; the same target in two
/// blocks must have the identical rule (else an error names it); all blocks must share d0 and the free variable.
pub fn load_blocks(spec: &str) -> Result<(RuleSet, String), String> {
    let paths: Vec<String> = if spec.ends_with(".json") && !spec.contains(',') {
        let v = read_json(spec)?;
        match v.get("blocks").and_then(|b| b.as_array()) {
            Some(b) => {
                let dir = Path::new(spec).parent().map(|d| d.to_path_buf()).unwrap_or_default();
                b.iter().map(|x| { let p = x.as_str().unwrap_or(""); if Path::new(p).is_absolute() { p.to_string() } else { dir.join(p).to_string_lossy().into_owned() } }).collect()
            }
            None => vec![spec.to_string()],
        }
    } else {
        spec.split(',').map(|x| x.to_string()).collect()
    };
    let mut out: Option<(RuleSet, String)> = None;
    for p in &paths {
        let de = read_json(p)?;
        let var = de["provenance"]["free"].as_str().unwrap_or("z").to_string();
        let rs = RuleSet::from_json(&de)?;
        match &mut out {
            None => out = Some((rs, var)),
            Some((acc, v0)) => {
                if acc.d0 != rs.d0 || *v0 != var {
                    return Err(format!("{p}: d0 {} / variable {var} differ from the first block ({} / {v0})", rs.d0, acc.d0));
                }
                for m in rs.masters {
                    if !acc.masters.contains(&m) {
                        acc.masters.push(m);
                    }
                }
                for (t, r) in rs.rules {
                    match acc.rules.get(&t) {
                        Some(r0) if *r0 != r => {
                            // same rule up to term order?
                            let mut a = r0.clone();
                            let mut b = r.clone();
                            a.sort();
                            b.sort();
                            if a != b {
                                return Err(format!("{p}: the rule for {t:?} differs from an earlier block (not one common basis?)"));
                            }
                        }
                        Some(_) => {}
                        None => {
                            acc.rules.insert(t, r);
                        }
                    }
                }
            }
        }
    }
    let (rs, var) = out.ok_or("no DE block")?;
    // one common basis: no master of one block may be reduced by a rule of another
    if let Some(m) = rs.masters.iter().find(|m| rs.rules.contains_key(*m)) {
        return Err(format!("master {m:?} of one block is reduced by a rule of another block: the blocks are not on one common basis"));
    }
    eprintln!("[pflow] DE from {} block(s): {} masters, {} rules", paths.len(), rs.masters.len(), rs.rules.len());
    Ok((rs, var))
}

/// run one d0 of the propagator configuration; `argv` without the program name
pub fn run(argv: &[String]) -> Result<Value, String> {
    let (cfgv, dir) = match argv.iter().position(|x| x == "--class-config") {
        Some(i) => {
            let p = argv.get(i + 1).ok_or("--class-config needs a value")?;
            (read_json(p)?, Path::new(p).parent().map(|d| d.to_path_buf()).unwrap_or_default())
        }
        None => (json!({}), Default::default()),
    };
    let o = Opts { argv, cfg: cfgv, dir };
    let arg = |k: &str| o.get(k);

        let (rs, var) = load_blocks(&o.path("--de").ok_or("--de")?)?;
        let e: usize = arg("--e").ok_or("--e")?.parse().map_err(|_| "--e")?;
        let sys = EBlockSystem::from_rules_with(&rs, e, &|c: &str| parse_coefficient(c, &var))?;
        let g = match o.cfg.get("graph") { Some(v) if v.is_object() && !argv.iter().any(|x| x == "--graph") => graph_from(v)?, _ => graph_from(&read_json(&o.path("--graph").ok_or("--graph")?)?)? };
        let num = |k: &str, d: &str| arg(k).unwrap_or(d.into());
        let mut cfg = PConfig::new(num("--w", "60").parse().unwrap(), num("--wf", "80").parse().unwrap(), num("--m", "80").parse().unwrap(), num("--loops", &g.loops.to_string()).parse().unwrap());
        cfg.s0 = num("--s0", "1/4").parse::<Rational>().map_err(|e| e.to_string())?;
        cfg.s_end = num("--S", "4").parse::<Rational>().map_err(|e| e.to_string())?;
        cfg.r = num("--r", "1/4").parse::<Rational>().map_err(|e| e.to_string())?;
        cfg.threads = num("--threads", "1").parse().unwrap();
        cfg.moments = num("--moments", "1").split(',').map(|x| x.parse().unwrap()).collect();
        cfg.nlo = num("--nlo", "-10").parse().unwrap();
        if let Some(w) = o.path("--workdir") {
            std::fs::create_dir_all(&w).map_err(|e| format!("{w}: {e}"))?;
            cfg.workdir = Some(w);
        }
        if let Some(f) = o.path("--provided") {
            // lower-loop massless values (E-normalisation, P^2 = 1) for graphs the evaluator cannot reduce
            let v = read_json(&f)?;
            cfg.provided = Provided::from_json(&v, crate::dispersive::run::bits(cfg.wf))?;
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
        if let Some(ex) = arg("--extra") {
            let efree: std::collections::HashMap<Vec<i64>, rug::Float> = match o.path("--efree") {
                Some(f) => {
                    let v: Value = serde_json::from_str(&std::fs::read_to_string(f).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
                    v.as_object().unwrap().iter().map(|(k, x)| (serde_json::from_str::<Vec<i64>>(k).unwrap(), rug::Float::with_val(crate::dispersive::run::bits(cfg.wf), rug::Float::parse(x.as_str().unwrap()).unwrap()))).collect()
                }
                None => Default::default(),
            };
            for t in ex.split(';') {
                let xv: Vec<i64> = t.split(',').map(|y| y.trim().parse().unwrap()).collect();
                let prec = crate::dispersive::run::bits(cfg.wf);
                let r1 = super::engine::extra_value(&rs.rules, &sys, &xv, &res.moments, &efree, &var, prec);
                let dig = match (&r1, &res2) {
                    (Ok((v, _)), Some(r2)) => super::engine::extra_value(&rs.rules, &sys, &xv, &r2.moments, &efree, &var, prec).ok().map(|(w, _)| {
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
}
