//! Propagator configuration of the RustFlow dispersive driver (route B; docs/pflow.md). One entry point for both front ends:
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

/// process exit code for a finished run: 0 if accepted (or --accept-diagnostics), 3 if an acceptance gate failed
pub fn exit_code(v: &Value) -> i32 {
    if v["accepted"] == true || v["accept_diagnostics"] == true { 0 } else { 3 }
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
        let mut rs = RuleSet::from_json(&de)?;
        // RustRed writes an identity rule (master = 1 x master) when a master is also a target: not a reduction
        let own: Vec<Vec<i64>> = rs.masters.clone();
        rs.rules.retain(|t, r| !(own.contains(t) || (r.len() == 1 && r[0].0 == *t && r[0].1.trim() == "1")));
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
    let pnum = |k: &str, d: &str| -> Result<i64, String> { arg(k).unwrap_or(d.into()).parse::<i64>().map_err(|_| format!("{k}: not an integer")) };
    let prat = |k: &str, d: &str| -> Result<Rational, String> { arg(k).unwrap_or(d.into()).parse::<Rational>().map_err(|_| format!("{k}: not a rational")) };
    let e = pnum("--e", "-1")?;
    if e < 0 {
        return Err("--e (massive line slot) is required".into());
    }
    let e = e as usize;
    let sys = EBlockSystem::from_rules_with(&rs, e, &|c: &str| parse_coefficient(c, &var))?;
    let g = match o.cfg.get("graph") { Some(v) if v.is_object() && !argv.iter().any(|x| x == "--graph") => graph_from(v)?, _ => graph_from(&read_json(&o.path("--graph").ok_or("--graph")?)?)? };
    let (w, wf, m) = (pnum("--w", "60")?, pnum("--wf", "80")?, pnum("--m", "80")?);
    if w < 10 || wf < w || m < 10 {
        return Err(format!("need w >= 10, wf >= w, m >= 10 (w {w}, wf {wf}, m {m})"));
    }
    let mut cfg = PConfig::new(w as u32, wf as u32, m, pnum("--loops", &g.loops.to_string())? as u32);
    cfg.s0 = prat("--s0", "1/4")?;
    cfg.s_end = prat("--S", "4")?;
    cfg.r = prat("--r", "1/4")?;
    cfg.threads = pnum("--threads", "1")?.max(1) as usize;
    cfg.moments = arg("--moments").unwrap_or("1".into()).split(',').map(|x| x.trim().parse::<i32>().map_err(|_| format!("--moments: {x}"))).collect::<Result<_, _>>()?;
    cfg.nlo = pnum("--nlo", "-10")?;
    cfg.gate_margin = pnum("--gate-margin", "20")? as u32;
    cfg.allow_skipped = argv.iter().any(|x| x == "--allow-skipped") || o.cfg["allow_skipped"] == true;
    let accept_diag = argv.iter().any(|x| x == "--accept-diagnostics") || o.cfg["accept_diagnostics"] == true;
    let workdir = o.path("--workdir");
    if let Some(f) = o.path("--provided") {
        let v = read_json(&f)?;
        cfg.provided = Provided::from_json(&v, crate::dispersive::run::bits(cfg.wf), &rs.d0).map_err(|e| format!("{f}: {e}"))?;
    }
    // one checkpoint subdirectory per precision pass (the passes have different fingerprints)
    let pass_dir = |c: &PConfig| -> Result<Option<String>, String> {
        match &workdir {
            None => Ok(None),
            Some(w) => {
                let d = format!("{w}/p{}-{}-{}", c.w, c.wf, c.m0);
                std::fs::create_dir_all(&d).map_err(|e| format!("{d}: {e}"))?;
                Ok(Some(d))
            }
        }
    };
    cfg.workdir = pass_dir(&cfg)?;
    let tp = pnum("--two-precision", "0")?;
    if tp < 0 || tp as u32 + 10 > cfg.w || tp > cfg.m0 - 10 {
        return Err(format!("--two-precision {tp}: must be >= 0 and leave w >= 10 and m >= 10"));
    }
    let res = run_b(&sys, &g, e, &cfg)?;
    let res2 = if tp > 0 {
        let mut c2 = cfg.clone();
        c2.w -= tp as u32;
        c2.wf -= tp as u32;
        c2.m0 -= tp;
        c2.minf -= tp;
        c2.mthr -= tp;
        c2.workdir = pass_dir(&c2)?;
        Some(run_b(&sys, &g, e, &c2)?)
    } else {
        None
    };
    let accepted = res.accepted && res2.as_ref().map_or(true, |r| r.accepted);
    // what each moment is (module docs of engine.rs): m = 1 massless value; m = n + 1 >= 2: (a_e)_n / n! E[a + n 1_e];
    // m = -k <= 0: (-1)^(a_e - 1) (a_e)_n / n! E[a_e -> -n], n = k + 1 - a_e (zero if n < 0)
    let poch = |a: i64, n: i64| -> Rational { let mut r = Rational::from(1); for i in 0..n { r *= Rational::from(a + i); r /= Rational::from(i + 1); } r };
    let mut out = Vec::new();
    for (k, (m, v)) in res.moments.iter().enumerate() {
        for (j, mast) in sys.masters.iter().enumerate() {
            let sgn = if mast.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
            let prec = v[j].prec();
            let val = rug::Float::with_val(prec, &v[j] * sgn);
            let dig = res2.as_ref().map(|r2| {
                let w = &r2.moments[k].1[j];
                let d = rug::Float::with_val(prec, &v[j] - w);
                if d.is_zero() { 999.0 } else { -(d.abs() / rug::Float::with_val(prec, v[j].abs_ref())).log10().to_f64() }
            });
            let ndig = (cfg.wf as usize).max(20);
            let mut row = json!({"m": m, "master": mast, "value": val.to_string_radix(10, Some(ndig)), "digits_two_precision": dig});
            let ae = mast[e];
            if *m >= 2 {
                let n = (*m - 1) as i64;
                let mut idx = mast.clone();
                idx[e] += n;
                let f = poch(ae, n);
                row["integral"] = json!(idx);
                row["integral_value"] = json!(rug::Float::with_val(prec, &val / &f).to_string_radix(10, Some(ndig)));
                row["quantity"] = json!(format!("value = (a_e)_n/n! E[integral], a_e = {ae}, n = {n}"));
            } else if *m <= 0 {
                let n = (1 - *m) as i64 - ae;
                if n >= 0 {
                    let mut idx = mast.clone();
                    idx[e] = -n;
                    let f = Rational::from(&poch(ae, n) * if (ae - 1) % 2 != 0 { -1 } else { 1 });
                    row["integral"] = json!(idx);
                    row["integral_value"] = json!(rug::Float::with_val(prec, &val / &f).to_string_radix(10, Some(ndig)));
                    row["quantity"] = json!(format!("value = (-1)^(a_e-1) (a_e)_n/n! E[integral] (line e as a numerator of power n), a_e = {ae}, n = {n}"));
                } else {
                    row["quantity"] = json!("no integer term at this order: the value should vanish");
                }
            }
            out.push(row);
        }
    }
    // non-master targets: --extra "a1,a2,...;b1,..." (pinched-master values from --efree FILE {"[..]": "value"} if needed)
    let mut extras = Vec::new();
    if let Some(ex) = arg("--extra") {
        let prec = crate::dispersive::run::bits(cfg.wf);
        let efree: std::collections::HashMap<Vec<i64>, rug::Float> = match o.path("--efree") {
            Some(f) => {
                let v = read_json(&f)?;
                let mut m = std::collections::HashMap::new();
                for (k, x) in v.as_object().ok_or("--efree: an object {\"[indices]\": \"value\"}")? {
                    let idx: Vec<i64> = serde_json::from_str(k).map_err(|e| format!("--efree key {k}: {e}"))?;
                    let val = rug::Float::parse(x.as_str().ok_or("--efree: values are strings")?).map_err(|e| format!("--efree {k}: {e}"))?;
                    m.insert(idx, rug::Float::with_val(prec, val));
                }
                m
            }
            None => Default::default(),
        };
        for t in ex.split(';') {
            let xv: Vec<i64> = t.split(',').map(|y| y.trim().parse::<i64>().map_err(|_| format!("--extra: {t}"))).collect::<Result<_, _>>()?;
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
    Ok(json!({"schema": "pflow-route-b.v1", "accepted": accepted, "accept_diagnostics": accept_diag, "digits_cap": res.diag["digits_cap"], "extras": extras, "d0": rs.d0.to_string(), "e_slot": e, "config": {"w": cfg.w, "wf": cfg.wf, "m": cfg.m0, "s0": cfg.s0.to_string(), "S": cfg.s_end.to_string(), "r": cfg.r.to_string(), "moments": cfg.moments, "nlo": cfg.nlo, "gate_margin": cfg.gate_margin},
            "values": out, "diag": res.diag, "diag_second": res2.map(|r| r.diag)}))
}
