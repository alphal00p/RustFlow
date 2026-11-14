//! Command-line driver of one dispersive sample, shared by `examples/dispersive.rs` (RustFlow, harness reader) and the dev tools.
//!
//! dispersive --de DE.json --class-config CONFIG.json --m0 M0 --minf MINF --w W --wf WF --workdir DIR --out OUT.json
//!            --v read:V.json | native:PLAN.json,SOFTSOL.json,SAMPLE.json
//!            [--pr-values LANE_OUT.json] [--pivoting MODE] [--threads N] [--loops L] [--two-precision dW,dWF]
//!            [--code-version STR] [--ref ORACLE_OUT.json]... [--kernel patch|standin] [--stop-at-patch K]
//!            [--checkpoint-every-seconds T] [--frobenius-checkpoint-every-seconds T]
//! The transport checkpoints after EVERY patch (exact, ~0.1-1 MB); the zero and infinity families every T seconds (default 900)
//! between two Frobenius orders (exact binary, ~1 GB at 5 loops); a restart resumes from the last one (a shed allows 2 s: no flush).
//! dispersive --progress WORKDIR: one JSON progress line from the checkpoint files (see `progress`).
//! --kernel patch (default) = RustFlow's fixed-precision patch kernel (checkpointed after every patch); standin = the W3a
//! reference recurrence (tests only, no per-patch checkpoint). --two-precision dW,dWF may be negative (a cheaper check run).
//! --stop-at-patch K stops the transport resumably at patch K (short-segment runs; the same command without it continues).
//! The e slot comes from the class config ("e"); the loop number from the config ("loops") or the class family file's loop momenta.
//! --pr-values replaces the 39 'PR'-lane endpoint values of SAMPLE.json by those of a RustFlow lane run (raw = J / e^(4 gamma eps)).
use super::boundary;
use super::run::{self, Config};
use super::system::EBlockSystem;
use super::transport::{TaylorRecurrence, Transport};
use rug::Float;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::time::Instant;

pub struct Args {
    pub de: String,
    pub class_config: String,
    pub v: String,
    pub m0: i64,
    pub minf: i64,
    pub w: u32,
    pub wf: u32,
    pub workdir: String,
    pub out: String,
    pub pr_values: Option<String>,
    pub pivoting: String,
    pub threads: usize,
    pub loops: Option<u32>,
    pub two_precision: Option<(i64, i64)>,
    pub kernel: String,
    pub stop_at_patch: Option<usize>,
    pub checkpoint_every_seconds: f64,
    pub frobenius_checkpoint_every_seconds: f64,
    pub code_version: String,
    pub refs: Vec<String>,
}

pub fn rd(p: &str) -> Result<Value, String> {
    serde_json::from_slice(&std::fs::read(p).map_err(|e| format!("{p}: {e}"))?).map_err(|e| format!("{p}: {e}"))
}

impl Args {
    pub fn parse(argv: &[String]) -> Result<Args, String> {
        let mut m: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut i = 0;
        while i < argv.len() {
            let k = argv[i].strip_prefix("--").ok_or_else(|| format!("unexpected argument {}", argv[i]))?.to_string();
            let v = argv.get(i + 1).ok_or_else(|| format!("--{k} needs a value"))?.clone();
            m.entry(k).or_default().push(v);
            i += 2;
        }
        let one = |k: &str| -> Result<String, String> { m.get(k).and_then(|v| v.last().cloned()).ok_or_else(|| format!("--{k} is required")) };
        let opt = |k: &str| m.get(k).and_then(|v| v.last().cloned());
        let num = |k: &str| -> Result<i64, String> { one(k)?.parse().map_err(|_| format!("--{k}: not a number")) };
        Ok(Args {
            de: one("de")?,
            class_config: one("class-config")?,
            v: one("v")?,
            m0: num("m0")?,
            minf: num("minf")?,
            w: num("w")? as u32,
            wf: num("wf")? as u32,
            workdir: one("workdir")?,
            out: one("out")?,
            pr_values: opt("pr-values"),
            pivoting: opt("pivoting").unwrap_or_else(|| super::frobenius::DEFAULT_PIVOTING.into()),
            threads: opt("threads").map(|x| x.parse().unwrap_or(1)).unwrap_or(1),
            loops: opt("loops").and_then(|x| x.parse().ok()),
            two_precision: opt("two-precision").map(|x| {
                let (a, b) = x.split_once(',').expect("--two-precision dW,dWF");
                (a.parse().unwrap(), b.parse().unwrap())
            }),
            code_version: opt("code-version").unwrap_or_else(|| "unknown".into()),
            kernel: opt("kernel").unwrap_or_else(|| "patch".into()),
            stop_at_patch: opt("stop-at-patch").map(|x| x.parse().map_err(|_| "--stop-at-patch: not a number".to_string())).transpose()?,
            checkpoint_every_seconds: opt("checkpoint-every-seconds").map(|x| x.parse().unwrap_or(600.0)).unwrap_or(600.0),
            frobenius_checkpoint_every_seconds: opt("frobenius-checkpoint-every-seconds").map(|x| x.parse().unwrap_or(900.0)).unwrap_or(900.0),
            refs: m.get("ref").cloned().unwrap_or_default(),
        })
    }

    /// (e slot, loops) from the class config
    pub fn class_settings(&self) -> Result<(usize, u32, Value), String> {
        let cfg = rd(&self.class_config)?;
        let e = cfg["e"].as_u64().ok_or("class config has no e slot")? as usize;
        let loops = match (self.loops, cfg["loops"].as_u64()) {
            (Some(l), _) => l,
            (None, Some(l)) => l as u32,
            (None, None) => {
                let dir = std::path::Path::new(&self.class_config).parent().unwrap().to_path_buf();
                let fam = dir.join(cfg["family"].as_str().unwrap_or("family-eta-s0.json"));
                rd(fam.to_str().unwrap())?["loop_momenta"].as_array().ok_or("class family has no loop_momenta")?.len() as u32
            }
        };
        Ok((e, loops, cfg))
    }
}

/// endpoint list of a 4L sample, optionally with the 'PR' rows replaced by a RustFlow lane run
pub fn endpoints(sample: &Value, pr_values: Option<&Value>, prec: u32) -> Result<(Vec<(String, Vec<i64>, Float)>, Value), String> {
    let vi = |x: &Value| -> Vec<i64> { x.as_array().unwrap().iter().map(|v| v.as_i64().unwrap()).collect() };
    let mut lane: Option<BTreeMap<Vec<i64>, Float>> = None;
    let mut info = json!({"pr_source": "sample file (in-house dispersive lane)"});
    if let Some(o) = pr_values {
        // raw = J / e^(4 gamma eps) = -(sin(pi d0/2)/pi) I / e^(4 gamma eps)  (values/d5l/v9a/tools/pr_samples.py)
        let d0: rug::Rational = o["d0"].as_str().ok_or("lane out has no d0")?.parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
        let pre = super::dispersion::prefactor(&d0, prec);
        let eps = Float::with_val(prec, (Float::with_val(prec, 4) - Float::with_val(prec, &d0)) / 2u32);
        let gam = Float::with_val(prec, Float::with_val(prec, rug::float::Constant::Euler) * &eps * 4u32).exp();
        let mut mp = BTreeMap::new();
        for (m, iv) in o["masters"].as_array().ok_or("lane out has no masters")?.iter().zip(o["I"].as_array().unwrap()) {
            let i = Float::with_val(prec, Float::parse(iv["I"].as_str().unwrap()).map_err(|e| e.to_string())?);
            mp.insert(vi(m), Float::with_val(prec, &pre * &i) / &gam);
        }
        info = json!({"pr_source": "RustFlow dispersive lane run", "lane_out_fingerprint": o["fingerprint"], "lane_provenance": o["provenance"]});
        lane = Some(mp);
    }
    if let Some(o) = pr_values {
        if o["eps"].as_str() != sample["epsilon"].as_str() {
            return Err(format!("lane run (eps {}) and 4L sample (eps {}) are at different epsilon", o["eps"], sample["epsilon"]));
        }
    }
    let mut out = Vec::new();
    let mut replaced = 0;
    for x in sample["endpoint"].as_array().ok_or("sample has no endpoint list")? {
        let t = vi(&x["target"]);
        let fam = x.get("family").and_then(|f| f.as_str()).unwrap_or("");
        let key = format!("{fam}{}", serde_json::to_string(&t).unwrap().replace(',', ", "));
        let mut v = Float::with_val(prec, Float::parse(x["value"]["re"].as_str().unwrap()).map_err(|e| e.to_string())?);
        if fam == "PR" {
            if let Some(l) = &lane {
                v = l.get(&t).ok_or_else(|| format!("lane run has no master {t:?}"))?.clone();
                replaced += 1;
            }
        }
        out.push((key, t, v));
    }
    if lane.is_some() && replaced == 0 {
        return Err("--pr-values given but the sample has no PR rows".into());
    }
    info["pr_rows_replaced"] = json!(replaced);
    Ok((out, info))
}

/// boundary rows + provenance for the V spec
pub fn boundary_rows(args: &Args, sys: &EBlockSystem, loops: u32, prec: u32) -> Result<(Vec<boundary::BoundaryRow>, Value), String> {
    let (v, src) = if let Some(spec) = args.v.strip_prefix("native:") {
        let f: Vec<&str> = spec.split(',').collect();
        if f.len() != 3 {
            return Err("--v native:PLAN,SOFTSOL,SAMPLE".into());
        }
        let pl = rd(f[0])?;
        let pl = if pl.get("plan").is_some() { pl["plan"].clone() } else { pl };
        let lane = args.pr_values.as_deref().map(rd).transpose()?;
        let (ep, info) = endpoints(&rd(f[2])?, lane.as_ref(), prec + 64)?;
        let (v, d) = boundary::solve_v(&pl, &rd(f[1])?, &ep, sys.e_slot, prec + 64)?;
        (v, json!({"V": "native (soft_gen solve port)", "plan": f[0], "soft_solution": f[1], "sample": f[2], "endpoints": info, "diag": d}))
    } else if let Some(p) = args.v.strip_prefix("read:") {
        let (v, m) = boundary::read_v_file(p)?;
        (v, json!({"V": "READ from a V file (see meta.producer; in-house soft_gen files have none)", "file": p, "meta": m}))
    } else {
        return Err("--v read:V.json | native:PLAN,SOFTSOL,SAMPLE".into());
    };
    Ok((boundary::rows_from_v(sys, &v, loops, prec)?, src))
}

fn digits_vs(out: &Value, oracle: &Value) -> Vec<Option<f64>> {
    // compare by master exponent vector when the oracle lists masters, else by index (the e-block order is asserted equal to the DE's)
    let om: Option<Vec<Value>> = oracle["masters"].as_array().cloned();
    let mine = out["masters"].as_array().unwrap();
    mine.iter()
        .enumerate()
        .map(|(i, m)| {
            let j = match &om {
                Some(list) => list.iter().position(|x| x == m)?,
                None => i,
            };
            let x = Float::with_val(4000, Float::parse(out["I"][i]["I"].as_str()?).ok()?);
            let y = Float::with_val(4000, Float::parse(oracle["I"][j]["I"].as_str()?).ok()?);
            if y.is_zero() {
                return Some(if x.is_zero() { 99.0 } else { -x.abs().log10().to_f64() });
            }
            let d = Float::with_val(4000, (x - &y) / y).abs();
            Some(if d.is_zero() { 99.0 } else { -d.log10().to_f64() })
        })
        .collect()
}

/// run (one or two precisions), export, compare. `cancel` / `transport_active`: see [`Config`] (the driver's signal handling)
pub fn run(args: &Args, sys: &EBlockSystem) -> Result<Value, String> {
    run_with(args, sys, &Default::default(), &Default::default())
}

pub fn run_with(
    args: &Args,
    sys: &EBlockSystem,
    cancel: &crate::options::CancellationToken,
    transport_active: &std::sync::Arc<std::sync::atomic::AtomicBool>,
) -> Result<Value, String> {
    let t0 = Instant::now();
    let (e_slot, loops, class_cfg) = args.class_settings()?;
    if e_slot != sys.e_slot {
        return Err(format!("system built with e slot {} but the class config says {e_slot}", sys.e_slot));
    }
    let mk = |w: u32, wf: u32, dir: String| -> Config {
        let mut c = Config::new(w, wf, args.m0, args.minf, loops);
        c.nthreads = args.threads;
        c.workdir = dir;
        c.pivoting = args.pivoting.clone();
        c.code_version = args.code_version.clone();
        c.cancel = cancel.clone();
        c.transport_active = transport_active.clone();
        c.stop_at_patch = args.stop_at_patch;
        c.checkpoint_every_seconds = args.checkpoint_every_seconds;
        c.frobenius_checkpoint_every_seconds = args.frobenius_checkpoint_every_seconds;
        c
    };
    let kernel: Box<dyn Transport> = match args.kernel.as_str() {
        "patch" => Box::new(super::kernel::PatchKernel::new(sys, args.threads)?),
        "standin" => Box::new(TaylorRecurrence::new(sys, args.threads)),
        k => return Err(format!("--kernel {k}: patch | standin")),
    };
    let mut runs = vec![(args.w, args.wf)];
    if let Some((dw, dwf)) = args.two_precision {
        let add = |a: u32, d: i64| -> Result<u32, String> { u32::try_from(a as i64 + d).map_err(|_| "--two-precision: negative precision".to_string()) };
        runs.push((add(args.w, dw)?, add(args.wf, dwf)?));
    }
    let mut outs = Vec::new();
    for (k, &(w, wf)) in runs.iter().enumerate() {
        let dir = if runs.len() == 1 { args.workdir.clone() } else { format!("{}/p{}", args.workdir, k + 1) };
        std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
        let cfg = mk(w, wf, dir);
        let meter = run::Meter::start("plan (singular points, family structure, path) + boundary rows");
        let plan = run::plan(sys, &cfg)?;
        let (bnd, src) = boundary_rows(args, sys, loops, run::bits(wf))?;
        meter.finish(&cfg.workdir, cfg.nthreads, true);
        let prov = json!({"de": args.de, "class_config": args.class_config, "boundary": src, "transport": kernel.id()});
        outs.push(run::run_sample(sys, &cfg, &plan, &bnd, kernel.as_ref(), prov)?);
    }
    let mut out = outs[0].clone();
    if outs.len() == 2 {
        let dg = run::two_precision_digits(&outs[0], &outs[1]);
        for (i, d) in dg.iter().enumerate() {
            out["I"][i]["digits"] = json!((d * 10.0).round() / 10.0);
        }
        out["digits_source"] = json!(format!("two-precision: (W, WF) = ({}, {}) vs ({}, {}); 'digits' = agreement of the two runs", runs[0].0, runs[0].1, runs[1].0, runs[1].1));
        out["second_run"] = json!({"w": runs[1].0, "wf": runs[1].1, "boundary_fit": outs[1]["boundary_fit"], "matching": outs[1]["matching"],
            "residuals": outs[1]["residuals"], "transport": outs[1]["transport"], "cost": outs[1]["cost"], "fingerprint": outs[1]["fingerprint"]});
        let mut s: Vec<f64> = dg.clone();
        s.sort_by(|a, b| a.total_cmp(b));
        out["digits_summary"] = json!({"min": s.first(), "median": s.get(s.len() / 2), "named_idx_564_568": (564..569).filter_map(|i| dg.get(i).map(|d| (d * 10.0).round() / 10.0)).collect::<Vec<_>>()});
    }
    // named masters (class config "names": {name: exponent vector})
    if let Some(names) = class_cfg["names"].as_object() {
        let mut named = serde_json::Map::new();
        for (n, m) in names {
            if let Some(i) = sys.masters.iter().position(|x| json!(x) == *m) {
                named.insert(n.clone(), json!({"idx": i, "I": out["I"][i]["I"], "J": out["J"][i], "digits": out["I"][i].get("digits"), "digits_estimate": out["I"][i]["digits_estimate"]}));
            } else {
                return Err(format!("named master {n} {m} is not in the e-block"));
            }
        }
        out["named"] = Value::Object(named);
    }
    let mut cmp = serde_json::Map::new();
    for r in &args.refs {
        let o = rd(r)?;
        let dg = digits_vs(&out, &o);
        let all: Vec<f64> = dg.iter().flatten().cloned().collect();
        let mut s = all.clone();
        s.sort_by(|a, b| a.total_cmp(b));
        let named: serde_json::Map<String, Value> = out
            .get("named")
            .and_then(|n| n.as_object())
            .map(|n| n.iter().map(|(k, v)| (k.clone(), json!(dg[v["idx"].as_u64().unwrap() as usize].map(|x| (x * 10.0).round() / 10.0)))).collect())
            .unwrap_or_default();
        let per: Vec<Value> = dg.iter().map(|d| d.map(|x| json!((x * 10.0).round() / 10.0)).unwrap_or(Value::Null)).collect();
        cmp.insert(r.clone(), json!({"entries": all.len(), "min": s.first(), "median": s.get(s.len() / 2), "named": named, "per_entry": per}));
    }
    out["oracle_comparison"] = Value::Object(cmp);
    out["wall_seconds"] = json!(t0.elapsed().as_secs_f64());
    let tmp = format!("{}.part", args.out);
    std::fs::write(&tmp, serde_json::to_string_pretty(&out).unwrap()).map_err(|e| e.to_string())?;
    std::fs::rename(&tmp, &args.out).map_err(|e| e.to_string())?;
    Ok(out)
}

/// Progress of a sample from its CHECKPOINT FILES only (so it can never run ahead of what a restart resumes): one JSON object
/// {"stage", "patch", "npatches", "detail", "ckpt_utc", "units_note"}. For drivers that compare (stage, patch) between attempts the
/// stage is the constant "transport" until every run is final ("done"), and "patch" is a monotone work counter over all runs of the
/// workdir: per run 150 units for the zero family (scaled by its Frobenius checkpoint), 1 per transport patch, 60 per infinity family
/// (scaled likewise), 5 for matching/assembly (in-house cost ratios at W165). `{}` when nothing is checkpointed yet.
pub fn progress(workdir: &str) -> Value {
    use std::path::Path;
    let dirs: Vec<String> = if Path::new(&format!("{workdir}/p1")).is_dir() {
        (1..=9).map(|k| format!("{workdir}/p{k}")).filter(|d| Path::new(d).is_dir()).collect()
    } else {
        vec![workdir.to_string()]
    };
    let nruns = if Path::new(&format!("{workdir}/p1")).is_dir() { dirs.len().max(2) } else { 1 };
    let mut newest = 0u64;
    let mut seen = |p: &str| -> bool {
        match std::fs::metadata(p) {
            Ok(m) => {
                if let Ok(t) = m.modified() {
                    newest = newest.max(t.duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0));
                }
                true
            }
            Err(_) => false,
        }
    };
    let (mut units, mut total) = (0.0f64, 0.0f64);
    let mut details = Vec::new();
    let mut all_final = true;
    let mut last_total = None;
    let mut any = false;
    for (k, d) in dirs.iter().enumerate() {
        let sum: Option<Value> = std::fs::read(format!("{d}/plan-summary.json")).ok().and_then(|b| serde_json::from_slice(&b).ok());
        let (np, loops) = match &sum {
            Some(s) => (s["npatches"].as_u64().unwrap_or(0) as f64, s["loops"].as_u64().unwrap_or(5) as f64),
            None => (0.0, 5.0),
        };
        let run_total = if sum.is_some() { 150.0 + np + 60.0 * loops + 5.0 } else { last_total.unwrap_or(0.0) };
        last_total = Some(run_total);
        total += run_total;
        let frac = |p: String| -> f64 {
            match super::ckpt::peek(Path::new(&p)) {
                Some((_, a, b, m)) if b > a => ((m - a) as f64 / (b - a) as f64).clamp(0.0, 1.0),
                _ => 0.0,
            }
        };
        let mut u = 0.0;
        let mut what = "not started".to_string();
        if seen(&format!("{d}/stage-zero.json")) {
            u += 150.0;
            what = "zero family done".into();
        } else if seen(&format!("{d}/frob-zero.bin")) {
            u += 142.0 * frac(format!("{d}/frob-zero.bin"));
            what = format!("zero family {:.0} %", 100.0 * frac(format!("{d}/frob-zero.bin")));
        }
        if seen(&format!("{d}/stage-transport.json")) {
            u += np;
            what = "transport done".into();
        } else if seen(&format!("{d}/stage-transport-kernel.json")) {
            let nx = std::fs::read(format!("{d}/stage-transport-kernel.json")).ok().and_then(|b| serde_json::from_slice::<Value>(&b).ok()).and_then(|v| v["next_patch"].as_u64()).unwrap_or(0);
            u += nx as f64;
            what = format!("transport patch {nx} of {np}");
        }
        for kk in 1..=(loops as usize) {
            if seen(&format!("{d}/stage-inf-{kk}.json")) {
                u += 60.0;
                what = format!("infinity family k={kk} done");
            } else if seen(&format!("{d}/frob-inf-{kk}.bin")) {
                let f = frac(format!("{d}/frob-inf-{kk}.bin"));
                u += 57.0 * f;
                what = format!("infinity family k={kk} {:.0} %", 100.0 * f);
            }
        }
        let fin = seen(&format!("{d}/stage-final.json"));
        if fin {
            u += 5.0;
            what = "final".into();
        }
        all_final &= fin;
        any |= u > 0.0 || sum.is_some();
        units += u;
        details.push(format!("run {} of {nruns}: {what}", k + 1));
    }
    if dirs.len() < nruns {
        all_final = false;
        total += last_total.unwrap_or(0.0) * (nruns - dirs.len()) as f64;
    }
    if !any {
        return json!({});
    }
    let utc = chrono_like(newest);
    json!({"stage": if all_final { "done" } else { "transport" }, "patch": units.floor() as u64, "npatches": total.ceil() as u64,
           "detail": details.join("; "), "ckpt_utc": utc,
           "units_note": "patch/npatches = monotone work units over all runs (150 zero family, 1 per patch, 60 per infinity family, 5 matching), read from the checkpoint files"})
}

fn chrono_like(secs: u64) -> String {
    // UTC timestamp without a date crate (civil-from-days, H. Hinnant)
    let (d, s) = ((secs / 86400) as i64, secs % 86400);
    let z = d + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + if month <= 2 { 1 } else { 0 };
    format!("{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z", s / 3600, (s / 60) % 60, s % 60)
}
