//! Command-line driver of one dispersive sample, shared by `examples/dispersive.rs` (RustFlow, harness reader) and the dev tools.
//!
//! dispersive --de DE.json --class-config CONFIG.json --m0 M0 --minf MINF --w W --wf WF --workdir DIR --out OUT.json
//!            --v read:V.json | native:PLAN.json,SOFTSOL.json,SAMPLE.json | lane:CHAIN.json,SAMPLE.json
//!            [--pr-values LANE_OUT.json] [--pivoting MODE] [--threads N] [--loops L] [--two-precision dW,dWF]
//!            [--code-version STR] [--ref ORACLE_OUT.json]...
//! The e slot comes from the class config ("e"); the loop number from the config ("loops") or the class family file's loop momenta.
//! --pr-values replaces the 39 'PR'-lane endpoint values of SAMPLE.json by those of a RustFlow lane run (raw = J / e^(4 gamma eps)).
//! lane:CHAIN.json,SAMPLE.json builds the soft V of a LANE family (e.g. the 4L prism) natively from the P14 rows of SAMPLE.json
//! (lower.rs; chain spec documented at lower::lane_boundary).
//! --merge TEMPLATE.json,P14SAMPLE.json,OUT.json (lane runs): writes OUT.json = TEMPLATE (a merged 4L sample) with its 'PR' rows
//! replaced by this run's values (raw = J / e^(4 gamma eps)) and a 'provenance' field on every row (merged_sample).
//! --input-budget BUDGET.json (with --merge): per-row relative error budget {"source": .., "entries": [{"family": "PR" | "", "target": [..],
//! "rel_error": "1e-90"}]} (e.g. Jacobian x input accuracy). Row provenance then gets 'digits' = -log10(rel_error + 10^-two_precision)
//! for lane rows (-log10(rel_error) for P14 rows), 'digits_input_budget' and 'digits_two_precision' separately. Without it 'digits'
//! is the two-precision agreement only, and 'digits_kind' says so.
//! --p14-producer STR (with --merge): names the producer of the P14 rows (e.g. the RustFlow commit) in their provenance.
use super::boundary;
use super::run::{self, Config};
use super::system::EBlockSystem;
use super::transport::TaylorRecurrence;
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
    pub two_precision: Option<(u32, u32)>,
    pub code_version: String,
    pub refs: Vec<String>,
    pub merge: Option<String>,
    pub input_budget: Option<String>,
    pub p14_producer: Option<String>,
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
            refs: m.get("ref").cloned().unwrap_or_default(),
            merge: opt("merge"),
            input_budget: opt("input-budget"),
            p14_producer: opt("p14-producer"),
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
    } else if let Some(spec) = args.v.strip_prefix("lane:") {
        let f: Vec<&str> = spec.split(',').collect();
        if f.len() != 2 {
            return Err("--v lane:CHAIN.json,SAMPLE.json".into());
        }
        let (_, _, cfg) = args.class_settings()?;
        let (v, prov) = super::lower::lane_boundary(f[0], f[1], &cfg, &sys.masters, &sys.d0, prec + 64)?;
        let worst = prov["lower_diag"]["endpoint_system"]["endpoint_consistency_worst_rel_log10"].as_f64().unwrap_or(0.0);
        if worst > -20.0 {
            return Err(format!("lower endpoint system inconsistent: worst relative residual 1e{worst:.1}"));
        }
        (v, prov)
    } else if let Some(p) = args.v.strip_prefix("read:") {
        let (v, m) = boundary::read_v_file(p)?;
        (v, json!({"V": "READ from a V file (see meta.producer; in-house soft_gen files have none)", "file": p, "meta": m}))
    } else {
        return Err("--v read:V.json | native:PLAN,SOFTSOL,SAMPLE | lane:CHAIN,SAMPLE".into());
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

/// run (one or two precisions), export, compare
pub fn run(args: &Args, sys: &EBlockSystem) -> Result<Value, String> {
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
        c
    };
    let kernel = TaylorRecurrence::new(sys, args.threads);
    let mut runs = vec![(args.w, args.wf)];
    if let Some((dw, dwf)) = args.two_precision {
        runs.push((args.w + dw, args.wf + dwf));
    }
    let mut outs = Vec::new();
    for (k, &(w, wf)) in runs.iter().enumerate() {
        let dir = if runs.len() == 1 { args.workdir.clone() } else { format!("{}/p{}", args.workdir, k + 1) };
        std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
        let cfg = mk(w, wf, dir);
        let plan = run::plan(sys, &cfg)?;
        let (bnd, src) = boundary_rows(args, sys, loops, run::bits(wf))?;
        let prov = json!({"de": args.de, "class_config": args.class_config, "boundary": src, "transport": "TaylorRecurrence stand-in (W3a)"});
        outs.push(run::run_sample(sys, &cfg, &plan, &bnd, &kernel, prov)?);
    }
    let mut out = outs[0].clone();
    if outs.len() == 2 {
        let dg = run::two_precision_digits(&outs[0], &outs[1]);
        for (i, d) in dg.iter().enumerate() {
            out["I"][i]["digits"] = json!((d * 10.0).round() / 10.0);
        }
        out["digits_source"] = json!(format!("two-precision: (W, WF) = ({}, {}) vs ({}, {}); 'digits' = agreement of the two runs", runs[0].0, runs[0].1, runs[1].0, runs[1].1));
        out["second_run"] = json!({"w": runs[1].0, "wf": runs[1].1, "boundary_fit": outs[1]["boundary_fit"], "matching": outs[1]["matching"], "fingerprint": outs[1]["fingerprint"]});
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
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let named: serde_json::Map<String, Value> = out
            .get("named")
            .and_then(|n| n.as_object())
            .map(|n| n.iter().map(|(k, v)| (k.clone(), json!(dg[v["idx"].as_u64().unwrap() as usize].map(|x| (x * 10.0).round() / 10.0)))).collect())
            .unwrap_or_default();
        cmp.insert(r.clone(), json!({"entries": all.len(), "min": s.first(), "median": s.get(s.len() / 2), "named": named}));
    }
    out["oracle_comparison"] = Value::Object(cmp);
    out["wall_seconds"] = json!(t0.elapsed().as_secs_f64());
    if let Some(m) = &args.merge {
        let f: Vec<&str> = m.split(',').collect();
        if f.len() != 3 {
            return Err("--merge TEMPLATE.json,P14SAMPLE.json,OUT.json".into());
        }
        let budget = args.input_budget.as_deref().map(rd).transpose()?;
        let merged = merged_sample(&rd(f[0])?, &rd(f[1])?, f[1], &out, &args.out, run::bits(args.w) + 64, budget.as_ref(), args.p14_producer.as_deref())?;
        let tmp = format!("{}.part", f[2]);
        std::fs::write(&tmp, serde_json::to_string(&merged).unwrap()).map_err(|e| e.to_string())?;
        std::fs::rename(&tmp, f[2]).map_err(|e| e.to_string())?;
        out["merged_sample"] = json!({"path": f[2], "summary": merged["provenance"]});
    }
    let tmp = format!("{}.part", args.out);
    std::fs::write(&tmp, serde_json::to_string_pretty(&out).unwrap()).map_err(|e| e.to_string())?;
    std::fs::rename(&tmp, &args.out).map_err(|e| e.to_string())?;
    Ok(out)
}

/// A merged 4L sample (schema of values/d5l/v9a/pr/samples: epsilon, attempt, index, lane_tag, endpoint[{family?, target, value{re}}])
/// whose 'PR' rows are this lane run's values, raw = J / e^(4 gamma eps) = -(sin(pi d0/2)/pi) I / e^(4 gamma eps) (pr_samples.py).
/// The P14 rows are taken from the RustFlow P14 sample (same targets and order as the template; the number of rows whose values differ
/// from the template is recorded). Every row gets 'provenance'.
pub fn merged_sample(template: &Value, p14: &Value, p14_path: &str, lane_out: &Value, lane_out_path: &str, prec: u32, budget: Option<&Value>, p14_producer: Option<&str>) -> Result<Value, String> {
    let p14_tool = match p14_producer {
        Some(p) => format!("RustFlow eta-flows + spectral gluing (P14 sample; {p})"),
        None => "RustFlow eta-flows + spectral gluing (P14 sample)".to_string(),
    };
    // relative error budget per (family, target)
    let mut bmap: BTreeMap<(String, Vec<i64>), f64> = BTreeMap::new();
    if let Some(b) = budget {
        for e in b["entries"].as_array().ok_or("input budget has no entries")? {
            let fam = e["family"].as_str().unwrap_or("").to_string();
            let t: Vec<i64> = e["target"].as_array().ok_or("budget entry without target")?.iter().map(|v| v.as_i64().unwrap()).collect();
            let r = e["rel_error"].as_str().ok_or("budget entry without rel_error")?;
            let x = Float::with_val(64, Float::parse(r).map_err(|e| e.to_string())?);
            bmap.insert((fam, t), if x.is_zero() { -1000.0 } else { x.abs().log10().to_f64() });
        }
    }
    if template["epsilon"] != p14["epsilon"] || lane_out["eps"] != template["epsilon"] {
        return Err(format!("epsilon mismatch: template {}, P14 sample {}, lane run {}", template["epsilon"], p14["epsilon"], lane_out["eps"]));
    }
    let d0: rug::Rational = lane_out["d0"].as_str().ok_or("lane out has no d0")?.parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
    let pre = super::dispersion::prefactor(&d0, prec);
    let eps = Float::with_val(prec, (Float::with_val(prec, 4) - Float::with_val(prec, &d0)) / 2u32);
    let gam = Float::with_val(prec, Float::with_val(prec, rug::float::Constant::Euler) * &eps * 4u32).exp();
    let vi = |x: &Value| -> Vec<i64> { x.as_array().unwrap().iter().map(|v| v.as_i64().unwrap()).collect() };
    let mut lane: BTreeMap<Vec<i64>, (String, Value)> = BTreeMap::new();
    for (m, iv) in lane_out["masters"].as_array().ok_or("lane out has no masters")?.iter().zip(lane_out["I"].as_array().ok_or("lane out has no I")?) {
        let s = iv["I"].as_str().ok_or("I entry without I")?;
        let i = Float::with_val(prec, Float::parse(s).map_err(|e| e.to_string())?);
        let raw = Float::with_val(prec, &pre * &i) / &gam;
        let nd = s.trim_start_matches('-').chars().filter(|c| c.is_ascii_digit()).count().max(20);
        lane.insert(vi(m), (raw.to_string_radix(10, Some(nd)), iv.get("digits").cloned().unwrap_or(Value::Null)));
    }
    let p14rows = p14["endpoint"].as_array().ok_or("P14 sample has no endpoints")?;
    let tool_lane = json!({"tool": "RustFlow dispersive mode (lane run)", "out": lane_out_path, "fingerprint": lane_out["fingerprint"],
        "code_version": lane_out["config"]["code_version"], "digits_source": lane_out["digits_source"]});
    let mut ep = Vec::new();
    let (mut n14, mut npr, mut n14_differ) = (0, 0, 0);
    for x in template["endpoint"].as_array().ok_or("template has no endpoints")? {
        let fam = x.get("family").and_then(|f| f.as_str()).unwrap_or("");
        let t = vi(&x["target"]);
        let mut row = x.clone();
        match fam {
            "" => {
                let y = p14rows.get(n14).ok_or("template has more P14 rows than the P14 sample")?;
                if y["target"] != x["target"] {
                    return Err(format!("P14 row {n14}: template target {} but P14 sample target {}", x["target"], y["target"]));
                }
                if y["value"]["re"] != x["value"]["re"] {
                    n14_differ += 1;
                }
                // the P14 values always come from the P14 sample (the template only fixes the schema and the row order)
                row["value"] = y["value"].clone();
                row["provenance"] = json!({"tool": p14_tool, "file": p14_path, "row": n14});
                if budget.is_some() {
                    let lb = *bmap.get(&(String::new(), t.clone())).ok_or_else(|| format!("input budget has no P14 entry {t:?}"))?;
                    row["provenance"]["digits"] = json!(((-lb) * 10.0).round() / 10.0);
                    row["provenance"]["digits_kind"] = json!("input accuracy (input budget)");
                }
                n14 += 1;
            }
            "PR" => {
                let (v, dg) = lane.get(&t).ok_or_else(|| format!("lane run has no master {t:?}"))?;
                row["value"] = json!({"re": v});
                let mut pv = tool_lane.clone();
                pv["digits_two_precision"] = dg.clone();
                if budget.is_some() {
                    let lb = *bmap.get(&("PR".to_string(), t.clone())).ok_or_else(|| format!("input budget has no PR entry {t:?}"))?;
                    let tp = dg.as_f64().unwrap_or(0.0);
                    let tot = -(10f64.powf(lb) + 10f64.powf(-tp)).log10();
                    pv["digits_input_budget"] = json!(((-lb) * 10.0).round() / 10.0);
                    pv["digits"] = json!((tot * 10.0).round() / 10.0);
                    pv["digits_kind"] = json!("-log10(input budget + two-precision error)");
                } else {
                    pv["digits"] = dg.clone();
                    pv["digits_kind"] = json!("two-precision only (propagation of input errors NOT included)");
                }
                row["provenance"] = pv;
                npr += 1;
            }
            f => return Err(format!("unknown endpoint family {f}")),
        }
        ep.push(row);
    }
    if n14 != p14rows.len() || npr != lane.len() {
        return Err(format!("row counts: P14 {n14} of {}, PR {npr} of {}", p14rows.len(), lane.len()));
    }
    let mut out = template.clone();
    out["endpoint"] = json!(ep);
    out["lane_tag"] = json!(format!("rustflow:{lane_out_path}"));
    let bsrc = budget.map(|b| b["source"].clone()).unwrap_or(Value::Null);
    out["provenance"] = json!({"input_budget_source": bsrc, "P14_rows": n14, "P14_rows_differing_from_template": n14_differ, "P14_tool": p14_tool, "P14_file": p14_path,
        "PR_rows": npr, "PR_tool": "RustFlow dispersive mode", "PR_out": lane_out_path, "PR_boundary": lane_out["provenance"]["boundary"]["V"]});
    Ok(out)
}
