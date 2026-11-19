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
use rug::ops::Pow;
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

/// one summary line of a finished run: accepted, digits_estimate and the check that limits it
pub fn summary_line(v: &Value) -> String {
    let est = v["digits_estimate"].as_f64().map_or("-".to_string(), |x| format!("{x:.1}"));
    let weakest = v["digits_by_check"].as_object().and_then(|m| m.iter().filter_map(|(k, x)| x.as_f64().map(|x| (k.clone(), x))).min_by(|a, b| a.1.total_cmp(&b.1)));
    format!("accepted {}, digits_estimate {est}{}", v["accepted"], weakest.map_or(String::new(), |(k, x)| format!(" (limited by {k}: {x:.1})")))
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
/// sha256 (hex) of a file's bytes
pub fn sha256_file(p: &str) -> Result<String, String> {
    use sha2::Digest;
    let b = std::fs::read(p).map_err(|e| format!("{p}: {e}"))?;
    Ok(format!("{:x}", sha2::Sha256::digest(&b)))
}

/// canonical serialisation (tools/pflow/pflowio.py `canonical` writes the same bytes): object keys sorted at every level, no
/// spaces, strings JSON-escaped (UTF-8 kept), every NUMBER as its exact rational ("p/q" or "p"), so that the float formatting of
/// the writer does not matter. The input of `sha256_json`.
pub fn canonical_json(v: &Value) -> String {
    match v {
        Value::Object(m) => {
            let mut ks: Vec<&String> = m.keys().collect();
            ks.sort();
            format!("{{{}}}", ks.iter().map(|k| format!("{}:{}", Value::String((*k).clone()), canonical_json(&m[*k]))).collect::<Vec<_>>().join(","))
        }
        Value::Array(a) => format!("[{}]", a.iter().map(canonical_json).collect::<Vec<_>>().join(",")),
        Value::Number(n) => {
            if let Some(i) = n.as_i64() {
                i.to_string()
            } else if let Some(u) = n.as_u64() {
                u.to_string()
            } else {
                Rational::from_f64(n.as_f64().unwrap_or(f64::NAN)).map_or("null".to_string(), |r| r.to_string())
            }
        }
        other => other.to_string(),
    }
}

/// sha256 of the canonical serialisation of a document without its `content_sha256` field
pub fn content_sha256(doc: &Value) -> String {
    let mut d = doc.clone();
    if let Some(m) = d.as_object_mut() {
        m.remove("content_sha256");
    }
    sha256_json(&d)
}

/// an integer that does not fit 64 bits: serde_json reads it as a float, so it is not exactly representable here
fn has_wide_integer(v: &Value) -> bool {
    match v {
        Value::Object(m) => m.values().any(has_wide_integer),
        Value::Array(a) => a.iter().any(has_wide_integer),
        Value::Number(n) => n.as_i64().is_none() && n.as_u64().is_none() && n.as_f64().is_some_and(|x| x.is_finite() && x.fract() == 0.0 && x.abs() >= 9.2e18),
        _ => false,
    }
}

/// refuse a document whose `content_sha256` is missing or does not match its content
pub fn check_seal(what: &str, doc: &Value) -> Result<(), String> {
    let got = doc["content_sha256"].as_str().ok_or(format!("{what}: no content_sha256 (unsealed file)"))?;
    if got != content_sha256(doc) {
        if has_wide_integer(doc) {
            return Err(format!("{what}: contains an integer beyond 64 bits, which RustFlow cannot represent exactly (write such numbers as strings); refused"));
        }
        return Err(format!("{what}: content_sha256 does not match the content (the file was edited)"));
    }
    Ok(())
}

/// every {"file": F, "sha256": H} record inside v
pub fn file_refs(v: &Value, out: &mut Vec<(String, String)>) {
    match v {
        Value::Object(m) => {
            if let (Some(f), Some(h)) = (m.get("file").and_then(|x| x.as_str()), m.get("sha256").and_then(|x| x.as_str())) {
                out.push((f.to_string(), h.to_string()));
            }
            for x in m.values() {
                file_refs(x, out);
            }
        }
        Value::Array(a) => a.iter().for_each(|x| file_refs(x, out)),
        _ => {}
    }
}

/// recompute the sha256 of every referenced source that is reachable (a mismatch is an error: stale or edited source); an
/// unreachable source is returned as an explicit downgrade
pub fn verify_sources(what: &str, v: &Value) -> Result<Vec<String>, String> {
    let mut refs = Vec::new();
    file_refs(v, &mut refs);
    let mut downs = Vec::new();
    for (f, h) in refs {
        if Path::new(&f).is_file() {
            let got = sha256_file(&f)?;
            if got != h {
                return Err(format!("{what}: source {f} has sha256 {}..., the file records {}... (stale or edited)", &got[..16], &h[..h.len().min(16)]));
            }
        } else {
            downs.push(format!("source not reachable, digest unverified: {f}"));
        }
    }
    Ok(downs)
}

/// sha256 (hex) of the canonical JSON of v
pub fn sha256_json(v: &Value) -> String {
    use sha2::Digest;
    format!("{:x}", sha2::Sha256::digest(canonical_json(v).as_bytes()))
}

/// the frame of a sample: everything that defines WHAT is computed, independent of d0 and of the numerics (normalisation, the
/// graph, the massive line, the master list). Samples of one grid must share its digest.
/// `boundary` = the d0-independent identity of the provided table (its master definitions and answer meaning digests and its
/// normalisation), or null without a table: samples built from different master definitions never share a frame digest.
pub fn frame_of(g: &Graph, e: usize, masters: &[Vec<i64>], boundary: Value) -> Value {
    json!({"normalisation": "E, P^2=1", "graph": {"edges": g.edges.iter().map(|x| [x.0, x.1]).collect::<Vec<_>>(), "ext": [g.ext.0, g.ext.1], "loops": g.loops}, "e_slot": e, "masters": masters, "boundary": boundary})
}

/// the DE files of a --de specification (one file, a comma list, or a manifest {"blocks": [...]})
pub fn block_paths(spec: &str) -> Result<Vec<String>, String> {
    Ok(if spec.ends_with(".json") && !spec.contains(',') {
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
    })
}

pub fn load_blocks(spec: &str) -> Result<(RuleSet, String), String> {
    let paths = block_paths(spec)?;
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

    let de_spec = o.path("--de").ok_or("--de")?;
    let (rs, var) = load_blocks(&de_spec)?;
    // provenance of every input (sha256 of the bytes read)
    let mut inputs = json!({"de": block_paths(&de_spec)?.iter().map(|p| Ok(json!({"file": p, "sha256": sha256_file(p)?}))).collect::<Result<Vec<_>, String>>()?});
    if de_spec.ends_with(".json") && !de_spec.contains(',') && block_paths(&de_spec)?.first() != Some(&de_spec) {
        inputs["de_manifest"] = json!({"file": de_spec, "sha256": sha256_file(&de_spec)?});
    }
    if let Some(i) = argv.iter().position(|x| x == "--class-config") {
        inputs["class_config"] = json!({"file": argv[i + 1], "sha256": sha256_file(&argv[i + 1])?});
    }
    let pnum = |k: &str, d: &str| -> Result<i64, String> { arg(k).unwrap_or(d.into()).parse::<i64>().map_err(|_| format!("{k}: not an integer")) };
    let prat = |k: &str, d: &str| -> Result<Rational, String> { arg(k).unwrap_or(d.into()).parse::<Rational>().map_err(|_| format!("{k}: not a rational")) };
    let e = pnum("--e", "-1")?;
    if e < 0 {
        return Err("--e (massive line slot) is required".into());
    }
    let e = e as usize;
    let sys = EBlockSystem::from_rules_with(&rs, e, &|c: &str| parse_coefficient(c, &var))?;
    let g = match o.cfg.get("graph") {
        Some(v) if v.is_object() && !argv.iter().any(|x| x == "--graph") => {
            inputs["graph"] = json!({"inline_sha256": sha256_json(v)});
            graph_from(v)?
        }
        _ => {
            let gp = o.path("--graph").ok_or("--graph")?;
            inputs["graph"] = json!({"file": gp, "sha256": sha256_file(&gp)?});
            graph_from(&read_json(&gp)?)?
        }
    };
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
    cfg.target_digits = match arg("--target-digits") { Some(t) => Some(t.parse::<f64>().map_err(|_| "--target-digits: a number")?), None => None };
    cfg.allow_skipped = argv.iter().any(|x| x == "--allow-skipped") || o.cfg["allow_skipped"] == true;
    cfg.allow_unverified_provided = argv.iter().any(|x| x == "--allow-unverified-provided") || o.cfg["allow_unverified_provided"] == true;
    let accept_diag = argv.iter().any(|x| x == "--accept-diagnostics") || o.cfg["accept_diagnostics"] == true;
    let workdir = o.path("--workdir");
    if let Some(f) = o.path("--provided") {
        let v = read_json(&f)?;
        cfg.provided = Provided::from_json(&v, crate::dispersive::run::bits(cfg.wf), &rs.d0).map_err(|e| format!("{f}: {e}"))?;
        inputs["provided"] = json!({"file": f, "sha256": sha256_file(&f)?, "master_definitions_sha256": v["master_definitions_sha256"], "sources": v["sources"]});
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
    // --third-pass: a third pass at w - 2K measures the convergence rate, so that the higher pass gets its own bound
    let third = argv.iter().any(|x| x == "--third-pass") || o.cfg["third_pass"] == true;
    if third && (tp == 0 || 2 * tp as u32 + 10 > cfg.w || 2 * tp > cfg.m0 - 10) {
        return Err(format!("--third-pass needs --two-precision K > 0 with w - 2K >= 10 and m - 2K >= 10"));
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
    let res3 = if third {
        let mut c3 = cfg.clone();
        c3.w -= 2 * tp as u32;
        c3.wf -= 2 * tp as u32;
        c3.m0 -= 2 * tp;
        c3.minf -= 2 * tp;
        c3.mthr -= 2 * tp;
        c3.workdir = pass_dir(&c3)?;
        Some(run_b(&sys, &g, e, &c3)?)
    } else {
        None
    };
    // the achieved absolute uncertainty of every output (of the HIGHEST pass, the one whose values are written), the larger of
    // - 10^-digits_estimate max(|v|, 10^-30 scale of the moment) (digits_estimate = min over every gate quantity), and
    // - with two passes: d1 = |v(w) - v(w-K)|. It measures the error of the LOWER pass, so it is a conservative bound for v(w);
    //   it is also written as abs_err_lower_pass;
    // - with a third pass (w - 2K): d2 = |v(w-K) - v(w-2K)|, q = d1/d2. If the measured rate log10(1/q) is at least 2 digits per
    //   step (fast geometric convergence; measured 5 per K = 8 at 2L, since w, wf and the series orders all drop together), the
    //   higher pass is bounded by the geometric tail d1 q/(1-q) with a safety factor 10; otherwise d1 is kept.
    // abs_err_kind names the rule used.
    let est = res.diag["digits_estimate"].as_f64().ok_or("digits_estimate missing")?;
    let rel = |prec: u32| rug::Float::with_val(prec, 10f64).pow(-est);
    let mut merrs: Vec<(i32, Vec<rug::Float>)> = Vec::new();
    let mut lower: Vec<Vec<Option<rug::Float>>> = Vec::new();
    let mut kinds: Vec<Vec<String>> = Vec::new();
    for (k, (m, v)) in res.moments.iter().enumerate() {
        let prec = v[0].prec();
        let sc = v.iter().map(|x| rug::Float::with_val(prec, x.abs_ref())).fold(rug::Float::new(prec), |a, b| if b > a { b } else { a });
        let floor = rug::Float::with_val(prec, &sc * rug::Float::with_val(prec, 10f64).pow(-30));
        let mut errs = Vec::new();
        let mut low = Vec::new();
        let mut kind = Vec::new();
        for (j, x) in v.iter().enumerate() {
            let a = rug::Float::with_val(prec, x.abs_ref());
            let mut e = rug::Float::with_val(prec, if a > floor { a } else { floor.clone() } * rel(prec));
            let mut kd = "run-wide digits_estimate".to_string();
            let mut lo = None;
            if let Some(r2) = &res2 {
                let d1 = rug::Float::with_val(prec, x - &r2.moments[k].1[j]).abs();
                lo = Some(d1.clone());
                let mut bound = d1.clone();
                kd = "two-precision difference (bound of the lower pass; conservative for this value)".into();
                if let Some(r3) = &res3 {
                    let d2 = rug::Float::with_val(prec, &r2.moments[k].1[j] - &r3.moments[k].1[j]).abs();
                    if !d1.is_zero() && !d2.is_zero() {
                        let rate = (rug::Float::with_val(prec, &d2 / &d1)).log10().to_f64();
                        if rate >= 2.0 {
                            let q = rug::Float::with_val(prec, &d1 / &d2);
                            let one_minus = rug::Float::with_val(prec, 1u32 - &q);
                            bound = rug::Float::with_val(prec, &d1 * &q) / one_minus * 10u32;
                            kd = format!("three-pass geometric tail (measured rate {rate:.1} digits per step of {tp}, safety 10)");
                        } else {
                            kd = format!("two-precision difference (third pass: rate {rate:.1} < 2 digits per step, not extrapolated)");
                        }
                    } else if d1.is_zero() && !d2.is_zero() {
                        bound = rug::Float::with_val(prec, &d2 * rug::Float::with_val(prec, 10f64).pow(-(tp as i32))) * 10u32;
                        kd = "three-pass: higher passes agree; d2 10^-K x 10".into();
                    }
                }
                if bound > e {
                    e = bound;
                } else {
                    kd = format!("run-wide digits_estimate (larger than the {kd})");
                }
            }
            errs.push(e);
            low.push(lo);
            kind.push(kd);
        }
        merrs.push((*m, errs));
        lower.push(low);
        kinds.push(kind);
    }
    let fe = |x: &rug::Float| x.to_string_radix(10, Some(3));
    let mut failures: Vec<Value> = res.diag["acceptance"]["failures"].as_array().cloned().unwrap_or_default();
    if let Some(r2) = &res2 {
        failures.extend(r2.diag["acceptance"]["failures"].as_array().cloned().unwrap_or_default().into_iter().map(|f| json!(format!("second pass: {}", f.as_str().unwrap_or("")))));
    }
    if let Some(r3) = &res3 {
        failures.extend(r3.diag["acceptance"]["failures"].as_array().cloned().unwrap_or_default().into_iter().map(|f| json!(format!("third pass: {}", f.as_str().unwrap_or("")))));
    }
    // what each moment is (module docs of engine.rs): m = 1 massless value; m = n + 1 >= 2: (a_e)_n / n! E[a + n 1_e];
    // m = -k <= 0: (-1)^(a_e - 1) (a_e)_n / n! E[a_e -> -n], n = k + 1 - a_e (zero if n < 0)
    let poch = |a: i64, n: i64| -> Rational { let mut r = Rational::from(1); for i in 0..n { r *= Rational::from(a + i); r /= Rational::from(i + 1); } r };
    let mut out = Vec::new();
    for (k, (m, v)) in res.moments.iter().enumerate() {
        for (j, mast) in sys.masters.iter().enumerate() {
            let sgn = if mast.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
            let prec = v[j].prec();
            let val = rug::Float::with_val(prec, &v[j] * sgn);
            let ae = &merrs[k].1[j];
            let dig = res2.as_ref().map(|r2| {
                let w = &r2.moments[k].1[j];
                let d = rug::Float::with_val(prec, &v[j] - w);
                if d.is_zero() { 999.0 } else { -(d.abs() / rug::Float::with_val(prec, v[j].abs_ref())).log10().to_f64() }
            });
            let ndig = (cfg.wf as usize).max(20);
            let mut row = json!({"m": m, "master": mast, "value": val.to_string_radix(10, Some(ndig)), "abs_err": fe(ae), "abs_err_kind": kinds[k][j], "digits_two_precision": dig});
            if let Some(lo) = &lower[k][j] {
                row["abs_err_lower_pass"] = json!(fe(lo));
            }
            let a_e = mast[e];
            if *m >= 2 {
                let n = (*m - 1) as i64;
                let mut idx = mast.clone();
                idx[e] += n;
                let f = poch(a_e, n);
                row["integral"] = json!(idx);
                row["integral_value"] = json!(rug::Float::with_val(prec, &val / &f).to_string_radix(10, Some(ndig)));
                row["integral_abs_err"] = json!(fe(&rug::Float::with_val(prec, ae / rug::Float::with_val(prec, &f).abs())));
                row["quantity"] = json!(format!("value = (a_e)_n/n! E[integral], a_e = {a_e}, n = {n}"));
            } else if *m <= 0 {
                let n = (1 - *m) as i64 - a_e;
                if n >= 0 {
                    let mut idx = mast.clone();
                    idx[e] = -n;
                    let f = Rational::from(&poch(a_e, n) * if (a_e - 1) % 2 != 0 { -1 } else { 1 });
                    row["integral"] = json!(idx);
                    row["integral_value"] = json!(rug::Float::with_val(prec, &val / &f).to_string_radix(10, Some(ndig)));
                    row["integral_abs_err"] = json!(fe(&rug::Float::with_val(prec, ae / rug::Float::with_val(prec, &f).abs())));
                    row["quantity"] = json!(format!("value = (-1)^(a_e-1) (a_e)_n/n! E[integral] (line e as a numerator of power n), a_e = {a_e}, n = {n}"));
                } else {
                    row["quantity"] = json!("no integer term at this order: the value should vanish");
                }
            }
            out.push(row);
        }
    }
    // explicit downgrades of this sample (a production contract refuses them)
    let mut downgrades: Vec<String> = cfg.provided.downgrades.iter().map(|d| format!("provided: {d}")).collect();
    if tp == 0 {
        downgrades.push("no two-precision pass: abs_err is a run-wide estimate".into());
    }
    if cfg.allow_unverified_provided {
        downgrades.push("--allow-unverified-provided".into());
    }
    if cfg.allow_skipped {
        downgrades.push("--allow-skipped".into());
    }
    if accept_diag {
        downgrades.push("--accept-diagnostics".into());
    }
    // non-master targets: --extra "a1,a2,...;b1,..." (pinched-master values from --efree FILE, schema pflow-efree.v1:
    // {"schema", "d0", "values": {"[..]": {"value", "abs_err"}}, "sources": [...]}, if needed)
    let mut extras = Vec::new();
    if let Some(ex) = arg("--extra") {
        let prec = crate::dispersive::run::bits(cfg.wf);
        let mut efree: std::collections::HashMap<Vec<i64>, rug::Float> = Default::default();
        let mut efree_downgrades: Vec<String> = Vec::new();
        let mut eferr: std::collections::HashMap<Vec<i64>, rug::Float> = Default::default();
        if let Some(f) = o.path("--efree") {
            let v = read_json(&f)?;
            if v["schema"] != "pflow-efree.v1" {
                return Err(format!("--efree {f}: schema must be \"pflow-efree.v1\" (values with abs_err, d0, sources)"));
            }
            check_seal(&format!("--efree {f}"), &v)?;
            efree_downgrades = verify_sources(&format!("--efree {f}"), &v["sources"])?;
            let fd0: Rational = v["d0"].as_str().ok_or("--efree: d0 missing")?.parse().map_err(|_| "--efree: d0")?;
            if fd0 != rs.d0 {
                return Err(format!("--efree {f}: values for d0 = {fd0}, run at d0 = {}", rs.d0));
            }
            for (k, x) in v["values"].as_object().ok_or("--efree: values {\"[indices]\": {\"value\", \"abs_err\"}}")? {
                let idx: Vec<i64> = serde_json::from_str(k).map_err(|e| format!("--efree key {k}: {e}"))?;
                let val = rug::Float::parse(x["value"].as_str().ok_or("--efree: value strings")?).map_err(|e| format!("--efree {k}: {e}"))?;
                let er = rug::Float::parse(x["abs_err"].as_str().ok_or("--efree: every value needs abs_err")?).map_err(|e| format!("--efree {k}: {e}"))?;
                efree.insert(idx.clone(), rug::Float::with_val(prec, val));
                eferr.insert(idx, rug::Float::with_val(prec, er));
            }
            inputs["efree"] = json!({"file": f, "sha256": sha256_file(&f)?, "sources": v["sources"]});
        }
        for t in ex.split(';') {
            let xv: Vec<i64> = t.split(',').map(|y| y.trim().parse::<i64>().map_err(|_| format!("--extra: {t}"))).collect::<Result<_, _>>()?;
            let r1 = super::engine::extra_value_err(&rs.rules, &sys, &xv, &res.moments, Some(&merrs), &efree, Some(&eferr), &var, prec);
            let dig = match (&r1, &res2) {
                (Ok((v, _, _)), Some(r2)) => super::engine::extra_value(&rs.rules, &sys, &xv, &r2.moments, &efree, &var, prec).ok().map(|(w, _)| {
                    let d = rug::Float::with_val(prec, v - &w);
                    if d.is_zero() { 999.0 } else { -(d.abs() / rug::Float::with_val(prec, v.abs_ref())).log10().to_f64() }
                }),
                _ => None,
            };
            downgrades.extend(efree_downgrades.drain(..).map(|d| format!("efree: {d}")));
            extras.push(match r1 {
                Ok((v, n, er)) => json!({"target": xv, "value": v.to_string_radix(10, Some((cfg.wf as usize).max(20))), "abs_err": fe(&er), "max_moment_n": n, "digits_two_precision": dig}),
                Err(err) => {
                    // a requested target that cannot be delivered: the sample is not accepted
                    failures.push(json!(format!("extra {xv:?}: {err}")));
                    json!({"target": xv, "error": err})
                }
            });
        }
    }
    let accepted = failures.is_empty();
    let boundary = if cfg.provided.values.is_empty() && cfg.provided.header.is_null() {
        Value::Null
    } else {
        json!({"master_definitions_sha256": cfg.provided.header["master_definitions_sha256"], "answer_meaning_sha256": cfg.provided.header["answer_meaning_sha256"], "normalisation": cfg.provided.header["normalisation"]})
    };
    let frame = frame_of(&g, e, &sys.masters, boundary);
    downgrades.dedup();
    let target_digits = -res.diag["gate_log10"].as_f64().unwrap_or(f64::NAN);
    let mut doc = json!({"schema": "pflow-route-b.v2", "accepted": accepted, "downgrades": downgrades, "acceptance": {"accepted": accepted, "failures": failures}, "accept_diagnostics": accept_diag,
            "d0": rs.d0.to_string(), "frame_sha256": sha256_json(&frame), "frame": frame, "inputs": inputs, "target_digits": target_digits,
            "allow_unverified_provided": cfg.allow_unverified_provided, "allow_skipped": cfg.allow_skipped,
            "digits_cap": res.diag["digits_cap"], "digits_estimate": res.diag["digits_estimate"], "digits_by_check": res.diag["digits_by_check"], "extras": extras, "e_slot": e,
            "config": {"w": cfg.w, "wf": cfg.wf, "m": cfg.m0, "s0": cfg.s0.to_string(), "S": cfg.s_end.to_string(), "r": cfg.r.to_string(), "moments": cfg.moments, "nlo": cfg.nlo, "gate_margin": cfg.gate_margin, "target_digits": cfg.target_digits, "two_precision": tp, "third_pass": third},
            "values": out, "diag": res.diag, "diag_second": res2.map(|r| r.diag), "diag_third": res3.map(|r| r.diag)});
    // the producer's content seal over everything (values, uncertainties, d0, frame, inputs, verdict)
    doc["content_sha256"] = json!(content_sha256(&doc));
    Ok(doc)
}
