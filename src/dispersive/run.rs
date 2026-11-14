//! Stage driver of one dispersive sample: system -> singular points -> family structure (mod p) -> s0, S, path ->
//! zero family + soft boundary fit + head -> transport (any `Transport`) -> infinity families + matching + tail ->
//! assembly and export. Expensive stages write checkpoint files (exact hex floats) into `workdir`.
//!
//! CHECKPOINT RULE (REVIEW-B BLOCKER 2): every stage file carries one content fingerprint over the DE terms, the plan (path, family
//! specs with ranks and pivots), the boundary rows (hex values), the configuration (precisions, orders, pivoting), the transport
//! kernel id and the code version. A stage file with another fingerprint is NEVER reused: the run stops and asks for a fresh workdir.
use super::boundary::BoundaryRow;
use super::dispersion::{self, log10abs, FamilyNum};
use super::field::{ModP, Scalar, Shadowed};
use super::frobenius::{self, FamilySpec, PivotMode};
use super::laurent::{self, Point};
use super::path::{self, Patch};
use super::singular;
use super::system::EBlockSystem;
use super::transport::Transport;
use rug::{Complex, Float, Rational};
use serde_json::{json, Value};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Config {
    /// transport working digits (W of the oracle: s0, S and Taylor orders follow from it)
    pub w: u32,
    /// digits for families, boundary fit, matching, head and tail (W + guard)
    pub wf: u32,
    pub ratio: f64,
    /// Frobenius orders at 0 and at infinity (class choices: required, no class defaults)
    pub m0: i64,
    pub minf: i64,
    pub nlo0: i64,
    pub nlo_inf: i64,
    /// number of infinity families k = 1..loops (the loop number of the family)
    pub loops: u32,
    pub nthreads: usize,
    pub workdir: String,
    /// frobenius::PivotMode::parse syntax; production default frobenius::DEFAULT_PIVOTING
    pub pivoting: String,
    /// code version (git describe --always --dirty of the build tree), part of the fingerprint
    pub code_version: String,
    /// basis columns per family checked against the recursion
    pub residual_cols: usize,
    /// the run stops when the soft fit or the matching is worse than this (log10, relative): catches garbage families
    pub fail_residual_above_log10: f64,
}

impl Config {
    /// the class-specific choices are explicit; everything else has the oracle's defaults
    pub fn new(w: u32, wf: u32, m0: i64, minf: i64, loops: u32) -> Config {
        Config {
            w,
            wf,
            ratio: 3.0,
            m0,
            minf,
            nlo0: -6,
            nlo_inf: -2,
            loops,
            nthreads: 1,
            workdir: ".".into(),
            pivoting: frobenius::DEFAULT_PIVOTING.into(),
            code_version: "unknown".into(),
            residual_cols: 16,
            fail_residual_above_log10: -10.0,
        }
    }
}

pub fn bits(digits: u32) -> u32 {
    (digits as f64 * 3.3219281).ceil() as u32 + 32
}

pub fn fhex(x: &Float) -> String {
    x.to_string_radix(16, None)
}
pub fn fparse(s: &str, prec: u32) -> Float {
    Float::with_val(prec, Float::parse_radix(s, 16).unwrap())
}

pub struct Plan {
    pub sing: singular::SingularSet,
    pub zero: FamilySpec,
    pub inf: Vec<FamilySpec>,
    pub s0: Rational,
    pub s_end: Rational,
    pub patches: Vec<Patch>,
}

fn log(msg: &str) {
    eprintln!("[dispersive {:?}] {msg}", std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs());
}

/// structure (exact) and path
pub fn plan(sys: &EBlockSystem, cfg: &Config) -> Result<Plan, String> {
    let t = Instant::now();
    let sing = singular::singular_points(&sys.distinct_denominators(), 128);
    log(&format!("singular points {} (base {}), R0 {:e} Rmax {} ({:.1} s)", sing.points.len(), sing.base.len(), sing.r0(), sing.rmax(), t.elapsed().as_secs_f64()));
    let (rho0, rhoinf) = path::radii(sing.r0(), sing.rmax());
    let e = sys.dim();
    let d0 = sys.d0.clone();
    let lam0 = Rational::from(&d0 / 2);
    let b0 = laurent::laurent::<ModP>(sys, Point::Zero, cfg.m0 + 8, 0, cfg.nthreads);
    let out0 = frobenius::frobenius::<ModP>(&b0, e, &ModP::from_rational(&lam0, 0), cfg.nlo0, cfg.m0, None, 0, cfg.nthreads);
    let zero = frobenius::family_spec(&out0, e, cfg.nlo0, cfg.m0, &lam0, &rho0, "soft D/2")?;
    drop(out0);
    log(&format!("zero family structure: dim {} k_stable {} first {} ({:.1} s)", zero.dim, zero.k_stable, zero.first_order, t.elapsed().as_secs_f64()));
    let binf = laurent::laurent::<ModP>(sys, Point::Infinity, cfg.minf + 8, 0, cfg.nthreads);
    let mut inf = Vec::new();
    for k in 1..=cfg.loops as i64 {
        let lam = Rational::from(-&d0) * k / 2;
        let out = frobenius::frobenius::<ModP>(&binf, e, &ModP::from_rational(&lam, 0), cfg.nlo_inf, cfg.minf, None, 0, cfg.nthreads);
        let spec = frobenius::family_spec(&out, e, cfg.nlo_inf, cfg.minf, &lam, &rhoinf, &format!("inf k={k}"))?;
        log(&format!("inf k={k} structure: dim {} k_stable {} first {} ({:.1} s)", spec.dim, spec.k_stable, spec.first_order, t.elapsed().as_secs_f64()));
        inf.push(spec);
    }
    let n_st0 = zero.k_stable - zero.first_order + 1;
    let n_inf = inf.iter().map(|f| f.k_stable - f.first_order + 1).min().unwrap();
    let (s0, s_end) = path::endpoints(&rho0, sing.rmax(), cfg.w, n_st0, n_inf);
    let patches = path::plan(&sing.points, &sing.axis(), &s0, &s_end, cfg.w, cfg.ratio);
    log(&format!("path: s0 {} S {} patches {} (complex {}, sum orders {})", s0, s_end, patches.len(), patches.iter().filter(|p| !p.is_real()).count(), patches.iter().map(|p| p.order).sum::<usize>()));
    Ok(Plan { sing, zero, inf, s0, s_end, patches })
}

fn point_of(spec: &FamilySpec) -> Point {
    if spec.name.starts_with("inf") { Point::Infinity } else { Point::Zero }
}

/// MPFR family (plain MPFR or shadowed, per `mode`) + its recursion residual on `residual_cols` columns
pub fn family_num(sys: &EBlockSystem, b: &laurent::Lmat<Float>, spec: &FamilySpec, mode: PivotMode, residual_cols: usize, prec: u32, nthreads: usize) -> Result<(FamilyNum, Value), String> {
    let lam = Float::with_val(prec, &spec.lam);
    let bs = laurent::scaled(b, &spec.rho, prec);
    let out = match mode {
        PivotMode::Shadow(p) => {
            let bsh = laurent::scaled(&laurent::laurent::<Shadowed>(sys, point_of(spec), spec.mm + 8, prec, nthreads), &spec.rho, prec);
            let o = frobenius::frobenius_with::<Shadowed>(&bsh, sys.dim(), &Shadowed::from_rational(&spec.lam, prec), spec.n_lo, spec.mm, None, None, p, prec, nthreads);
            if o.ranks != spec.ranks {
                return Err(format!("{}: shadowed ranks differ from the structure run", spec.name));
            }
            frobenius::shadowed_to_float(o)
        }
        PivotMode::Float(p) => frobenius::frobenius_with::<Float>(&bs, sys.dim(), &lam, spec.n_lo, spec.mm, Some(&spec.ranks), Some(&spec.pivots), p, prec, nthreads),
    };
    let basis = frobenius::stable_basis(&out, spec, sys.dim(), prec)?;
    let d = &out.diag;
    let fam = FamilyNum { spec: spec.clone(), basis };
    let tr = Instant::now();
    let res = dispersion::family_residual(&bs, &fam, sys.dim(), residual_cols, prec);
    let nz = |x: f64| -> Value { if x.is_finite() { json!(x) } else { json!(null) } };
    let diag = json!({"name": spec.name, "dim": spec.dim, "pivoting": mode.label(), "pole_rank": d.pole_rank, "steps": d.steps, "q_final": d.q_final,
        "top": d.top, "T_nnz_max": d.t_nnz_max, "seconds": d.seconds, "rows_kept": fam.basis.len(),
        "min_pivot_log10": nz(d.min_pivot_log10),
        "min_rank_gap_log10": if matches!(mode, PivotMode::Shadow(_)) { json!("exact rank (shadow residues)") } else { nz(d.min_rank_gap_log10) },
        "min_structural_pivot_rel_log10": nz(d.min_structural_pivot_rel_log10), "structural_fallback_step": d.structural_fallback_step,
        "de_residual_log10": res, "de_residual_cols": residual_cols.min(spec.dim), "de_residual_seconds": tr.elapsed().as_secs_f64()});
    Ok((fam, diag))
}

// ---------------- fingerprint (dependency-free 128-bit FNV-1a, two lanes) ----------------
struct Fp(u64, u64);
impl Fp {
    fn new() -> Fp {
        Fp(0xcbf29ce484222325, 0x84222325cbf29ce4)
    }
    fn bytes(&mut self, b: &[u8]) {
        for &x in b {
            self.0 = (self.0 ^ x as u64).wrapping_mul(0x100000001b3);
            self.1 = (self.1 ^ (x as u64).rotate_left(5)).wrapping_mul(0x100000001b3 ^ 0x9e3779b97f4a7c15);
        }
        self.0 = (self.0 ^ 0xff).wrapping_mul(0x100000001b3);
    }
    fn s(&mut self, x: &str) {
        self.bytes(x.as_bytes());
    }
    fn hex(&self) -> String {
        format!("{:016x}{:016x}", self.0, self.1)
    }
}

/// content fingerprint of a sample: DE, plan, boundary, configuration, kernel, code
pub fn fingerprint(sys: &EBlockSystem, cfg: &Config, plan: &Plan, bnd: &[BoundaryRow], kernel: &dyn Transport) -> String {
    let mut h = Fp::new();
    h.s("w3a-dispersive-fingerprint-v2");
    h.s(&sys.d0.to_string());
    h.s(&sys.e_slot.to_string());
    for m in &sys.masters {
        h.s(&format!("{m:?}"));
    }
    for t in &sys.terms {
        h.s(&format!("{} {}", t.i, t.j));
        for p in [&t.num, &t.den] {
            h.s(&p.content.to_string());
            for c in &p.prim {
                h.s(&c.to_string());
            }
        }
    }
    h.s(&format!("s0 {} S {}", plan.s0, plan.s_end));
    for p in &plan.patches {
        h.s(&format!("{} {} {} {} {}", p.centre.0, p.centre.1, p.step.0, p.step.1, p.order));
    }
    for f in std::iter::once(&plan.zero).chain(&plan.inf) {
        h.s(&serde_json::to_string(&f.to_json()).unwrap());
        h.s(&format!("{:?}", f.pivots));
    }
    for b in bnd {
        h.s(&format!("{} {} {}", b.comp, b.n, fhex(&b.value)));
    }
    h.s(&format!("w {} wf {} ratio {} m0 {} minf {} nlo0 {} nlo_inf {} L {} pivoting {} kernel {} code {}", cfg.w, cfg.wf, cfg.ratio, cfg.m0, cfg.minf,
        cfg.nlo0, cfg.nlo_inf, cfg.loops, cfg.pivoting, kernel.id(), cfg.code_version));
    h.hex()
}

fn save(path: &str, v: &Value) -> Result<(), String> {
    let tmp = format!("{path}.part");
    std::fs::write(&tmp, serde_json::to_vec(v).unwrap()).map_err(|e| e.to_string())?;
    std::fs::rename(&tmp, path).map_err(|e| e.to_string())
}

/// Ok(Some) = resume, Ok(None) = no stage file, Err = a stage file of ANOTHER configuration/input (never reused)
fn load(path: &str, fp: &str) -> Result<Option<Value>, String> {
    let Ok(bytes) = std::fs::read(path) else { return Ok(None) };
    let v: Value = serde_json::from_slice(&bytes).map_err(|e| format!("{path}: unreadable stage file ({e}); remove it"))?;
    if v["fingerprint"] == fp {
        Ok(Some(v))
    } else {
        Err(format!("{path}: stage file of another input/configuration (fingerprint {} != {fp}); refusing to resume: use a fresh workdir", v["fingerprint"]))
    }
}

/// the full sample. `bnd` = soft boundary rows (computed at >= wf digits by the caller); `provenance` must name the V source
pub fn run_sample(sys: &EBlockSystem, cfg: &Config, plan: &Plan, bnd: &[BoundaryRow], kernel: &dyn Transport, provenance: Value) -> Result<Value, String> {
    let t0 = Instant::now();
    let e = sys.dim();
    let pf = bits(cfg.wf);
    let pt = bits(cfg.w);
    let mode = PivotMode::parse(&cfg.pivoting)?;
    let fp = fingerprint(sys, cfg, plan, bnd, kernel);
    let boundary_source = provenance.get("boundary").cloned().unwrap_or(json!("unspecified"));
    log(&format!("fingerprint {fp} (pivoting {}, kernel {}, code {})", mode.label(), kernel.id(), cfg.code_version));
    let mut out = json!({"schema": "v6a.disp.v1", "producer": "RustFlow dispersive (W3a)", "prec_bits": pt, "prec_bits_families": pf, "E": e,
                         "config": {"w": cfg.w, "wf": cfg.wf, "ratio": cfg.ratio, "m0": cfg.m0, "minf": cfg.minf, "loops": cfg.loops, "e_slot": sys.e_slot,
                                    "pivoting": mode.label(), "threads": cfg.nthreads, "kernel": kernel.id(), "code_version": cfg.code_version},
                         "fingerprint": fp, "provenance": provenance});
    // ---- zero family, soft fit, head ----
    let st0 = format!("{}/stage-zero.json", cfg.workdir);
    let (d0v, headv) = if let Some(v) = load(&st0, &fp)? {
        log("stage zero: resumed from checkpoint");
        out["zero_diag"] = v["zero_diag"].clone();
        out["boundary_fit"] = v["boundary_fit"].clone();
        let rd = |k: &str| -> Vec<Float> { v[k].as_array().unwrap().iter().map(|x| fparse(x.as_str().unwrap(), pf)).collect() };
        (rd("d_s0"), rd("head"))
    } else {
        let tz = Instant::now();
        let b0 = laurent::laurent::<Float>(sys, Point::Zero, cfg.m0 + 8, pf, cfg.nthreads);
        let (fam0, diag0) = family_num(sys, &b0, &plan.zero, mode, cfg.residual_cols, pf, cfg.nthreads)?;
        drop(b0);
        let fit = dispersion::soft_fit(&fam0, bnd, e, pf);
        if fit.rows < fam0.spec.dim {
            return Err(format!("soft fit: {} boundary rows inside the kept orders for a {}-dimensional family", fit.rows, fam0.spec.dim));
        }
        let hterms = ((cfg.w as f64 + 10.0) / (1.0 / plan.s0.to_f64()).log10()) as i64 + 6;
        let (dv, hd) = dispersion::head(&fam0, &fit.csoft, &plan.s0, e, hterms, pf);
        let mut zd = diag0;
        zd["seconds_total"] = json!(tz.elapsed().as_secs_f64());
        out["zero_diag"] = zd;
        out["boundary_fit"] = json!({"rows": fit.rows, "dim": fam0.spec.dim, "residual_log10": fit.residual_log10, "worst_rel_log10": fit.worst_rel_log10});
        log(&format!("zero family + fit: {} {} ({:.1} s)", out["zero_diag"], out["boundary_fit"], t0.elapsed().as_secs_f64()));
        save(&st0, &json!({"fingerprint": fp, "boundary_source": boundary_source, "zero_diag": out["zero_diag"], "boundary_fit": out["boundary_fit"],
            "d_s0": dv.iter().map(fhex).collect::<Vec<_>>(), "head": hd.iter().map(fhex).collect::<Vec<_>>()}))?;
        (dv, hd)
    };
    let fitw = out["boundary_fit"]["worst_rel_log10"].as_f64().unwrap_or(0.0);
    if fitw > cfg.fail_residual_above_log10 {
        return Err(format!("soft fit worst relative residual 1e{fitw:.1} above the bar 1e{}: family or boundary unusable", cfg.fail_residual_above_log10));
    }
    // ---- transport, in chunks with a partial checkpoint after each chunk ----
    let stt = format!("{}/stage-transport.json", cfg.workdir);
    let cj = |v: &[Complex]| -> Vec<Value> { v.iter().map(|z| json!([fhex(z.real()), fhex(z.imag())])).collect() };
    let rcj = |v: &Value| -> Vec<Complex> {
        v.as_array().unwrap().iter().map(|x| Complex::with_val(pt, (fparse(x[0].as_str().unwrap(), pt), fparse(x[1].as_str().unwrap(), pt)))).collect()
    };
    let (end, mid) = if let Some(v) = load(&stt, &fp)? {
        log("stage transport: resumed from checkpoint");
        out["transport"] = v["transport"].clone();
        (rcj(&v["end"]), rcj(&v["mid"]))
    } else {
        let stp = format!("{}/stage-transport-partial.json", cfg.workdir);
        let chunk: usize = std::env::var("DISPERSIVE_TRANSPORT_CHUNK").ok().and_then(|v| v.parse().ok()).unwrap_or(20);
        let (mut k0, mut y, mut mid, mut secs, mut setup, mut mx) = match load(&stp, &fp)? {
            Some(v) => {
                log(&format!("stage transport: resuming at patch {}", v["next"]));
                (v["next"].as_u64().unwrap() as usize, rcj(&v["y"]), rcj(&v["mid"]), v["seconds"].as_f64().unwrap(), v["setup_seconds"].as_f64().unwrap(), v["max_dropped"].as_f64().unwrap())
            }
            None => (0, d0v.iter().map(|x| Complex::with_val(pt, (x, 0))).collect(), vec![Complex::new(pt); e], 0.0, 0.0, -1e9f64),
        };
        while k0 < plan.patches.len() {
            let hi = (k0 + chunk).min(plan.patches.len());
            let tr = kernel.transport(sys, &y, &plan.patches[k0..hi], pt)?;
            y = tr.end;
            for (m, x) in mid.iter_mut().zip(&tr.weighted_integral) {
                *m += x;
            }
            if y.iter().chain(&mid).any(|z| !z.real().is_finite() || !z.imag().is_finite()) {
                return Err(format!("transport produced a non-finite value in patches {k0}..{hi}"));
            }
            secs += tr.patches.iter().map(|p| p.seconds).sum::<f64>();
            setup += tr.patches.iter().map(|p| p.setup_seconds).sum::<f64>();
            mx = tr.patches.iter().map(|p| p.dropped_im_over_re_log10).fold(mx, f64::max);
            k0 = hi;
            save(&stp, &json!({"fingerprint": fp, "next": k0, "y": cj(&y), "mid": cj(&mid), "seconds": secs, "setup_seconds": setup, "max_dropped": mx}))?;
            log(&format!("transport: {k0} of {} patches ({secs:.0} s in patches)", plan.patches.len()));
        }
        let midim = (0..e).map(|c| log10abs(mid[c].imag()) - log10abs(mid[c].real())).fold(-1e9f64, f64::max);
        out["transport"] = json!({"patches": plan.patches.len(), "seconds": secs, "setup_seconds": setup,
            "max_log10_im_over_re_after_detour": mx, "mid_integral_max_log10_im_over_re": midim,
            // the real part is taken after a detour: valid only if the discontinuity is analytic through the axis points (gross check)
            "detour_check": if mx < -10.0 && midim < -10.0 { "pass" } else { "FAIL: discontinuity not real after the detours" }});
        log(&format!("transport: {} ({:.1} s)", out["transport"], t0.elapsed().as_secs_f64()));
        save(&stt, &json!({"fingerprint": fp, "transport": out["transport"], "end": cj(&y), "mid": cj(&mid)}))?;
        (y, mid)
    };
    // ---- infinity families (one checkpoint per family: matching columns and tails), matching, tail ----
    let ds: Vec<Float> = end.iter().map(|z| Float::with_val(pf, z.real())).collect();
    let tterms = ((cfg.w as f64 + 10.0) / plan.s_end.to_f64().log10()) as i64 + 6;
    let mut binf: Option<laurent::Lmat<Float>> = None;
    let mut cols: Vec<Vec<Float>> = Vec::new();
    let mut tails: Vec<Vec<Float>> = Vec::new();
    let mut fdiag = Vec::new();
    for (k, spec) in plan.inf.iter().enumerate() {
        let sf = format!("{}/stage-inf-{}.json", cfg.workdir, k + 1);
        let rv = |v: &Value| -> Vec<Vec<Float>> { v.as_array().unwrap().iter().map(|c| c.as_array().unwrap().iter().map(|x| fparse(x.as_str().unwrap(), pf)).collect()).collect() };
        if let Some(v) = load(&sf, &fp)? {
            log(&format!("{}: resumed from checkpoint", spec.name));
            cols.extend(rv(&v["cols"]));
            tails.extend(rv(&v["tails"]));
            fdiag.push(v["diag"].clone());
            continue;
        }
        if binf.is_none() {
            binf = Some(laurent::laurent::<Float>(sys, Point::Infinity, cfg.minf + 8, pf, cfg.nthreads));
        }
        let (f, d) = family_num(sys, binf.as_ref().unwrap(), spec, mode, cfg.residual_cols, pf, cfg.nthreads)?;
        let (c, t) = dispersion::family_columns(&f, &plan.s_end, e, tterms, pf);
        log(&format!("{}: {} ({:.1} s)", spec.name, d, t0.elapsed().as_secs_f64()));
        let wv = |m: &[Vec<Float>]| -> Value { json!(m.iter().map(|c| c.iter().map(fhex).collect::<Vec<_>>()).collect::<Vec<_>>()) };
        save(&sf, &json!({"fingerprint": fp, "diag": d, "cols": wv(&c), "tails": wv(&t)}))?;
        cols.extend(c);
        tails.extend(t);
        fdiag.push(d);
    }
    drop(binf);
    out["inf_families"] = json!(fdiag);
    let mt = dispersion::match_columns(&cols, &tails, &ds, e, pf);
    out["matching"] = json!({"S": plan.s_end.to_string(), "rows": mt.rows, "cols": mt.cols, "residual_log10": mt.residual_log10,
                             "worst_rel_log10": mt.worst_rel_log10, "max_log10_coef": mt.max_log10_coef});
    log(&format!("matching: {}", out["matching"]));
    if mt.worst_rel_log10 > cfg.fail_residual_above_log10 {
        return Err(format!("matching worst relative residual 1e{:.1} above the bar", mt.worst_rel_log10));
    }
    // ---- assembly ----
    let pre = dispersion::prefactor(&sys.d0, pf);
    let digits = (pt as f64 * 0.30103) as usize;
    // residual-based ESTIMATE only: the fit/matching residuals overstate the per-entry accuracy by ~25-40 digits (REVIEW-B MAJOR 3);
    // the real per-entry digits come from a two-precision comparison (driver option) or from an oracle
    let est = (-(fitw.max(mt.worst_rel_log10)) - 40.0).max(0.0).min(cfg.w as f64);
    let mut res = Vec::new();
    let mut jv = Vec::new();
    for c in 0..e {
        let midr = Float::with_val(pf, mid[c].real());
        let mut tot = Float::with_val(pf, &headv[c] + &midr);
        tot += &mt.tail[c];
        if !tot.is_finite() {
            return Err(format!("entry {c}: non-finite I"));
        }
        let parts = [log10abs(&headv[c]), log10abs(&midr), log10abs(&mt.tail[c]), log10abs(&tot)];
        let canc = parts[..3].iter().cloned().fold(-1e9f64, f64::max) - parts[3];
        res.push(json!({"I": tot.to_string_radix(10, Some(digits)), "I_im_mid": mid[c].imag().to_string_radix(10, Some(20)),
            "log10_parts_head_mid_tail_total": parts, "digits_estimate": (est - canc.max(0.0)).max(0.0)}));
        jv.push(Float::with_val(pf, &tot * &pre).to_string_radix(10, Some(digits)));
    }
    out["I"] = json!(res);
    out["digits_note"] = json!("I strings carry the full transport precision; 'digits_estimate' = min(fit, matching) residual digits - 40 - head/mid/tail cancellation (an ESTIMATE that may still overstate); 'digits' (when present) comes from a two-precision comparison");
    out["J"] = json!(jv);
    out["J_note"] = json!("J_j(z = 1) = -(sin(pi d0/2)/pi) I_j (campaign normalisation e^(L gamma eps), Minkowski storage, as post.py)");
    out["masters"] = json!(sys.masters);
    out["d0"] = json!(sys.d0.to_string());
    out["eps"] = json!(sys.epsilon().to_string());
    out["S_end"] = json!(plan.s_end.to_string());
    out["s0"] = json!(plan.s0.to_string());
    out["seconds"] = json!(t0.elapsed().as_secs_f64());
    // post.py reads masters and d0 from inp/<tag>.json next to the out file
    save(&format!("{}/inp.json", cfg.workdir), &json!({"d0": sys.d0.to_string(), "eps": sys.epsilon().to_string(), "E": e, "masters": sys.masters}))?;
    Ok(out)
}

/// per-entry digits of two runs of the same sample at different precisions: -log10 |I_a - I_b| / |I_b|
pub fn two_precision_digits(a: &Value, b: &Value) -> Vec<f64> {
    let n = a["I"].as_array().map(|v| v.len()).unwrap_or(0);
    (0..n)
        .map(|i| {
            let x = Float::with_val(4000, Float::parse(a["I"][i]["I"].as_str().unwrap()).unwrap());
            let y = Float::with_val(4000, Float::parse(b["I"][i]["I"].as_str().unwrap()).unwrap());
            let d = Float::with_val(4000, (x - &y) / y).abs();
            if d.is_zero() { 99.0 } else { -d.log10().to_f64() }
        })
        .collect()
}
