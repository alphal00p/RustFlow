//! Stage driver of one dispersive sample: system -> singular points -> family structure (mod p) -> s0, S, path ->
//! zero family + soft boundary fit + head -> transport (any `Transport`) -> infinity families + matching + tail ->
//! assembly and export. Expensive stages write checkpoint files (exact hex floats) into `workdir` and are skipped
//! on a rerun with the same configuration fingerprint.
use super::boundary::BoundaryRow;
use super::dispersion::{self, log10abs, FamilyNum};
use super::field::{ModP, Scalar, Shadowed};
use super::frobenius::{self, FamilySpec};
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
    pub m0: i64,
    pub minf: i64,
    pub nlo0: i64,
    pub nlo_inf: i64,
    pub loops: u32,
    pub nthreads: usize,
    pub workdir: String,
}

impl Default for Config {
    fn default() -> Self {
        Config { w: 60, wf: 160, ratio: 3.0, m0: 24, minf: 50, nlo0: -6, nlo_inf: -2, loops: 5, nthreads: 1, workdir: ".".into() }
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

pub fn family_num(sys: &EBlockSystem, b: &laurent::Lmat<Float>, spec: &FamilySpec, prec: u32, nthreads: usize) -> Result<(FamilyNum, Value), String> {
    let lam = Float::with_val(prec, &spec.lam);
    // DISPERSIVE_PIVOTING = shadow (complete pivoting on the exact nonzeros) | shadow:threshold:<bits> | shadow:first (exact pivots)
    let mode = std::env::var("DISPERSIVE_PIVOTING").unwrap_or_default();
    let shadow = mode.starts_with("shadow");
    let pivoting = match mode.strip_prefix("shadow") {
        Some("") => frobenius::Pivoting::Complete,
        Some(":first") => frobenius::Pivoting::Structural,
        Some(x) if x.starts_with(":threshold:") => frobenius::Pivoting::Threshold(x[11..].parse().map_err(|_| format!("bad {mode}"))?),
        Some(_) => return Err(format!("DISPERSIVE_PIVOTING={mode}")),
        None => frobenius::Pivoting::from_env(),
    };
    let out = if shadow {
        // exact zero pattern from a 2^61 - 1 shadow (no spurious fill), exact pivots, MPFR values
        let bsh = laurent::laurent::<Shadowed>(sys, point_of(spec), spec.mm + 8, prec, nthreads);
        let bsh = laurent::scaled(&bsh, &spec.rho, prec);
        let lsh = Shadowed::from_rational(&spec.lam, prec);
        let o = frobenius::frobenius_with::<Shadowed>(&bsh, sys.dim(), &lsh, spec.n_lo, spec.mm, None, None, pivoting, prec, nthreads);
        if o.ranks != spec.ranks {
            return Err(format!("{}: shadowed ranks differ from the structure run", spec.name));
        }
        frobenius::shadowed_to_float(o)
    } else {
        let bs = laurent::scaled(b, &spec.rho, prec);
        frobenius::frobenius_with::<Float>(&bs, sys.dim(), &lam, spec.n_lo, spec.mm, Some(&spec.ranks), Some(&spec.pivots), pivoting, prec, nthreads)
    };
    let _ = &lam;
    let basis = frobenius::stable_basis(&out, spec, sys.dim(), prec)?;
    let d = &out.diag;
    let diag = json!({"name": spec.name, "dim": spec.dim, "pivoting": if shadow { format!("Shadowed {pivoting:?}") } else { format!("{pivoting:?}") },
        "min_structural_pivot_rel_log10": d.min_structural_pivot_rel_log10, "sum_structural_pivot_rel_log10": d.sum_structural_pivot_rel_log10, "structural_fallback_step": d.structural_fallback_step, "pole_rank": d.pole_rank, "steps": d.steps, "q_final": d.q_final, "top": d.top,
        "min_rank_gap_log10": d.min_rank_gap_log10, "min_pivot_log10": d.min_pivot_log10, "T_nnz_max": d.t_nnz_max, "seconds": d.seconds, "rows_kept": basis.len()});
    Ok((FamilyNum { spec: spec.clone(), basis }, diag))
}

fn save(path: &str, v: &Value) -> Result<(), String> {
    let tmp = format!("{path}.part");
    std::fs::write(&tmp, serde_json::to_vec(v).unwrap()).map_err(|e| e.to_string())?;
    std::fs::rename(&tmp, path).map_err(|e| e.to_string())
}
fn load(path: &str, fp: &str) -> Option<Value> {
    let v: Value = serde_json::from_slice(&std::fs::read(path).ok()?).ok()?;
    if v["fingerprint"] == fp { Some(v) } else { None }
}

/// the full sample. `bnd` = soft boundary rows (computed at >= wf digits by the caller from V)
pub fn run_sample(sys: &EBlockSystem, cfg: &Config, plan: &Plan, bnd: &[BoundaryRow], kernel: &dyn Transport, provenance: Value) -> Result<Value, String> {
    let t0 = Instant::now();
    let e = sys.dim();
    let pf = bits(cfg.wf);
    let pt = bits(cfg.w);
    let fp = format!("d0={} w={} wf={} ratio={} m0={} minf={} nlo0={} nlo_inf={} L={} nbnd={}", sys.d0, cfg.w, cfg.wf, cfg.ratio, cfg.m0, cfg.minf, cfg.nlo0, cfg.nlo_inf, cfg.loops, bnd.len());
    let mut out = json!({"schema": "v6a.disp.v1", "producer": "RustFlow dispersive (W3a)", "prec_bits": pt, "prec_bits_families": pf, "E": e,
                         "config": {"w": cfg.w, "wf": cfg.wf, "ratio": cfg.ratio, "threads": cfg.nthreads}, "provenance": provenance});
    // ---- zero family, soft fit, head ----
    let st0 = format!("{}/stage-zero.json", cfg.workdir);
    let (d0v, headv) = if let Some(v) = load(&st0, &fp) {
        log("stage zero: resumed from checkpoint");
        out["zero_diag"] = v["zero_diag"].clone();
        out["boundary_fit"] = v["boundary_fit"].clone();
        let rd = |k: &str| -> Vec<Float> { v[k].as_array().unwrap().iter().map(|x| fparse(x.as_str().unwrap(), pf)).collect() };
        (rd("d_s0"), rd("head"))
    } else {
        let tz = Instant::now();
        let b0 = laurent::laurent::<Float>(sys, Point::Zero, cfg.m0 + 8, pf, cfg.nthreads);
        let (fam0, diag0) = family_num(sys, &b0, &plan.zero, pf, cfg.nthreads)?;
        drop(b0);
        let fit = dispersion::soft_fit(&fam0, bnd, e, pf);
        let hterms = ((cfg.w as f64 + 10.0) / (1.0 / plan.s0.to_f64()).log10()) as i64 + 6;
        let (dv, hd) = dispersion::head(&fam0, &fit.csoft, &plan.s0, e, hterms, pf);
        let mut zd = diag0;
        zd["seconds_total"] = json!(tz.elapsed().as_secs_f64());
        out["zero_diag"] = zd;
        out["boundary_fit"] = json!({"rows": fit.rows, "dim": fam0.spec.dim, "residual_log10": fit.residual_log10, "worst_rel_log10": fit.worst_rel_log10});
        log(&format!("zero family + fit: {} {} ({:.1} s)", out["zero_diag"], out["boundary_fit"], t0.elapsed().as_secs_f64()));
        save(&st0, &json!({"fingerprint": fp, "zero_diag": out["zero_diag"], "boundary_fit": out["boundary_fit"],
            "d_s0": dv.iter().map(fhex).collect::<Vec<_>>(), "head": hd.iter().map(fhex).collect::<Vec<_>>()}))?;
        (dv, hd)
    };
    // ---- transport ----
    let stt = format!("{}/stage-transport.json", cfg.workdir);
    let (end, mid) = if let Some(v) = load(&stt, &fp) {
        log("stage transport: resumed from checkpoint");
        out["transport"] = v["transport"].clone();
        let rc = |k: &str| -> Vec<Complex> {
            v[k].as_array().unwrap().iter().map(|x| Complex::with_val(pt, (fparse(x[0].as_str().unwrap(), pt), fparse(x[1].as_str().unwrap(), pt)))).collect()
        };
        (rc("end"), rc("mid"))
    } else {
        // chunks of patches with a partial checkpoint after each chunk (resume at the next chunk)
        let cj = |v: &[Complex]| -> Vec<Value> { v.iter().map(|z| json!([fhex(z.real()), fhex(z.imag())])).collect() };
        let rcj = |v: &Value| -> Vec<Complex> {
            v.as_array().unwrap().iter().map(|x| Complex::with_val(pt, (fparse(x[0].as_str().unwrap(), pt), fparse(x[1].as_str().unwrap(), pt)))).collect()
        };
        let stp = format!("{}/stage-transport-partial.json", cfg.workdir);
        let chunk: usize = std::env::var("DISPERSIVE_TRANSPORT_CHUNK").ok().and_then(|v| v.parse().ok()).unwrap_or(20);
        let (mut k0, mut y, mut mid, mut secs, mut setup, mut mx) = match load(&stp, &fp) {
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
            secs += tr.patches.iter().map(|p| p.seconds).sum::<f64>();
            setup += tr.patches.iter().map(|p| p.setup_seconds).sum::<f64>();
            mx = tr.patches.iter().map(|p| p.dropped_im_over_re_log10).fold(mx, f64::max);
            k0 = hi;
            save(&stp, &json!({"fingerprint": fp, "next": k0, "y": cj(&y), "mid": cj(&mid), "seconds": secs, "setup_seconds": setup, "max_dropped": mx}))?;
            log(&format!("transport: {k0} of {} patches ({secs:.0} s in patches)", plan.patches.len()));
        }
        let midim = (0..e).map(|c| log10abs(mid[c].imag()) - log10abs(mid[c].real())).fold(-1e9f64, f64::max);
        out["transport"] = json!({"patches": plan.patches.len(), "seconds": secs, "setup_seconds": setup,
            "max_log10_im_over_re_after_detour": mx, "mid_integral_max_log10_im_over_re": midim});
        log(&format!("transport: {} ({:.1} s)", out["transport"], t0.elapsed().as_secs_f64()));
        save(&stt, &json!({"fingerprint": fp, "transport": out["transport"], "end": cj(&y), "mid": cj(&mid)}))?;
        (y, mid)
    };
    // ---- infinity families, matching, tail ----
    let ds: Vec<Float> = end.iter().map(|z| Float::with_val(pf, z.real())).collect();
    let binf = laurent::laurent::<Float>(sys, Point::Infinity, cfg.minf + 8, pf, cfg.nthreads);
    let mut fams = Vec::new();
    let mut fdiag = Vec::new();
    for spec in &plan.inf {
        let (f, d) = family_num(sys, &binf, spec, pf, cfg.nthreads)?;
        log(&format!("{}: {} ({:.1} s)", spec.name, d, t0.elapsed().as_secs_f64()));
        fdiag.push(d);
        fams.push(f);
    }
    drop(binf);
    out["inf_families"] = json!(fdiag);
    let tterms = ((cfg.w as f64 + 10.0) / plan.s_end.to_f64().log10()) as i64 + 6;
    let mt = dispersion::match_at_infinity(&fams, &ds, &plan.s_end, e, tterms, pf);
    out["matching"] = json!({"S": plan.s_end.to_string(), "rows": mt.rows, "cols": mt.cols, "residual_log10": mt.residual_log10,
                             "worst_rel_log10": mt.worst_rel_log10, "max_log10_coef": mt.max_log10_coef});
    log(&format!("matching: {}", out["matching"]));
    // ---- assembly ----
    let pre = dispersion::prefactor(&sys.d0, pf);
    let digits = (pt as f64 * 0.30103) as usize;
    let mut res = Vec::new();
    let mut jv = Vec::new();
    for c in 0..e {
        let midr = Float::with_val(pf, mid[c].real());
        let mut tot = Float::with_val(pf, &headv[c] + &midr);
        tot += &mt.tail[c];
        let parts = [log10abs(&headv[c]), log10abs(&midr), log10abs(&mt.tail[c]), log10abs(&tot)];
        res.push(json!({"I": tot.to_string_radix(10, Some(digits)), "I_im_mid": mid[c].imag().to_string_radix(10, Some(20)), "log10_parts_head_mid_tail_total": parts}));
        jv.push(Float::with_val(pf, &tot * &pre).to_string_radix(10, Some(digits)));
    }
    out["I"] = json!(res);
    out["J"] = json!(jv);
    out["J_note"] = json!("J_j(z = 1) = -(sin(pi d0/2)/pi) I_j (campaign normalisation e^(L gamma eps), Minkowski storage, as post.py)");
    out["masters"] = json!(sys.masters);
    out["S_end"] = json!(plan.s_end.to_string());
    out["s0"] = json!(plan.s0.to_string());
    out["seconds"] = json!(t0.elapsed().as_secs_f64());
    Ok(out)
}
