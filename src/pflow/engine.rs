//! Route (B) for a one-mass massless propagator family at a fixed rational d0 (PFLOW-S1.md sections 1.4-1.7):
//!   rho(s) = -(1/pi) Im J(-s + i0),  moments  M_m = FP int_0^inf rho(s) s^-m ds  (m = 1: the massless value; m = n + 1: the
//!   z^n analytic coefficient times (-1)^n; m <= 0: the e-pinched massless integrals E[a_e -> m - 1]).
//! rho on (0, 1) = disc of the non-integer families at z = 0 (classes k h, k = 1..L), on (1, inf) = disc of the non-integer
//! families at infinity (classes -k h); both obey the homogeneous e-block. Two half-line transports (patch kernel, real
//! projection, weight (0 - z)^-m), local families at the threshold z = -1 (all classes j h mod 1, j = 0..2L), termwise end
//! pieces; head [0, s0] and tail [S, inf) termwise. The integer-family projections at z = -1 from both sides must agree
//! (continuity of Im R): reported as a check.
//! Normalisation: the DE is RustRed's (Minkowski, D_e = k^2 - z); its solution vector is I = (-1)^(sum a) E with E the Euclidean
//! value (measure d^dk/pi^(d/2), p^2 = 1). Boundary rows are converted with that sign; outputs are given in E.
use super::graph::{self, Graph};
use super::local::{self, Pt};
use crate::dispersive::boundary::BoundaryRow;
use crate::dispersive::dispersion::{self, FamilyNum, log10abs};
use crate::dispersive::field::{ModP, Scalar};
use crate::dispersive::frobenius::{self, PivotMode};
use crate::dispersive::run;
use crate::dispersive::system::EBlockSystem;
use crate::patch_kernel::transport::{Transport, TransportConfig};
use crate::patch_kernel::{PatchSpec, PlannedPath};
use rug::float::Constant;
use rug::ops::Pow;
use rug::{Complex, Float, Integer, Rational};
use serde_json::{Value, json};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct PConfig {
    /// transport working digits; families, fits and end pieces run at `wf`
    pub w: u32,
    pub wf: u32,
    pub m0: i64,
    pub minf: i64,
    pub mthr: i64,
    pub nlo: i64,
    pub loops: u32,
    /// split points: head [0, s0], threshold neighbourhoods |s - 1| < r, tail [S, inf)
    pub s0: Rational,
    pub s_end: Rational,
    pub r: Rational,
    /// patch step = distance to the nearest singular point / ratio
    pub ratio: f64,
    pub threads: usize,
    /// moments m (weights s^-m) to compute
    pub moments: Vec<i32>,
    /// all-order rows of the simple masters up to this order
    pub exact_rows: usize,
}

impl PConfig {
    pub fn new(w: u32, wf: u32, m: i64, loops: u32) -> PConfig {
        PConfig {
            w,
            wf,
            m0: m,
            minf: m,
            mthr: m,
            nlo: -10,
            loops,
            s0: Rational::from((1, 4)),
            s_end: Rational::from(4),
            r: Rational::from((1, 4)),
            ratio: 3.0,
            threads: 1,
            moments: vec![1],
            exact_rows: 40,
        }
    }
}

fn pi(prec: u32) -> Float {
    Float::with_val(prec, Constant::Pi)
}
fn frac(x: &Rational) -> Rational {
    let f = x.clone().floor();
    Rational::from(x - &f)
}

pub struct Fam {
    pub side: &'static str,
    pub class: i64,
    pub num: FamilyNum,
    pub coef: Option<Vec<Float>>, // fitted scaled coefficients (rows (n - n_lo) E + comp)
    pub diag: Value,
}

fn family(sys: &EBlockSystem, pt: &Pt, lam: &Rational, rho: &Rational, cfg: &PConfig, mm: i64, name: &str) -> Result<Option<(FamilyNum, Value)>, String> {
    let e = sys.dim();
    let off = Rational::new();
    let bp = local::lmat::<ModP>(sys, &off, pt, mm + 8, 0);
    let out = frobenius::frobenius::<ModP>(&bp, e, &ModP::from_rational(lam, 0), cfg.nlo, mm, None, 0, cfg.threads);
    if out.ranks.iter().all(|&r| r == 0) || out.q == 0 {
        return Ok(None);
    }
    let spec = match frobenius::family_spec(&out, e, cfg.nlo, mm, lam, rho, name) {
        Ok(s) => s,
        // ranks zero up to the truncation edge: no solution of this class (the threshold dimension check guards completeness)
        Err(err) if err.contains("no rank plateau") && out.t[..((((mm - cfg.nlo) / 2 + 1) as usize) * e).min(out.t.len())].iter().all(|r| r.is_empty()) => return Ok(None),
        Err(err) => return Err(format!("{name}: {err}")),
    };
    if spec.dim == 0 {
        return Ok(None);
    }
    let prec = run::bits(cfg.wf);
    let bf = local::lmat::<Float>(sys, &off, pt, mm + 8, prec);
    let mode = PivotMode::parse(frobenius::DEFAULT_PIVOTING)?;
    let (num, diag) = run::family_num(sys, &bf, &spec, mode, 8, prec, cfg.threads)?;
    Ok(Some((num, diag)))
}

/// boundary rows of one class from the graph regions (I-normalised)
fn rows(sys: &EBlockSystem, g: &Graph, e: usize, h: &Rational, at_zero: bool, k: i64, fam: &FamilyNum, cfg: &PConfig, prec: u32) -> Result<(Vec<BoundaryRow>, usize, usize), String> {
    let n_lo = fam.spec.n_lo;
    let mut out = Vec::new();
    let (mut nonzero, mut skipped) = (0usize, 0usize);
    for (j, a) in sys.masters.iter().enumerate() {
        let sgn = if a.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
        let push = |out: &mut Vec<BoundaryRow>, n: i64, v: Float| {
            if v.is_nan() {
                eprintln!("[pflow] NaN boundary row: master {a:?} class {k} at_zero {at_zero} n {n}");
                return;
            }
            out.push(BoundaryRow { comp: j, n, value: Float::with_val(prec, &v * sgn) })
        };
        if let Some(ex) = graph::exact_series(g, e, a, h, at_zero, k, cfg.exact_rows, prec) {
            if ex.is_empty() {
                for m in n_lo..n_lo + 12 {
                    push(&mut out, m, Float::new(prec));
                }
            } else {
                let nmin = ex.iter().map(|x| x.0).min().unwrap();
                for m in n_lo..nmin {
                    push(&mut out, m, Float::new(prec));
                }
                for (n, v) in ex {
                    push(&mut out, n, v);
                    nonzero += 1;
                }
            }
            continue;
        }
        match graph::regions(g, e, a, h, prec, at_zero) {
            Err(_) | Ok(None) => skipped += 1,
            Ok(Some(rs)) => {
                let mine: Vec<&(i64, i64, Float)> = rs.iter().filter(|x| x.0 == k).collect();
                if mine.is_empty() {
                    for m in n_lo..n_lo + 12 {
                        push(&mut out, m, Float::new(prec));
                    }
                } else {
                    // only the LOWEST order is clean
                    let nmin = mine.iter().map(|x| x.1).min().unwrap();
                    let mut v = Float::new(prec);
                    for x in mine.iter().filter(|x| x.1 == nmin) {
                        v += &x.2;
                    }
                    for m in n_lo..nmin {
                        push(&mut out, m, Float::new(prec));
                    }
                    push(&mut out, nmin, v);
                    nonzero += 1;
                }
            }
        }
    }
    Ok((out, nonzero, skipped))
}

/// rank-revealing least squares: Householder QR with column pivoting; columns whose remaining norm falls below
/// 2^(-prec/2) of the largest are left free (set to 0). Returns (x, rank).
pub fn lstsq_rr(a: &[Vec<Float>], b: &[Float], prec: u32) -> (Vec<Float>, usize) {
    let m = a.len();
    let n = if m > 0 { a[0].len() } else { 0 };
    let mut r: Vec<Vec<Float>> = a.to_vec();
    let mut y: Vec<Float> = b.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    let tol_exp = -(prec as i32) / 2;
    let colnorm = |r: &Vec<Vec<Float>>, k: usize, from: usize| -> Float {
        let mut s = Float::new(prec);
        for i in from..m {
            s += Float::with_val(prec, r[i][k].square_ref());
        }
        s.sqrt()
    };
    let mut big = Float::new(prec);
    for k in 0..n {
        let c = colnorm(&r, k, 0);
        if c > big {
            big = c;
        }
    }
    let mut rank = 0;
    for k in 0..n.min(m) {
        // pivot: largest remaining column
        let (mut best, mut bv) = (k, Float::new(prec));
        for j in k..n {
            let c = colnorm(&r, j, k);
            if c > bv {
                bv = c;
                best = j;
            }
        }
        if bv.is_zero() || big.is_zero() || Float::with_val(prec, &bv / &big).get_exp().unwrap_or(i32::MIN) < tol_exp {
            break;
        }
        if best != k {
            for row in r.iter_mut() {
                row.swap(k, best);
            }
            perm.swap(k, best);
        }
        // Householder on column k, rows k..m
        let mut alpha = bv.clone();
        if r[k][k] > 0 {
            alpha = -alpha;
        }
        let mut v: Vec<Float> = (k..m).map(|i| r[i][k].clone()).collect();
        v[0] -= &alpha;
        let mut vn = Float::new(prec);
        for x in &v {
            vn += Float::with_val(prec, x.square_ref());
        }
        if vn.is_zero() {
            rank += 1;
            continue;
        }
        for j in k..n {
            let mut s = Float::new(prec);
            for (t, i) in (k..m).enumerate() {
                s += Float::with_val(prec, &v[t] * &r[i][j]);
            }
            let f = Float::with_val(prec, &s * 2u32) / &vn;
            for (t, i) in (k..m).enumerate() {
                let d = Float::with_val(prec, &f * &v[t]);
                r[i][j] -= d;
            }
        }
        let mut s = Float::new(prec);
        for (t, i) in (k..m).enumerate() {
            s += Float::with_val(prec, &v[t] * &y[i]);
        }
        let f = Float::with_val(prec, &s * 2u32) / &vn;
        for (t, i) in (k..m).enumerate() {
            let d = Float::with_val(prec, &f * &v[t]);
            y[i] -= d;
        }
        rank += 1;
    }
    let mut xp = vec![Float::new(prec); n];
    for k in (0..rank).rev() {
        let mut s = y[k].clone();
        for j in k + 1..rank {
            s -= Float::with_val(prec, &r[k][j] * &xp[j]);
        }
        xp[k] = Float::with_val(prec, &s / &r[k][k]);
    }
    let mut x = vec![Float::new(prec); n];
    for k in 0..n {
        x[perm[k]] = xp[k].clone();
    }
    (x, rank)
}

/// relative-weighted fit of the boundary rows on a family basis (as dispersion::soft_fit, rank revealing)
fn fit_rows(fam: &FamilyNum, bnd: &[BoundaryRow], e: usize, prec: u32) -> (Vec<Float>, usize, usize, f64) {
    let rho = Float::with_val(prec, &fam.spec.rho);
    let mut amat = Vec::new();
    let mut bvec = Vec::new();
    for br in bnd {
        let row = ((br.n - fam.spec.n_lo) as usize) * e + br.comp;
        if br.n < fam.spec.n_lo || row >= fam.basis.len() {
            continue;
        }
        let sc = Float::with_val(prec, (&rho).pow(br.n as i32));
        let vs = Float::with_val(prec, &br.value * &sc);
        let w = if vs.is_zero() { Float::with_val(prec, 1) } else { Float::with_val(prec, 1 / Float::with_val(prec, vs.clone().abs())) };
        amat.push(fam.basis[row].iter().map(|x| Float::with_val(prec, x * &w)).collect::<Vec<_>>());
        bvec.push(Float::with_val(prec, &vs * &w));
    }
    let (x, rank) = lstsq_rr(&amat, &bvec, prec);
    let mut worst = -1e9f64;
    for (i, row) in amat.iter().enumerate() {
        let mut v = Float::with_val(prec, Float::dot(row.iter().zip(&x)));
        v -= &bvec[i];
        worst = worst.max(log10abs(&v));
    }
    let c = fam.basis.iter().map(|row| Float::with_val(prec, Float::dot(row.iter().zip(&x)))).collect();
    (c, amat.len(), rank, worst)
}

/// value of family `f` (fitted coefficients `c`) or one basis column at x (real; |x|^lam on x < 0)
fn eval_rows(f: &FamilyNum, rowsv: &dyn Fn(usize) -> Float, e: usize, x: &Float, prec: u32) -> Vec<Float> {
    let rho = Float::with_val(prec, &f.spec.rho);
    let lam = Float::with_val(prec, &f.spec.lam);
    let ax = Float::with_val(prec, x.abs_ref());
    let pre = Float::with_val(prec, (&ax).pow(&lam));
    let t = Float::with_val(prec, x / &rho);
    let nrow = f.basis.len() / e;
    let mut out = vec![Float::new(prec); e];
    let mut tp = Float::with_val(prec, (&t).pow(f.spec.n_lo as i32));
    for b in 0..nrow {
        for c in 0..e {
            out[c] += Float::with_val(prec, &rowsv(b * e + c) * &tp);
        }
        tp *= &t;
    }
    for v in out.iter_mut() {
        *v *= &pre;
    }
    out
}

/// sum over families of rho at s (0-side: families at 0 evaluated at z = -s; inf-side: at w = -1/s)
fn rho_at(fams: &[Fam], e: usize, s: &Float, prec: u32) -> Vec<Float> {
    let mut out = vec![Float::new(prec); e];
    for f in fams {
        let c = f.coef.as_ref().unwrap();
        let lam = Float::with_val(prec, &f.num.spec.lam);
        let sn = Float::with_val(prec, Float::with_val(prec, &lam * pi(prec)).sin() / pi(prec));
        let (x, sign) = if f.side == "0" {
            (Float::with_val(prec, -s), Float::with_val(prec, -&sn))
        } else {
            (Float::with_val(prec, -1 / Float::with_val(prec, s)), sn)
        };
        // family = |x|^lam sum c (x/rho)^n ;  0-side rho = -(sin/pi) s^lam f(-s) ; inf-side rho = (sin/pi) s^-mu g(-1/s)
        let v = eval_rows(&f.num, &|i| c[i].clone(), e, &x, prec);
        for i in 0..e {
            out[i] += Float::with_val(prec, &v[i] * &sign);
        }
    }
    out
}

/// generalised binomial binom(-m, j) for integer m
fn binom_neg(m: i32, j: u32) -> Rational {
    let mut r = Rational::from(1);
    for i in 0..j {
        r *= Rational::from(-(m as i64) - i as i64);
        r /= Rational::from(i as i64 + 1);
    }
    r
}

/// head [0, s0] (0-side) or tail [S, inf) (inf-side) of FP int rho s^-m ds, termwise
fn end_piece(fams: &[Fam], e: usize, at: &Rational, m: i32, prec: u32) -> Vec<Float> {
    let mut out = vec![Float::new(prec); e];
    let sf = Float::with_val(prec, at);
    for f in fams {
        let c = f.coef.as_ref().unwrap();
        let lam = Float::with_val(prec, &f.num.spec.lam);
        let rho = Float::with_val(prec, &f.num.spec.rho);
        let sn = Float::with_val(prec, Float::with_val(prec, &lam * pi(prec)).sin() / pi(prec));
        let nrow = f.basis_len() / e;
        for b in 0..nrow {
            let n = f.num.spec.n_lo + b as i64;
            // (-1/rho)^n
            let mut w = Float::with_val(prec, Float::with_val(prec, -1 / Float::with_val(prec, &rho)).pow(n as i32));
            if f.side == "0" {
                let ex = Float::with_val(prec, &lam + (n + 1 - m as i64));
                if ex.is_zero() {
                    continue;
                }
                w *= Float::with_val(prec, (&sf).pow(&ex)) / &ex;
                w *= Float::with_val(prec, -&sn);
            } else {
                let ex = Float::with_val(prec, &lam + (n + m as i64 - 1));
                if ex.is_zero() {
                    continue;
                }
                w *= Float::with_val(prec, (&sf).pow(&Float::with_val(prec, -&ex))) / &ex;
                w *= &sn;
            }
            for i in 0..e {
                out[i] += Float::with_val(prec, &c[b * e + i] * &w);
            }
        }
    }
    out
}

impl Fam {
    fn basis_len(&self) -> usize {
        self.num.basis.len()
    }
}

/// threshold neighbourhood: decompose rho at s = 1 - r (side 0) / 1 + r (side inf) in the local families at z = -1, integrate
/// FP int rho s^-m ds to s = 1 termwise. Returns (integral, decomposition residual log10, integer-family coefficients)
fn threshold_piece(tf: &[FamilyNum], y: &[Float], r: &Rational, inf_side: bool, m: i32, jmax: u32, prec: u32) -> (Vec<Float>, f64, Vec<Float>) {
    let e = y.len();
    let rf = Float::with_val(prec, r);
    let x = if inf_side { Float::with_val(prec, -&rf) } else { rf.clone() };
    let mut cols: Vec<Vec<Float>> = Vec::new();
    let mut spec: Vec<(usize, usize)> = Vec::new();
    for (fi, f) in tf.iter().enumerate() {
        for t in 0..f.spec.dim {
            cols.push(eval_rows(f, &|i| f.basis[i][t].clone(), e, &x, prec));
            spec.push((fi, t));
        }
    }
    let a: Vec<Vec<Float>> = (0..e).map(|i| cols.iter().map(|c| c[i].clone()).collect()).collect();
    let (p, _rank) = lstsq_rr(&a, y, prec);
    let mut res = Float::new(prec);
    for (i, row) in a.iter().enumerate() {
        let v = Float::with_val(prec, Float::dot(row.iter().zip(&p))) - &y[i];
        if v.clone().abs() > res { res = v.abs(); }
    }
    let ynorm = y.iter().map(log10abs).fold(-1e9f64, f64::max);
    let mut out = vec![Float::new(prec); e];
    let mut intc = Vec::new();
    for ((fi, t), pc) in spec.iter().zip(&p) {
        let f = &tf[*fi];
        if f.spec.lam == 0 {
            intc.push(pc.clone());
        }
        let lam = Float::with_val(prec, &f.spec.lam);
        let rho = Float::with_val(prec, &f.spec.rho);
        let nrow = f.basis.len() / e;
        for b in 0..nrow {
            let n = f.spec.n_lo + b as i64;
            let rowmax = (0..e).map(|i| log10abs(&f.basis[b * e + i][*t])).fold(-1e9f64, f64::max);
            if rowmax < -1e8 {
                continue;
            }
            // side 0: x^(lam) (x/rho)^n (1-x)^-m = sum_j binom(-m,j) (-1)^j x^(lam+n+j) rho^-n
            // side inf: u = -x: |x|^lam (x/rho)^n = u^lam (-u/rho)^n ; (1+u)^-m = sum_j binom(-m,j) u^j
            let base = Float::with_val(prec, Float::with_val(prec, 1 / Float::with_val(prec, &rho)).pow(n as i32));
            let sgn_n = if inf_side && n.rem_euclid(2) == 1 { -1 } else { 1 };
            let mut acc = Float::new(prec);
            for j in 0..jmax {
                let ex = Float::with_val(prec, &lam + (n + j as i64 + 1));
                if ex.is_zero() {
                    continue; // an x^-1 term: the physical weight must vanish (non-physical pole solution)
                }
                let mut term = Float::with_val(prec, binom_neg(m, j));
                if !inf_side && j % 2 == 1 {
                    term = -term;
                }
                term *= Float::with_val(prec, (&rf).pow(&ex)) / &ex;
                acc += term;
                if m >= 1 && j > 8 && log10abs(&acc) > -1e8 {
                    // geometric convergence ~ r^j
                }
            }
            let w = Float::with_val(prec, &base * &acc) * sgn_n;
            for i in 0..e {
                out[i] += Float::with_val(prec, &(Float::with_val(prec, pc * &f.basis[b * e + i][*t])) * &w);
            }
        }
    }
    (out, log10abs(&res) - ynorm, intc)
}

/// real segment a -> b (rationals), steps r/ratio with r the distance to the nearest singular point; detours (half 12-gons,
/// radius 0.4 x gap) above axis points strictly inside
fn plan_segment(a: &Rational, b: &Rational, points: &[(f64, f64)], w: u32, ratio: f64) -> Vec<PatchSpec> {
    let qf = |x: &Rational| x.to_f64();
    let rr = |x: f64| -> Rational {
        let s = (1i64 << 24) as f64;
        Rational::from(((x * s).round() as i64, 1i64 << 24))
    };
    let radius = |zr: f64, zi: f64| points.iter().map(|&(p, q)| (zr - p).hypot(zi - q)).fold(f64::INFINITY, f64::min);
    let nord = |r: f64, h: f64| -> usize { ((w as f64 + 8.0) / (r / h.abs()).log10()).ceil() as usize + 4 };
    let dir = if b > a { 1.0 } else { -1.0 };
    let axis: Vec<f64> = points.iter().filter(|p| p.1.abs() < 1e-12 && (p.0 - qf(a)) * dir > 0.0 && (qf(b) - p.0) * dir > 0.0).map(|p| p.0).collect();
    let mut axis = axis;
    axis.sort_by(|x, y| (x * dir).partial_cmp(&(y * dir)).unwrap());
    let mut out = Vec::new();
    let mut z = a.clone();
    let step_to = |z: &mut Rational, target: &Rational, out: &mut Vec<PatchSpec>| {
        while (qf(target) - qf(z)) * dir > 1e-300 {
            let zc = qf(z);
            let r = radius(zc, 0.0);
            let h = r / ratio;
            let hq = if (qf(target) - zc).abs() <= 1.5 * h { Rational::from(target - &*z) } else { rr(h * dir) };
            out.push(PatchSpec { centre: (z.clone(), Rational::new()), step: (hq.clone(), Rational::new()), order: nord(r, qf(&hq)) });
            *z += hq;
        }
    };
    for xp in axis {
        let gap = points.iter().filter(|p| (p.0 - xp).hypot(p.1) > 1e-9 * (1.0 + xp.abs())).map(|p| (p.0 - xp).hypot(p.1)).fold(f64::INFINITY, f64::min);
        let rd = 0.4 * gap;
        let start = rr(xp - dir * rd);
        let end = rr(xp + dir * rd);
        step_to(&mut z, &start, &mut out);
        let nv = 12;
        let mut cr = z.clone();
        let mut ci = Rational::new();
        for k in 1..=nv {
            let th = std::f64::consts::PI * k as f64 / nv as f64;
            let (vr, vi) = if k == nv { (end.clone(), Rational::new()) } else { (rr(xp - dir * rd * th.cos()), rr(rd * th.sin())) };
            let (hr, hi) = (Rational::from(&vr - &cr), Rational::from(&vi - &ci));
            let r = radius(qf(&cr), qf(&ci));
            let hn = qf(&hr).hypot(qf(&hi));
            out.push(PatchSpec { centre: (cr.clone(), ci.clone()), step: (hr, hi), order: nord(r, hn) });
            cr = vr;
            ci = vi;
        }
        z = cr;
    }
    step_to(&mut z, b, &mut out);
    out
}

pub struct PResult {
    pub moments: Vec<(i32, Vec<Float>)>,
    pub diag: Value,
}

/// singular points of the system (z-plane, f64) from the exact denominators (dispersive::singular works in eta = z - 1 and
/// always includes z = 0; the denominators are handed over in that variable)
pub fn singular_points(sys: &EBlockSystem) -> Vec<(f64, f64)> {
    let dens = local::denominators_z(sys, &Rational::new());
    let one = Rational::from(1);
    let ints: Vec<Vec<Integer>> = dens
        .iter()
        .map(|d| {
            let d = local::shift(d, &one); // p(z) -> p(eta + 1)
            let lcm = d.iter().fold(Integer::from(1), |acc, c| acc.lcm(c.denom()));
            d.iter().map(|c| Rational::from(c * &lcm).numer().clone()).collect::<Vec<Integer>>()
        })
        .collect();
    crate::dispersive::singular::singular_points(&ints, 256).points
}

/// route (B) at one d0
pub fn run_b(sys: &EBlockSystem, g: &Graph, e_slot: usize, cfg: &PConfig) -> Result<PResult, String> {
    let t0 = Instant::now();
    let e = sys.dim();
    let prec = run::bits(cfg.wf);
    let h = Rational::from(&sys.d0 / 2u32);
    let pts = singular_points(sys);
    let dist = |c: (f64, f64)| pts.iter().filter(|p| (p.0 - c.0).hypot(p.1 - c.1) > 1e-12).map(|p| (p.0 - c.0).hypot(p.1 - c.1)).fold(f64::INFINITY, f64::min);
    // exact dyadic rationals with 24 significant bits (never 0)
    let rq = |x: f64| -> Rational {
        let ex = x.log2().floor() as i32 - 23;
        let m = (x / 2f64.powi(ex)).floor() as i64;
        if ex >= 0 { Rational::from(m) * Rational::from(rug::Integer::from(1) << ex as u32) } else { Rational::from((m, 1)) / Rational::from(rug::Integer::from(1) << (-ex) as u32) }
    };
    let d0 = dist((0.0, 0.0));
    let dthr = dist((-1.0, 0.0));
    let rmax = pts.iter().map(|p| p.0.hypot(p.1)).fold(0.0f64, f64::max);
    let rho0 = rq(0.5 * d0);
    let rho_inf = rq(0.5 / rmax);
    let rho_thr = rq(0.5 * dthr);
    // split points adapted to the singular points closest to 0, -1 and infinity (spurious points of the Laporta basis may sit
    // at distance O(eps) from the physical ones): the families converge on the solutions' own radius, the transports detour
    let mut cfg = cfg.clone();
    let cap = |x: &Rational, lim: f64| -> Rational { if x.to_f64() > lim { rq(lim) } else { x.clone() } };
    cfg.s0 = cap(&cfg.s0, 0.4 * d0);
    cfg.r = cap(&cfg.r, 0.4 * dthr);
    if cfg.s_end.to_f64() < 2.5 * rmax {
        cfg.s_end = rq(2.5 * rmax);
    }
    let cfg = &cfg;
    let mut diag = json!({"s0": cfg.s0.to_string(), "r": cfg.r.to_string(), "S": cfg.s_end.to_string(), "singular_points": pts.iter().map(|p| vec![p.0, p.1]).collect::<Vec<_>>(), "rho0": rho0.to_string(), "rho_inf": rho_inf.to_string(), "rho_thr": rho_thr.to_string(), "families": []});
    let mut fams0 = Vec::new();
    let mut famsi = Vec::new();
    for k in 1..=cfg.loops as i64 {
        for (side, pt, lam, rho, mm) in [("0", Pt::Finite(Rational::new()), Rational::from(&h * k), &rho0, cfg.m0), ("inf", Pt::Infinity, -Rational::from(&h * k), &rho_inf, cfg.minf)] {
            let t = Instant::now();
            let Some((num, fdiag)) = family(sys, &pt, &lam, rho, cfg, mm, &format!("{side} k={k}"))? else { continue };
            let (br, nz, sk) = rows(sys, g, e_slot, &h, side == "0", k, &num, cfg, prec)?;
            let (csoft, nrows, rank, worst) = fit_rows(&num, &br, e, prec);
            diag["families"].as_array_mut().unwrap().push(json!({"side": side, "class": k, "dim": num.spec.dim, "rows": nrows, "rank": rank, "nonzero_rows": nz, "skipped_masters": sk,
                "fit_worst_rel_log10": worst, "frobenius": fdiag, "seconds": t.elapsed().as_secs_f64()}));
            eprintln!("[pflow] family {side} k={k}: dim {} rank {rank} n_lo {} k_stable {} basis rows {} fit rows {nrows} worst {worst:.1} ({:.1} s)", num.spec.dim, num.spec.n_lo, num.spec.k_stable, num.basis.len(), t.elapsed().as_secs_f64());
            let f = Fam { side: if side == "0" { "0" } else { "inf" }, class: k, num, coef: Some(csoft), diag: Value::Null };
            if side == "0" { fams0.push(f) } else { famsi.push(f) }
        }
    }
    // threshold families
    let mut tfams: Vec<FamilyNum> = Vec::new();
    let mut seen: Vec<Rational> = Vec::new();
    for j in 0..=(2 * cfg.loops as i64) {
        let lam = frac(&Rational::from(&h * j));
        if seen.contains(&lam) {
            continue;
        }
        seen.push(lam.clone());
        if let Some((num, _)) = family(sys, &Pt::Finite(Rational::from(-1)), &lam, &rho_thr, cfg, cfg.mthr, &format!("thr {j}"))? {
            tfams.push(num);
        }
    }
    let tdim: usize = tfams.iter().map(|f| f.spec.dim).sum();
    diag["threshold_dims"] = json!(tfams.iter().map(|f| (f.spec.lam.to_string(), f.spec.dim)).collect::<Vec<_>>());
    if tdim != e {
        diag["warning"] = json!(format!("threshold families span {tdim} of {e}"));
    }
    // transports
    let ksys = crate::dispersive::kernel::rational_system(sys)?;
    let tprec = run::bits(cfg.w);
    let one = Rational::from(1);
    let p0 = plan_segment(&Rational::from(-&cfg.s0), &Rational::from(&cfg.r - &one), &pts, cfg.w, cfg.ratio);
    let pi_ = plan_segment(&Rational::from(-&cfg.s_end), &Rational::from(-(Rational::from(&one + &cfg.r))), &pts, cfg.w, cfg.ratio);
    let s0f = Float::with_val(prec, &cfg.s0);
    let sef = Float::with_val(prec, &cfg.s_end);
    let y0 = rho_at(&fams0, e, &s0f, prec);
    let yi = rho_at(&famsi, e, &sef, prec);
    diag["patches"] = json!([p0.len(), pi_.len()]);
    eprintln!("[pflow] rho(s0) = {:?}", y0.iter().map(|x| x.to_f64()).collect::<Vec<_>>());
    eprintln!("[pflow] rho(S)  = {:?}", yi.iter().map(|x| x.to_f64()).collect::<Vec<_>>());
    eprintln!("[pflow] singular points {:?}; patches {} / {}", pts, p0.len(), pi_.len());
    let mut moments = Vec::new();
    let mut cont = -1e9f64;
    for &m in &cfg.moments {
        let mut tot = end_piece(&fams0, e, &cfg.s0, m, prec);
        let tl = end_piece(&famsi, e, &cfg.s_end, m, prec);
        for (half, y, patches) in [(0, &y0, &p0), (1, &yi, &pi_)] {
            let mut tc = TransportConfig::dispersive(tprec);
            tc.cauchy_pole = Some(Rational::new());
            tc.cauchy_power = m;
            tc.threads = cfg.threads;
            tc.pin_threads = false;
            let path = PlannedPath { offset: Rational::new(), patches: patches.clone() };
            let tr = Transport::new(&ksys, &path, tc).map_err(|x| x.to_string())?;
            let yc: Vec<Complex> = y.iter().map(|v| Complex::with_val(tprec, v)).collect();
            let mut st = tr.initial_state(&yc).map_err(|x| x.to_string())?;
            tr.run(&mut st, usize::MAX, &yc, &Default::default(), None, |_| {}).map_err(|x| x.to_string())?;
            let yend: Vec<Float> = st.y.iter().map(|c| Float::with_val(prec, c.real())).collect();
            let integ = st.integral.as_ref().unwrap();
            // kernel: int y (0 - z)^-m dz ; on z = -s: = -int rho s^-m ds (0-side, s0 -> 1-r) and = +int_{1+r}^{S} (inf-side)
            let mut mid = Vec::new();
            for i in 0..e {
                let v = Float::with_val(prec, integ[i].real());
                mid.push(if half == 0 { Float::with_val(prec, -&v) } else { v.clone() });
                if half == 0 { tot[i] -= v } else { tot[i] += v }
            }
            let (th, res, ic) = threshold_piece(&tfams, &yend, &cfg.r, half == 1, m, 400, prec);
            if m == cfg.moments[0] {
                let f = |v: &Vec<Float>| v.iter().map(|x| x.to_string_radix(10, Some(20))).collect::<Vec<_>>();
                diag[format!("pieces_mid{half}")] = json!(f(&mid));
                diag[format!("pieces_thr{half}")] = json!(f(&th));
                diag[format!("pieces_yend{half}")] = json!(f(&yend));
            }
            for i in 0..e {
                tot[i] += &th[i];
            }
            diag[format!("threshold_residual_{half}_m{m}")] = json!(res);
            if m == cfg.moments[0] {
                diag[format!("intcoef_{half}")] = json!(ic.iter().map(|x| x.to_string_radix(10, Some(30))).collect::<Vec<_>>());
            }
        }
        for i in 0..e {
            tot[i] += &tl[i];
        }
        if m == cfg.moments[0] {
            let f = |v: &Vec<Float>| v.iter().map(|x| x.to_string_radix(10, Some(20))).collect::<Vec<_>>();
            diag["pieces_m1"] = json!({"head": f(&end_piece(&fams0, e, &cfg.s0, m, prec)), "tail": f(&tl), "total": f(&tot)});
        }
        moments.push((m, tot));
    }
    // continuity of the integer-family projection at z = -1
    if let (Some(a), Some(b)) = (diag["intcoef_0"].as_array(), diag["intcoef_1"].as_array()) {
        for (x, y) in a.iter().zip(b) {
            let (x, y) = (Float::with_val(prec, Float::parse(x.as_str().unwrap()).unwrap()), Float::with_val(prec, Float::parse(y.as_str().unwrap()).unwrap()));
            let d = log10abs(&Float::with_val(prec, &x - &y)) - log10abs(&x).max(log10abs(&y));
            cont = cont.max(d);
        }
    }
    diag["threshold_continuity_log10"] = json!(cont);
    diag["seconds"] = json!(t0.elapsed().as_secs_f64());
    Ok(PResult { moments, diag })
}

/// massless value (E-normalisation) of a NON-master integral X of the one-mass family: X(z) = sum_j c_j(z) M_j(z) (RustRed rule,
/// Minkowski I-normalisation), analytic part of M_j: A_j(z) = sum_n a_{j,n} z^n with a_{j,n} = (-1)^n Mom_j(n + 1)
/// (moments m = n + 1 of the e-block masters, I-normalised). E_massless[X] = (-1)^(sum X) sum_j [c_j A_j]_{z^0}. Rule terms on
/// e-free (pinched) masters need their massless values (`efree`, E-normalised); missing ones are an error.
pub fn extra_value(
    rules: &std::collections::HashMap<Vec<i64>, Vec<(Vec<i64>, String)>>,
    sys: &EBlockSystem,
    x: &[i64],
    moments: &[(i32, Vec<Float>)],
    efree: &std::collections::HashMap<Vec<i64>, Float>,
    var: &str,
    prec: u32,
) -> Result<(Float, i64), String> {
    use crate::dispersive::system::parse_coefficient;
    let xv = x.to_vec();
    let terms: Vec<(Vec<i64>, String)> = if sys.masters.contains(&xv) { vec![(xv.clone(), "1".into())] } else { rules.get(&xv).cloned().ok_or("no rule for the target")? };
    let mut tot = Float::new(prec);
    let mut maxn = 0i64;
    let mom = |m: i32| moments.iter().find(|(k, _)| *k == m).map(|(_, v)| v);
    for (ti, c) in &terms {
        let (n, d) = parse_coefficient(c, var)?;
        // exact Laurent expansion of n/d at z = 0
        let mut nn = n.clone();
        let mut dd = d.clone();
        crate::dispersive::poly::trim(&mut nn);
        crate::dispersive::poly::trim(&mut dd);
        let Some(vn) = crate::dispersive::poly::valuation(&nn) else { continue };
        let vd = crate::dispersive::poly::valuation(&dd).unwrap();
        let shift = vn as i64 - vd as i64;
        let need = (-shift).max(0) as usize + 1;
        let ser = crate::dispersive::field::series_div::<Float>(
            &nn[vn..].iter().map(|q| Float::with_val(prec, q)).collect::<Vec<_>>(),
            &dd[vd..].iter().map(|q| Float::with_val(prec, q)).collect::<Vec<_>>(),
            need,
            prec,
        );
        if let Some(j) = sys.masters.iter().position(|m| m == ti) {
            // z^0 coefficient of c(z) A_j(z): sum_k g_k a_{j,-k}, k = shift + i <= 0
            for (i, g) in ser.iter().enumerate() {
                let k = shift + i as i64;
                if k > 0 {
                    break;
                }
                let nidx = -k;
                maxn = maxn.max(nidx);
                let v = mom((nidx + 1) as i32).ok_or(format!("moment m = {} needed", nidx + 1))?;
                let mut t = Float::with_val(prec, g * &v[j]);
                if nidx % 2 == 1 {
                    t = -t;
                }
                tot += t;
            }
        } else {
            // a z-independent (pinched) master: the z^0 coefficient of c(z)
            if shift > 0 {
                continue;
            }
            let g = &ser[(-shift) as usize];
            let sgn = if ti.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
            let ev = efree.get(ti).ok_or(format!("value of the pinched master {ti:?} needed"))?;
            tot += Float::with_val(prec, g * ev) * sgn; // I = (-1)^(sum) E
        }
    }
    let sx = if x.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
    Ok((Float::with_val(prec, &tot * sx), maxn))
}
