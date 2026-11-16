//! Route (B) for a one-mass massless propagator family at a fixed rational d0 (method: docs/pflow.md):
//!   rho(s) = -(1/pi) Im J(-s + i0),  moments  M_m = FP int_0^inf rho(s) s^-m ds.
//! What the moments are (E-normalisation, D_e = k^2 + z Euclidean, a = a_e the index of the massive line, J = E[a](z)):
//!   * m = 1: the massless value E[a](z = 0);
//!   * m = n + 1 >= 2: the z^n Taylor coefficient of the analytic part times (-1)^n, i.e. (a)_n / n! E[a + n 1_e] (NOT the dotted
//!     integral itself unless a = 1);
//!   * m = -k <= 0: the z^(-k-1) coefficient of the integer part at infinity times (-1)^k, i.e. (-1)^(a-1) (a)_n / n! E[a_e -> -n]
//!     with n = k + 1 - a (the e-pinched integral with the massive line as a numerator of power n); zero when n < 0.
//! rho on (0, 1) = disc of the non-integer families at z = 0 (classes k h, k = 1..L), on (1, inf) = disc of the non-integer
//! families at infinity (classes -k h); both obey the homogeneous e-block. Two half-line transports (patch kernel, real
//! projection, weight (0 - z)^-m), local families at the threshold z = -1 (all classes j h mod 1, j = 0..2L), termwise end
//! pieces; head [0, s0] and tail [S, inf) termwise. Boundary directions not fixed by the region rows are fixed by the continuity
//! of the integer-family projection at z = -1; directions that continuity cannot see are tested for their effect on the outputs.
//! Every consistency check is an acceptance GATE (see `Acceptance`): a failed gate makes the run "not accepted".
//! Normalisation: the DE is RustRed's (Minkowski, D_e = k^2 - z); its solution vector is I = (-1)^(sum a) E with E the Euclidean
//! value (measure d^dk/pi^(d/2), p^2 = 1). Boundary rows are converted with that sign; outputs are given in E.
use super::graph::{self, Graph};
use super::local::{self, Pt};
use crate::dispersive::boundary::BoundaryRow;
use crate::dispersive::dispersion::{FamilyNum, log10abs};
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
    /// highest allowed lowest Frobenius order of the region families: the n_lo of a class is min(nlo, lowest nonzero boundary row - 2);
    /// the threshold families (no rows) use min(nlo, -10)
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
    /// radius of convergence of the masters' expansions at 0 and (in 1/z) at infinity: distance to the nearest PHYSICAL singular
    /// point (1 for one mass, p^2 = 1: the threshold z = -1); closer DE singular points are spurious (row precision budget)
    pub radius: f64,
    /// lower-loop massless values the graph evaluator cannot reduce (from a lower-loop run of this mode)
    pub provided: graph::Provided,
    /// work directory for checkpoints (Frobenius families every `frobenius_every_seconds`, every transport after every patch,
    /// finished moments as stage files); a rerun of the same command resumes. None: no checkpoints.
    pub workdir: Option<String>,
    pub frobenius_every_seconds: f64,
    /// acceptance: gates hold to 10^-(w - gate_margin) relative (fits, threshold decomposition, continuity, undetermined
    /// directions, exactly-vanishing integrals, closed-form masters); a failed gate makes the result "not accepted"
    pub gate_margin: u32,
    /// series-tail bound of the transports (log10); None = the kernel default -(working digits)/2
    pub max_tail_log10: Option<f64>,
    /// digits the caller needs: every gate asks for 10^-target (default: min(w, provided digits) - gate_margin)
    pub target_digits: Option<f64>,
    /// masters without boundary rows (numerators, irreducible lower-loop parts) fail the acceptance unless allowed
    pub allow_skipped: bool,
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
            radius: 1.0,
            provided: graph::Provided::default(),
            workdir: None,
            frobenius_every_seconds: 900.0,
            gate_margin: 20,
            target_digits: None,
            max_tail_log10: None,
            allow_skipped: false,
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

#[derive(Clone)]
pub struct Fam {
    pub side: &'static str,
    pub class: i64,
    pub num: std::sync::Arc<FamilyNum>,
    pub coef: Option<Vec<Float>>, // fitted scaled coefficients (rows (n - n_lo) E + comp)
    /// directions of the coefficient space not fixed by the boundary rows (fixed later by the threshold continuity)
    pub null: Vec<Vec<Float>>,
    pub diag: Value,
}

fn family(sys: &EBlockSystem, pt: &Pt, lam: &Rational, rho: &Rational, cfg: &PConfig, nlo: i64, mm: i64, name: &str) -> Result<Option<(FamilyNum, Value)>, String> {
    let e = sys.dim();
    let off = Rational::new();
    let bp = local::lmat::<ModP>(sys, &off, pt, mm + 8, 0);
    let out = frobenius::frobenius::<ModP>(&bp, e, &ModP::from_rational(lam, 0), nlo, mm, None, 0, cfg.threads);
    if out.ranks.iter().all(|&r| r == 0) || out.q == 0 {
        return Ok(None);
    }
    let spec = match frobenius::family_spec(&out, e, nlo, mm, lam, rho, name) {
        Ok(s) => s,
        // ranks zero up to the truncation edge: no solution of this class (the threshold dimension check guards completeness)
        Err(err) if err.contains("no rank plateau") && out.t[..((((mm - nlo) / 2 + 1) as usize) * e).min(out.t.len())].iter().all(|r| r.is_empty()) => return Ok(None),
        Err(err) => return Err(format!("{name}: {err}")),
    };
    if spec.dim == 0 {
        return Ok(None);
    }
    let prec = run::bits(cfg.wf);
    let bf = local::lmat::<Float>(sys, &off, pt, mm + 8, prec);
    let mode = PivotMode::parse(frobenius::DEFAULT_PIVOTING)?;
    let ck = cfg.workdir.as_ref().map(|w| crate::dispersive::ckpt::FrobCheckpoint {
        path: std::path::PathBuf::from(format!("{w}/frob-{}.bin", name.replace(' ', "_").replace('=', ""))),
        fingerprint: format!("{} {name} nlo {nlo} mm {mm} rho {rho} prec {prec} pivoting {}", system_fingerprint(sys), mode.label()),
        every_seconds: cfg.frobenius_every_seconds,
        stop_at_m: None,
    });
    let (num, diag) = run::family_num_ck(sys, &bf, &spec, mode, 8, prec, cfg.threads, ck.as_ref())?;
    Ok(Some((num, diag)))
}

/// boundary data of one master on one side: all-order exact rows (simple masters), the leading row of every region class, or the
/// reason it has none
/// reason it has none (`Skipped`: not computable here, gated unless --allow-skipped; `Failed`: a value that cannot be a zero, e.g. a
/// Gamma pole at a non-generic d0; always a failed gate)
enum MasterRows {
    Exact,
    Regions(Vec<(i64, i64, Float)>),
    Skipped(String),
    Failed(String),
}

fn master_rows(g: &Graph, e: usize, a: &[i64], h: &Rational, at_zero: bool, cfg: &PConfig, prec: u32) -> MasterRows {
    match graph::exact_series(g, e, a, h, at_zero, 1, 1, prec) {
        Ok(Some(_)) => return MasterRows::Exact,
        Err(err) => return MasterRows::Failed(err),
        Ok(None) => {}
    }
    match graph::regions(g, e, a, h, prec, at_zero, &cfg.provided) {
        Err(err) if err.starts_with(graph::GAMMA_POLE) => MasterRows::Failed(err),
        Err(err) => MasterRows::Skipped(err),
        Ok(None) => MasterRows::Skipped("a lower-loop part is neither reducible nor provided".into()),
        Ok(Some(rs)) => MasterRows::Regions(rs),
    }
}

/// nonzero rows (component, n, value; I-normalised) and all-zero masters of one class
struct ClassRows {
    nonzero: Vec<(usize, i64, Float)>,
    zero: Vec<usize>,
    /// rows that could not be computed (an exact series that fails at this class, NaN values): never read as zero
    errors: Vec<String>,
}

fn class_rows(sys: &EBlockSystem, g: &Graph, e: usize, h: &Rational, at_zero: bool, k: i64, data: &[MasterRows], cfg: &PConfig, prec: u32) -> ClassRows {
    let mut out = ClassRows { nonzero: vec![], zero: vec![], errors: vec![] };
    for (j, (a, d)) in sys.masters.iter().zip(data).enumerate() {
        let sgn = if a.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
        let mut push = |n: i64, v: Float| out.nonzero.push((j, n, Float::with_val(prec, &v * sgn)));
        match d {
            MasterRows::Exact => {
                let ex = match graph::exact_series(g, e, a, h, at_zero, k, cfg.exact_rows, prec) {
                    Ok(Some(ex)) => ex,
                    Ok(None) => {
                        out.errors.push(format!("{a:?}: exact series not available at class {k}"));
                        continue;
                    }
                    Err(err) => {
                        out.errors.push(err);
                        continue;
                    }
                };
                if ex.iter().any(|x| x.1.is_nan()) {
                    out.errors.push(format!("{a:?}: NaN in the exact series at class {k}"));
                    continue;
                }
                if ex.iter().all(|x| x.1.is_zero()) {
                    out.zero.push(j);
                } else {
                    for (n, v) in ex {
                        push(n, v);
                    }
                }
            }
            MasterRows::Regions(rs) => {
                let mine: Vec<&(i64, i64, Float)> = rs.iter().filter(|x| x.0 == k).collect();
                if mine.is_empty() {
                    out.zero.push(j);
                } else {
                    // only the LOWEST order is clean
                    let nmin = mine.iter().map(|x| x.1).min().unwrap();
                    let mut v = Float::new(prec);
                    for x in mine.iter().filter(|x| x.1 == nmin) {
                        v += &x.2;
                    }
                    push(nmin, v);
                }
            }
            MasterRows::Skipped(_) | MasterRows::Failed(_) => {}
        }
    }
    out
}

/// the fit rows of one family: the nonzero rows, zero rows below each master's first row (from n_lo), zero rows at every order
/// for the masters without a region in this class; errors for nonzero rows the family cannot represent
fn boundary_rows(cr: &ClassRows, fam: &FamilyNum, e: usize, prec: u32) -> (Vec<BoundaryRow>, Vec<String>, usize) {
    let n_lo = fam.spec.n_lo;
    let nrow = (fam.basis.len() / e) as i64;
    let mut out = Vec::new();
    let mut errs = Vec::new();
    // a master without a region of this class has an identically vanishing component in it: zero rows at EVERY order of the
    // family (the precision budget cuts the top); padding only the lowest orders left higher-order solutions unconstrained
    for j in &cr.zero {
        for n in n_lo..n_lo + nrow {
            out.push(BoundaryRow { comp: *j, n, value: Float::new(prec) });
        }
    }
    let mut first: std::collections::BTreeMap<usize, i64> = Default::default();
    for (j, n, _) in &cr.nonzero {
        let f = first.entry(*j).or_insert(*n);
        *f = (*f).min(*n);
    }
    for (j, f) in &first {
        for n in n_lo..(*f).min(n_lo + nrow) {
            out.push(BoundaryRow { comp: *j, n, value: Float::new(prec) });
        }
    }
    let mut beyond = 0usize;
    for (j, n, v) in &cr.nonzero {
        if *n < n_lo {
            errs.push(format!("{}: nonzero row of master {j} at n = {n} below n_lo = {n_lo}", fam.spec.name));
            continue;
        }
        if *n >= n_lo + nrow {
            // higher orders of an all-order (exact) series than the family carries: not usable, counted; losing a master's
            // LEADING row is an error
            if first.get(j) == Some(n) {
                errs.push(format!("{}: the leading row of master {j} (n = {n}) lies beyond the family's {nrow} orders", fam.spec.name));
            }
            beyond += 1;
            continue;
        }
        out.push(BoundaryRow { comp: *j, n: *n, value: v.clone() });
    }
    (out, errs, beyond)
}

/// rank-revealing least squares: Householder QR with column pivoting; columns whose remaining norm falls below
/// 2^(-prec/2) of the largest are left free (set to 0). Returns (x, rank).
pub fn lstsq_rr(a: &[Vec<Float>], b: &[Float], prec: u32) -> (Vec<Float>, usize) {
    let (x, rank, _) = lstsq_rr_null(a, b, prec);
    (x, rank)
}

/// as `lstsq_rr`, also returning a basis of the null space of `a` (free columns = unit, pivot columns from R11 w = -R12 e_j)
pub fn lstsq_rr_null(a: &[Vec<Float>], b: &[Float], prec: u32) -> (Vec<Float>, usize, Vec<Vec<Float>>) {
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
    let mut null = Vec::new();
    for f in rank..n {
        let mut w = vec![Float::new(prec); rank];
        for k in (0..rank).rev() {
            let mut s = Float::with_val(prec, -&r[k][f]);
            for j in k + 1..rank {
                s -= Float::with_val(prec, &r[k][j] * &w[j]);
            }
            w[k] = Float::with_val(prec, &s / &r[k][k]);
        }
        let mut v = vec![Float::new(prec); n];
        v[perm[f]] = Float::with_val(prec, 1);
        for k in 0..rank {
            v[perm[k]] = w[k].clone();
        }
        null.push(v);
    }
    (x, rank, null)
}

/// relative-weighted fit of the boundary rows on a family basis (as dispersion::soft_fit, rank revealing)
/// Highest row order a boundary fit can use. A singular point of the DE at distance `delta` from the expansion point, closer than
/// the solutions' own radius `radius` (a spurious point of the master basis: the masters are holomorphic there), makes the local
/// basis columns grow like delta^-n while the physical rows grow like radius^-n: an order-n row needs n log10(radius/delta) digits
/// of cancellation in the family rows. Rows beyond the budget wf - w - margin would make a consistent fit look inconsistent.
pub fn order_budget(nmin: i64, wf: u32, w: u32, radius: f64, delta: f64) -> (i64, f64) {
    let loss = (radius / delta).log10().max(0.0);
    let budget = wf as f64 - w as f64 - 10.0;
    if loss < 1e-9 {
        return (i64::MAX, loss);
    }
    (nmin + (budget / loss).floor().max(1.0) as i64, loss)
}

/// result of one boundary fit
struct Fit {
    coef: Vec<Float>,
    rows: usize,
    dropped_budget: usize,
    rank: usize,
    /// worst relative residual (log10) of an over-determined fit; None without rows
    worst: Option<f64>,
    null: Vec<Vec<Float>>,
}

/// `nmax`: rows of higher order are left out (precision budget, see `order_budget`)
fn fit_rows(fam: &FamilyNum, bnd: &[BoundaryRow], e: usize, prec: u32, nmax: i64) -> Fit {
    let rho = Float::with_val(prec, &fam.spec.rho);
    let mut amat = Vec::new();
    let mut bvec = Vec::new();
    let mut dropped = 0usize;
    for br in bnd {
        if br.n > nmax {
            dropped += 1;
            continue;
        }
        let row = ((br.n - fam.spec.n_lo) as usize) * e + br.comp;
        let sc = Float::with_val(prec, (&rho).pow(br.n as i32));
        let vs = Float::with_val(prec, &br.value * &sc);
        let w = if vs.is_zero() { Float::with_val(prec, 1) } else { Float::with_val(prec, 1 / Float::with_val(prec, vs.clone().abs())) };
        amat.push(fam.basis[row].iter().map(|x| Float::with_val(prec, x * &w)).collect::<Vec<_>>());
        bvec.push(Float::with_val(prec, &vs * &w));
    }
    let (x, rank, null) = lstsq_rr_null(&amat, &bvec, prec);
    let mut worst: Option<f64> = None;
    for (i, row) in amat.iter().enumerate() {
        let mut v = Float::with_val(prec, Float::dot(row.iter().zip(&x)));
        v -= &bvec[i];
        worst = Some(worst.unwrap_or(-1e9).max(log10abs(&v)));
    }
    let coef = fam.basis.iter().map(|row| Float::with_val(prec, Float::dot(row.iter().zip(&x)))).collect();
    let null = null.iter().map(|v| fam.basis.iter().map(|row| Float::with_val(prec, Float::dot(row.iter().zip(v)))).collect()).collect();
    Fit { coef, rows: amat.len(), dropped_budget: dropped, rank, worst, null }
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

/// log10 of the condition number of `a` with unit-norm columns: |R_00| / |R_rr| of a column-pivoted QR (rank r)
pub fn cond_log10(a: &[Vec<Float>], prec: u32) -> f64 {
    let m = a.len();
    if m == 0 {
        return 0.0;
    }
    let n = a[0].len();
    let norms: Vec<Float> = (0..n).map(|k| { let mut s = Float::new(prec); for row in a { s += Float::with_val(prec, row[k].square_ref()); } s.sqrt() }).collect();
    let mut r: Vec<Vec<Float>> = a.iter().map(|row| row.iter().zip(&norms).map(|(x, c)| if c.is_zero() { Float::new(prec) } else { Float::with_val(prec, x / c) }).collect()).collect();
    let tol_exp = -(prec as i32) / 2;
    let (mut first, mut last) = (None::<f64>, None::<f64>);
    for k in 0..n.min(m) {
        let (mut best, mut bv) = (k, Float::new(prec));
        for j in k..n {
            let mut c = Float::new(prec);
            for i in k..m {
                c += Float::with_val(prec, r[i][j].square_ref());
            }
            let c = c.sqrt();
            if c > bv {
                bv = c;
                best = j;
            }
        }
        if bv.is_zero() || bv.get_exp().unwrap_or(i32::MIN) < tol_exp {
            break;
        }
        let lv = log10abs(&bv);
        first.get_or_insert(lv);
        last = Some(lv);
        for row in r.iter_mut() {
            row.swap(k, best);
        }
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
    }
    match (first, last) { (Some(a), Some(b)) => a - b, _ => 0.0 }
}

/// threshold neighbourhood: decompose rho at s = 1 - r (side 0) / 1 + r (side inf) in the local families at z = -1, integrate
/// FP int rho s^-m ds to s = 1 termwise. Returns (integral, decomposition residual log10, integer-family coefficients)
fn threshold_piece(tf: &[FamilyNum], y: &[Float], r: &Rational, inf_side: bool, m: i32, jmax: u32, prec: u32) -> (Vec<Float>, f64, Vec<Float>, f64) {
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
    let cond = cond_log10(&a, prec);
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
            }
            let w = Float::with_val(prec, &base * &acc) * sgn_n;
            for i in 0..e {
                out[i] += Float::with_val(prec, &(Float::with_val(prec, pc * &f.basis[b * e + i][*t])) * &w);
            }
        }
    }
    (out, log10abs(&res) - ynorm, intc, cond)
}

/// real segment a -> b (rationals), steps r/ratio with r the distance to the nearest singular point; detours (half 12-gons,
/// radius 0.4 x gap) above axis points strictly inside
pub(crate) fn plan_segment(a: &Rational, b: &Rational, points: &[(f64, f64)], w: u32, ratio: f64) -> Vec<PatchSpec> {
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
    /// every acceptance gate passed (diag["acceptance"])
    pub accepted: bool,
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

/// exact dyadic rational with 24 significant bits (never 0)
fn rq(x: f64) -> Rational {
    let ex = x.log2().floor() as i32 - 23;
    let m = (x / 2f64.powi(ex)).floor() as i64;
    if ex >= 0 { Rational::from(m) * Rational::from(rug::Integer::from(1) << ex as u32) } else { Rational::from((m, 1)) / Rational::from(rug::Integer::from(1) << (-ex) as u32) }
}

pub struct Splits {
    /// distances of the nearest other singular point to 0 and to -1, largest |singular point|
    pub d0: f64,
    pub dthr: f64,
    pub rmax: f64,
    pub s0: Rational,
    pub r: Rational,
    pub s_end: Rational,
}

/// split points adapted to the singular points closest to 0, -1 and infinity (spurious points of the Laporta basis may sit at
/// distance O(eps) from the physical ones): head [0, s0] with s0 <= 0.4 d0, threshold neighbourhoods r <= 0.4 dthr, tail from
/// S >= 2.5 rmax; the families converge on the solutions' own radius, the transports detour
pub fn split_points(cfg: &PConfig, pts: &[(f64, f64)]) -> Splits {
    let dist = |c: (f64, f64)| pts.iter().filter(|p| (p.0 - c.0).hypot(p.1 - c.1) > 1e-12).map(|p| (p.0 - c.0).hypot(p.1 - c.1)).fold(f64::INFINITY, f64::min);
    let d0 = dist((0.0, 0.0));
    let dthr = dist((-1.0, 0.0));
    let rmax = pts.iter().map(|p| p.0.hypot(p.1)).fold(0.0f64, f64::max);
    let cap = |x: &Rational, lim: f64| -> Rational { if x.to_f64() > lim { rq(lim) } else { x.clone() } };
    let s_end = if cfg.s_end.to_f64() < 2.5 * rmax { rq(2.5 * rmax) } else { cfg.s_end.clone() };
    Splits { d0, dthr, rmax, s0: cap(&cfg.s0, 0.4 * d0), r: cap(&cfg.r, 0.4 * dthr), s_end }
}

/// identity of the DE for the family checkpoints: d0, masters and the CONTENT of every term (a corrected DE with the same shape
/// never reuses stale state); the boundary data are not part of it (the families do not depend on them)
pub fn system_fingerprint(sys: &EBlockSystem) -> String {
    let mut h = blake3::Hasher::new();
    h.update(format!("pflow e-block d0 {} e {} dim {}", sys.d0, sys.e_slot, sys.dim()).as_bytes());
    for m in &sys.masters {
        h.update(format!("{m:?};").as_bytes());
    }
    for t in &sys.terms {
        h.update(format!("{t:?};").as_bytes());
    }
    h.finalize().to_hex()[..24].to_string()
}

fn is_tail_failure(e: &crate::Error) -> bool {
    matches!(e, crate::Error::Numerical(m) if m.contains("series tail"))
}

/// the tail heal of the dispersive kernel (main e19645f2, `dispersive::kernel::TAIL_HEAL_POLICY`) for one pflow transport: the
/// patch `state.next_patch` that failed the series-tail bound is retried at the policy's higher Taylor orders (same centre, step,
/// precision and bound) with the kernel's own step on a one-patch path; Ok(true) = `state` advanced by the healed patch
fn heal_patch(tr: &Transport, state: &mut crate::patch_kernel::transport::TransportState) -> Result<bool, String> {
    let k = state.next_patch;
    let planned = tr.path.patches[k].order;
    for order in crate::dispersive::kernel::tail_heal_orders(planned) {
        let mut spec = tr.path.patches[k].clone();
        spec.order = order;
        let one = PlannedPath { offset: tr.path.offset.clone(), patches: vec![spec] };
        let t1 = Transport::new(tr.system, &one, tr.config.clone()).map_err(|x| x.to_string())?;
        let mut st = state.clone();
        st.next_patch = 0;
        match t1.step(&mut st, &Default::default()) {
            Ok(()) => {
                st.records.last_mut().expect("record").index = k;
                st.next_patch = k + 1;
                eprintln!("[pflow] tail heal: patch {} order {planned} -> {order}: accepted", k + 1);
                *state = st;
                return Ok(true);
            }
            Err(e) if is_tail_failure(&e) => eprintln!("[pflow] tail heal: patch {} order {planned} -> {order}: {e}", k + 1),
            Err(e) => return Err(e.to_string()),
        }
    }
    Ok(false)
}

/// one patch-kernel transport, with the tail heal; with a work directory it checkpoints after every patch (and after a healed
/// patch) and resumes from `<workdir>/tr-<tag>.json`. `heals` counts the healed patches (one counter per run, so concurrent runs in
/// one process do not mix)
pub(crate) fn transport_ck(tr: &Transport, y0: &[Complex], cfg: &PConfig, tag: &str, heals: &std::sync::atomic::AtomicUsize) -> Result<crate::patch_kernel::transport::TransportState, String> {
    let file = cfg.workdir.as_ref().map(|w| std::path::PathBuf::from(format!("{w}/tr-{tag}.json")));
    let mut st = match &file {
        Some(f) if f.exists() => {
            let st = tr.resume(f, y0, false).map_err(|e| format!("{}: {e} (use a fresh workdir)", f.display()))?;
            eprintln!("[pflow] transport {tag}: resumed at patch {}", st.next_patch);
            st
        }
        _ => tr.initial_state(y0).map_err(|x| x.to_string())?,
    };
    let policy = file.as_ref().map(|f| crate::patch_kernel::checkpoint::CheckpointPolicy { every_patches: 1, every_seconds: 300.0, fresh: false, ..crate::patch_kernel::checkpoint::CheckpointPolicy::new(f) });
    loop {
        match tr.run(&mut st, usize::MAX, y0, &Default::default(), policy.as_ref(), |_| {}) {
            Ok(_) => return Ok(st),
            Err(e) if is_tail_failure(&e) && st.next_patch < tr.path.patches.len() => {
                let k = st.next_patch;
                eprintln!("[pflow] transport {tag}: patch {}: {e}; tail heal at orders {:?}", k + 1, crate::dispersive::kernel::tail_heal_orders(tr.path.patches[k].order));
                if !heal_patch(tr, &mut st)? {
                    return Err(format!("transport {tag}: {e}; tail heal failed ({})", crate::dispersive::kernel::TAIL_HEAL_POLICY));
                }
                heals.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if let Some(f) = &file {
                    let _lock = crate::patch_kernel::checkpoint::CheckpointLock::acquire(f).map_err(|e| e.to_string())?;
                    tr.write_checkpoint(f, &st, y0).map_err(|e| e.to_string())?;
                }
            }
            Err(e) => return Err(e.to_string()),
        }
    }
}

/// diagnostics of one final transport (M3): the largest overlap discrepancy, series tail and dropped imaginary part over its patches
/// (log10, relative), the patches whose overlap probe was skipped, and the healed patches (order above the planned one)
fn transport_summary(st: &crate::patch_kernel::transport::TransportState) -> Value {
    let worst = |f: &dyn Fn(&crate::patch_kernel::transport::PatchRecord) -> f64| {
        st.records.iter().map(f).fold(None, |m: Option<f64>, x| if x.is_nan() { Some(f64::NAN) } else { Some(m.map_or(x, |m| if m.is_nan() { m } else { m.max(x) })) })
    };
    let nz = |x: Option<f64>| match x { Some(v) if v.is_finite() => json!(v), Some(v) => json!(v.to_string()), None => Value::Null };
    json!({"patches": st.records.len(), "max_overlap_log10": nz(worst(&|r| r.overlap_log10)), "max_tail_log10": nz(worst(&|r| r.tail_log10)),
        "max_dropped_imag_log10": nz(worst(&|r| r.dropped_imag_log10)), "overlap_skipped": st.records.iter().filter(|r| r.overlap_skipped).count()})
}

/// records the largest error (log10) of one kind of gate quantity
fn note_err(m: &mut std::collections::BTreeMap<&'static str, f64>, kind: &'static str, x: f64) {
    if x.is_nan() {
        m.insert(kind, 1e300); // a NaN gate quantity: no digits
        return;
    }
    let e = m.entry(kind).or_insert(f64::NEG_INFINITY);
    *e = e.max(x);
}

/// acceptance record: every gate with its verdict
#[derive(Default)]
struct Acceptance {
    failures: Vec<String>,
    checks: Vec<Value>,
}

impl Acceptance {
    fn check(&mut self, name: &str, ok: bool, detail: String) {
        self.checks.push(json!({"gate": name, "passed": ok, "detail": detail}));
        if !ok {
            self.failures.push(format!("{name}: {detail}"));
        }
    }
}

/// route (B) at one d0. The result carries `accepted` (all gates passed) and the gate record in diag["acceptance"].
pub fn run_b(sys: &EBlockSystem, g: &Graph, e_slot: usize, cfg: &PConfig) -> Result<PResult, String> {
    let t0 = Instant::now();
    let heals = std::sync::atomic::AtomicUsize::new(0);
    let e = sys.dim();
    let prec = run::bits(cfg.wf);
    let h = Rational::from(&sys.d0 / 2u32);
    let pts = singular_points(sys);
    let sp = split_points(cfg, &pts);
    let (d0, dthr, rmax) = (sp.d0, sp.dthr, sp.rmax);
    let rho0 = rq(0.5 * d0);
    let rho_inf = rq(0.5 / rmax);
    let rho_thr = rq(0.5 * dthr);
    let mut cfg = cfg.clone();
    (cfg.s0, cfg.r, cfg.s_end) = (sp.s0, sp.r, sp.s_end);
    let cfg = &cfg;
    // the accuracy the gates ask for (relative, log10): the transports carry w digits, the provided values their own
    let digits_cap = (cfg.w as f64).min(cfg.provided.digits.unwrap_or(f64::INFINITY));
    let gate = -cfg.target_digits.unwrap_or(digits_cap - cfg.gate_margin as f64);
    let mut acc = Acceptance::default();
    // every gate quantity that measures an error (log10, relative); digits_estimate is the minimum over them
    let mut errq: std::collections::BTreeMap<&'static str, f64> = std::collections::BTreeMap::new();
    let mut diag = json!({"s0": cfg.s0.to_string(), "r": cfg.r.to_string(), "S": cfg.s_end.to_string(), "singular_points": pts.iter().map(|p| vec![p.0, p.1]).collect::<Vec<_>>(), "rho0": rho0.to_string(), "rho_inf": rho_inf.to_string(), "rho_thr": rho_thr.to_string(), "families": [], "gate_log10": gate, "digits_cap": digits_cap});
    // ---- boundary data of every master, both sides (once) ----
    let mut data = Vec::new();
    for at_zero in [true, false] {
        let d: Vec<MasterRows> = sys.masters.iter().map(|a| master_rows(g, e_slot, a, &h, at_zero, cfg, prec)).collect();
        let skipped: Vec<String> = d.iter().zip(&sys.masters).filter_map(|(x, a)| if let MasterRows::Skipped(why) = x { Some(format!("{a:?}: {why}")) } else { None }).collect();
        acc.check(&format!("masters with boundary rows ({})", if at_zero { "0" } else { "inf" }), skipped.is_empty() || cfg.allow_skipped, format!("{} skipped: {:?}", skipped.len(), skipped));
        let failed: Vec<String> = d.iter().zip(&sys.masters).filter_map(|(x, a)| if let MasterRows::Failed(why) = x { Some(format!("{a:?}: {why}")) } else { None }).collect();
        acc.check(&format!("boundary rows computable ({})", if at_zero { "0" } else { "inf" }), failed.is_empty(), format!("{} failed: {:?}", failed.len(), failed));
        data.push(d);
    }
    let mut fams0 = Vec::new();
    let mut famsi = Vec::new();
    for k in 1..=cfg.loops as i64 {
        for (side, pt, lam, rho, mm) in [("0", Pt::Finite(Rational::new()), Rational::from(&h * k), &rho0, cfg.m0), ("inf", Pt::Infinity, -Rational::from(&h * k), &rho_inf, cfg.minf)] {
            let t = Instant::now();
            let cr = class_rows(sys, g, e_slot, &h, side == "0", k, &data[if side == "0" { 0 } else { 1 }], cfg, prec);
            // n_lo from the rows: the family must contain the solution that starts at the lowest nonzero row
            acc.check(&format!("class rows computable {side} k={k}"), cr.errors.is_empty(), cr.errors.join("; "));
            let nmin_nz = cr.nonzero.iter().map(|x| x.1).min();
            let nlo = cfg.nlo.min(nmin_nz.map_or(cfg.nlo, |n| n - 2));
            let mm_c = mm; // absolute top order: a lower n_lo adds orders below
            let Some((num, fdiag)) = family(sys, &pt, &lam, rho, cfg, nlo, mm_c, &format!("{side} k={k}"))? else {
                acc.check(&format!("family {side} k={k}"), cr.nonzero.is_empty(), format!("empty class but {} nonzero boundary rows", cr.nonzero.len()));
                continue;
            };
            let (br, errs, beyond) = boundary_rows(&cr, &num, e, prec);
            acc.check(&format!("rows representable {side} k={k}"), errs.is_empty(), errs.join("; "));
            // precision budget counted from the first PHYSICAL (nonzero) row
            let delta = if side == "0" { d0 } else { 1.0 / rmax };
            let (nmax, loss) = match nmin_nz {
                Some(n) => order_budget(n, cfg.wf, cfg.w, cfg.radius, delta),
                None => (i64::MAX, 0.0),
            };
            if loss > 1e-9 {
                acc.check(&format!("order budget {side} k={k}"), cfg.wf as f64 - cfg.w as f64 - 10.0 > loss, format!("wf - w - 10 = {} for {loss:.2} digits/order", cfg.wf as i64 - cfg.w as i64 - 10));
            }
            let mut fit = fit_rows(&num, &br, e, prec, nmax);
            if fit.rank < num.spec.dim && fit.dropped_budget > 0 {
                let all = fit_rows(&num, &br, e, prec, i64::MAX);
                acc.check(&format!("order budget keeps the rank {side} k={k}"), all.rank <= fit.rank,
                    format!("rank {} with the budget (n <= {nmax}), {} with all rows: raise wf", fit.rank, all.rank));
            }
            if let Some(wst) = fit.worst {
                if fit.rows > fit.rank {
                    acc.check(&format!("boundary fit {side} k={k}"), wst < gate, format!("worst relative residual 10^{wst:.1}"));
                    note_err(&mut errq, "boundary_fits", wst);
                }
            }
            diag["families"].as_array_mut().unwrap().push(json!({"side": side, "class": k, "dim": num.spec.dim, "n_lo": num.spec.n_lo, "rows": fit.rows,
                "rows_offered": br.len(), "rows_beyond_family_orders": beyond, "rows_dropped_budget": fit.dropped_budget, "row_order_max": if nmax == i64::MAX { Value::Null } else { json!(nmax) },
                "digits_lost_per_order": loss, "rank": fit.rank, "rank_deficit": num.spec.dim - fit.rank, "nonzero_rows": cr.nonzero.len(),
                "fit_worst_rel_log10": fit.worst, "frobenius": fdiag, "seconds": t.elapsed().as_secs_f64()}));
            eprintln!("[pflow] family {side} k={k}: dim {} rank {} n_lo {} fit rows {} of {} ({} dropped by the budget n <= {nmax}, {loss:.2} digits/order) worst {:?} ({:.1} s)",
                num.spec.dim, fit.rank, num.spec.n_lo, fit.rows, br.len(), fit.dropped_budget, fit.worst.map(|x| (x * 10.0).round() / 10.0), t.elapsed().as_secs_f64());
            let null = std::mem::take(&mut fit.null);
            let f = Fam { side: if side == "0" { "0" } else { "inf" }, class: k, num: std::sync::Arc::new(num), coef: Some(fit.coef), null, diag: Value::Null };
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
        // no rows set n_lo at the threshold: at most -10 (the span gate below checks completeness)
        if let Some((num, _)) = family(sys, &Pt::Finite(Rational::from(-1)), &lam, &rho_thr, cfg, cfg.nlo.min(-10), cfg.mthr, &format!("thr {j}"))? {
            tfams.push(num);
        }
    }
    let tdim: usize = tfams.iter().map(|f| f.spec.dim).sum();
    diag["threshold_dims"] = json!(tfams.iter().map(|f| (f.spec.lam.to_string(), f.spec.dim)).collect::<Vec<_>>());
    acc.check("threshold families span the e-block", tdim == e, format!("{tdim} of {e}"));
    // binomial terms of the threshold pieces: (w + 10) / log10(1/r) + |m| + 10
    let jmax = |m: i32| -> u32 { ((cfg.w as f64 + 10.0) / (1.0 / cfg.r.to_f64()).log10()).ceil() as u32 + m.unsigned_abs() + 10 };
    // transports
    let ksys = crate::dispersive::kernel::rational_system(sys)?;
    let tprec = run::bits(cfg.w);
    let one = Rational::from(1);
    let p0 = plan_segment(&Rational::from(-&cfg.s0), &Rational::from(&cfg.r - &one), &pts, cfg.w, cfg.ratio);
    let pi_ = plan_segment(&Rational::from(-&cfg.s_end), &Rational::from(-(Rational::from(&one + &cfg.r))), &pts, cfg.w, cfg.ratio);
    let s0f = Float::with_val(prec, &cfg.s0);
    let sef = Float::with_val(prec, &cfg.s_end);
    let allzero = |y: &[Float]| y.iter().all(|x| x.is_zero());
    // ---- the moments of a given set of family coefficients (also used for the undetermined directions) ----
    let moments_of = |f0: &[Fam], fi: &[Fam], tag: &str, pieces: bool| -> Result<(Vec<(i32, Vec<Float>)>, Value), String> {
        let mut out = Vec::new();
        let mut d = json!({});
        let y0 = rho_at(f0, e, &s0f, prec);
        let yi = rho_at(fi, e, &sef, prec);
        for &m in &cfg.moments {
            let mut tot = end_piece(f0, e, &cfg.s0, m, prec);
            let tl = end_piece(fi, e, &cfg.s_end, m, prec);
            for (half, y, patches) in [(0, &y0, &p0), (1, &yi, &pi_)] {
                if allzero(y) {
                    continue;
                }
                let mut tc = TransportConfig::dispersive(tprec);
                tc.cauchy_pole = Some(Rational::new());
                tc.cauchy_power = m;
                if let Some(b) = cfg.max_tail_log10 {
                    tc.max_tail_log10 = b;
                }
                tc.threads = cfg.threads;
                tc.pin_threads = false;
                let path = PlannedPath { offset: Rational::new(), patches: patches.clone() };
                let tr = Transport::new(&ksys, &path, tc).map_err(|x| x.to_string())?;
                let yc: Vec<Complex> = y.iter().map(|v| Complex::with_val(tprec, v)).collect();
                let st = transport_ck(&tr, &yc, cfg, &format!("{tag}m{m}-h{half}"), &heals)?;
                if pieces {
                    d[format!("transport_m{m}_h{half}")] = transport_summary(&st);
                }
                let yend: Vec<Float> = st.y.iter().map(|c| Float::with_val(prec, c.real())).collect();
                let integ = st.integral.as_ref().unwrap();
                // kernel: int y (0 - z)^-m dz ; on z = -s: = -int rho s^-m ds (0-side, s0 -> 1-r) and = +int_{1+r}^{S} (inf-side)
                let mut mid = Vec::new();
                for i in 0..e {
                    let v = Float::with_val(prec, integ[i].real());
                    mid.push(if half == 0 { Float::with_val(prec, -&v) } else { v.clone() });
                    if half == 0 { tot[i] -= v } else { tot[i] += v }
                }
                let (th, res, ic, cond) = threshold_piece(&tfams, &yend, &cfg.r, half == 1, m, jmax(m), prec);
                for i in 0..e {
                    tot[i] += &th[i];
                }
                if pieces {
                    d[format!("threshold_residual_{half}_m{m}")] = json!(res);
                    d[format!("threshold_condition_{half}_log10")] = json!(cond);
                    if m == cfg.moments[0] {
                        let f = |v: &Vec<Float>| v.iter().map(|x| x.to_string_radix(10, Some(20))).collect::<Vec<_>>();
                        d[format!("pieces_mid{half}")] = json!(f(&mid));
                        d[format!("pieces_thr{half}")] = json!(f(&th));
                        d[format!("pieces_yend{half}")] = json!(f(&yend));
                        d[format!("intcoef_{half}")] = json!(ic.iter().map(|x| x.to_string_radix(10, None)).collect::<Vec<_>>());
                    }
                }
            }
            for i in 0..e {
                tot[i] += &tl[i];
            }
            if pieces && m == cfg.moments[0] {
                let f = |v: &Vec<Float>| v.iter().map(|x| x.to_string_radix(10, Some(20))).collect::<Vec<_>>();
                d["pieces_m1"] = json!({"head": f(&end_piece(f0, e, &cfg.s0, m, prec)), "tail": f(&tl), "total": f(&tot)});
            }
            out.push((m, tot));
        }
        Ok((out, d))
    };
    // ---- free boundary directions: fixed by the continuity of the integer-family projection at z = -1 ----
    let ndir0: usize = fams0.iter().map(|f| f.null.len()).sum();
    let ndiri: usize = famsi.iter().map(|f| f.null.len()).sum();
    // directions (combinations of family null vectors) that nothing determines: (side, family index, coefficient vector)
    let mut undetermined: Vec<Vec<(usize, usize, Vec<Float>)>> = Vec::new();
    if ndir0 + ndiri > 0 {
        let t = Instant::now();
        let transport_end = |y: &[Float], patches: &[PatchSpec], tag: String| -> Result<Vec<Float>, String> {
            let mut tc = TransportConfig::dispersive(tprec);
            tc.cauchy_pole = None;
            if let Some(b) = cfg.max_tail_log10 {
                tc.max_tail_log10 = b;
            }
            tc.threads = cfg.threads;
            tc.pin_threads = false;
            let path = PlannedPath { offset: Rational::new(), patches: patches.to_vec() };
            let tr = Transport::new(&ksys, &path, tc).map_err(|x| x.to_string())?;
            let yc: Vec<Complex> = y.iter().map(|v| Complex::with_val(tprec, v)).collect();
            let st = transport_ck(&tr, &yc, cfg, &format!("cont-{tag}"), &heals)?;
            Ok(st.y.iter().map(|c| Float::with_val(prec, c.real())).collect())
        };
        // integer-family projection and the decomposition residual (log10, relative)
        let proj = |y: &[Float], inf_side: bool| -> (Vec<Float>, f64) {
            let (_, res, ic, _) = threshold_piece(&tfams, y, &cfg.r, inf_side, 1, 1, prec);
            (ic, res)
        };
        let single = |fams: &[Fam], fi: usize, c: &[Float], s: &Float| -> Vec<Float> {
            let f = &fams[fi];
            let tmp = [Fam { side: f.side, class: f.class, num: f.num.clone(), coef: Some(c.to_vec()), null: vec![], diag: Value::Null }];
            rho_at(&tmp, e, s, prec)
        };
        // noise floor of a projection: the transports carry w digits (|noise| ~ 10^-(w-8) |transported vector|), and the
        // decomposition at the threshold adds its own residual
        let ymax = |y: &[Float]| y.iter().map(log10abs).fold(-1e9f64, f64::max);
        let floor = |y: &[Float], res: f64| (ymax(y) - (cfg.w as f64 - 8.0)).max(ymax(y) + res);
        let y0e = transport_end(&rho_at(&fams0, e, &s0f, prec), &p0, "base0".into())?;
        let yie = transport_end(&rho_at(&famsi, e, &sef, prec), &pi_, "basei".into())?;
        let (b0, r0) = proj(&y0e, false);
        let (bi, ri) = proj(&yie, true);
        let base_floor = floor(&y0e, r0).max(floor(&yie, ri));
        let mut cols: Vec<Vec<Float>> = Vec::new();
        let mut which: Vec<(usize, usize, usize)> = Vec::new(); // (side, family index, direction)
        let mut dropped = 0usize;
        let mut cols_detail = Vec::new();
        for (side, fams, s, patches) in [(0usize, &fams0, &s0f, &p0), (1, &famsi, &sef, &pi_)] {
            for (fi, f) in fams.iter().enumerate() {
                for (di, c) in f.null.iter().enumerate() {
                    let y = transport_end(&single(fams, fi, c, s), patches, format!("dir{side}-{fi}-{di}"))?;
                    let (pr, res) = proj(&y, side == 1);
                    let (col, fl) = (ymax(&pr), floor(&y, res));
                    cols_detail.push(json!({"side": side, "family": fi, "direction": di, "column_log10": col, "floor_log10": fl}));
                    // a direction whose continuity column is at the noise floor is not seen by continuity: undetermined (tested below)
                    if col < fl + 3.0 {
                        dropped += 1;
                        undetermined.push(vec![(side, fi, c.clone())]);
                        continue;
                    }
                    cols.push(if side == 0 { pr } else { pr.iter().map(|x| Float::with_val(prec, -x)).collect() });
                    which.push((side, fi, di));
                }
            }
        }
        let neq = b0.len();
        let rhs: Vec<Float> = (0..neq).map(|i| Float::with_val(prec, &bi[i] - &b0[i])).collect();
        let base_resid = ymax(&rhs);
        // nothing to fix if the base solution is already continuous to the noise floor (the visible directions are then 0)
        let active = !cols.is_empty() && base_resid > base_floor + 3.0;
        let (alpha, rank) = if active {
            // unit-norm columns, so that the rank decision compares like with like
            let norms: Vec<Float> = cols.iter().map(|c| { let mut n = Float::new(prec); for x in c { n += Float::with_val(prec, x.square_ref()); } n.sqrt() }).collect();
            let a: Vec<Vec<Float>> = (0..neq).map(|i| cols.iter().zip(&norms).map(|(c, n)| Float::with_val(prec, &c[i] / n)).collect()).collect();
            let (al, rank, nullc) = lstsq_rr_null(&a, &rhs, prec);
            // combinations of visible directions that continuity does not fix: undetermined as well
            for v in &nullc {
                let mut comb = Vec::new();
                for (k, ((side, fi, di), n)) in which.iter().zip(&norms).enumerate() {
                    if v[k].is_zero() {
                        continue;
                    }
                    let f = if *side == 0 { &fams0[*fi] } else { &famsi[*fi] };
                    let sc = Float::with_val(prec, &v[k] / n);
                    comb.push((*side, *fi, f.null[*di].iter().map(|x| Float::with_val(prec, x * &sc)).collect::<Vec<_>>()));
                }
                undetermined.push(comb);
            }
            (al.iter().zip(&norms).map(|(x, n)| Float::with_val(prec, x / n)).collect::<Vec<Float>>(), rank)
        } else {
            (vec![Float::new(prec); cols.len()], cols.len())
        };
        let mut worst = -1e9f64;
        let scale = rhs.iter().chain(b0.iter()).map(log10abs).fold(-1e9f64, f64::max);
        for i in 0..neq {
            let mut v = Float::with_val(prec, -&rhs[i]);
            for (c, al) in cols.iter().zip(&alpha) {
                v += Float::with_val(prec, &c[i] * al);
            }
            worst = worst.max(log10abs(&v) - scale);
        }
        if active {
            acc.check("continuity fit residual", worst < gate, format!("10^{worst:.1}"));
            note_err(&mut errq, "continuity_fit_residual", worst);
        }
        for ((side, fi, di), al) in which.iter().zip(&alpha) {
            let f = if *side == 0 { &mut fams0[*fi] } else { &mut famsi[*fi] };
            let d = f.null[*di].clone();
            let c = f.coef.as_mut().unwrap();
            for (x, y) in c.iter_mut().zip(&d) {
                *x += Float::with_val(prec, al * y);
            }
        }
        diag["continuity_fit_detail"] = json!({"base_residual_log10": base_resid, "noise_floor_log10": base_floor, "applied": active, "dropped_noise_directions": dropped, "columns": cols_detail});
        diag["continuity_fit"] = json!({"free_directions": [ndir0, ndiri], "equations": neq, "active_columns": cols.len(), "rank": rank, "worst_rel_residual_log10": worst, "seconds": t.elapsed().as_secs_f64()});
        eprintln!("[pflow] continuity fit: {} + {} free directions ({dropped} at the noise floor), {} equations, rank {rank} of {} columns, base residual {base_resid:.1} vs floor {base_floor:.1}, applied {active}, worst residual {worst:.1} ({:.1} s)", ndir0, ndiri, neq, cols.len(), t.elapsed().as_secs_f64());
    }
    diag["patches"] = json!([p0.len(), pi_.len()]);
    eprintln!("[pflow] singular points {:?}; patches {} / {}", pts, p0.len(), pi_.len());
    let (moments, pd) = moments_of(&fams0, &famsi, "", true)?;
    if let Value::Object(m) = pd {
        for (k, v) in m {
            diag[k] = v;
        }
    }
    for &m in &cfg.moments {
        for half in 0..2 {
            if let Some(r) = diag[format!("threshold_residual_{half}_m{m}")].as_f64() {
                acc.check(&format!("threshold decomposition half {half} m {m}"), r < gate, format!("10^{r:.1}"));
                note_err(&mut errq, "threshold_decomposition", r);
            }
        }
    }
    // continuity of the integer-family projection at z = -1
    let mut cont = -1e9f64;
    if let (Some(a), Some(b)) = (diag["intcoef_0"].as_array(), diag["intcoef_1"].as_array()) {
        for (x, y) in a.iter().zip(b) {
            let (x, y) = (Float::with_val(prec, Float::parse(x.as_str().unwrap()).unwrap()), Float::with_val(prec, Float::parse(y.as_str().unwrap()).unwrap()));
            let d = log10abs(&Float::with_val(prec, &x - &y)) - log10abs(&x).max(log10abs(&y));
            cont = cont.max(d);
        }
        // an ACCURACY ESTIMATE of the threshold join (measured at 4L: no1 at w45 vs w55 agrees to the continuity level, while the
        // decomposition's condition number, reported for information, is a far too pessimistic bound)
        let cond = diag["threshold_condition_0_log10"].as_f64().unwrap_or(0.0).max(diag["threshold_condition_1_log10"].as_f64().unwrap_or(0.0));
        acc.check("threshold continuity", cont < gate, format!("10^{cont:.1} (gate 10^{gate:.1}; decomposition condition 10^{cond:.1})"));
        note_err(&mut errq, "threshold_continuity", cont);
    }
    diag["threshold_continuity_log10"] = json!(cont);
    // ---- undetermined directions: their effect on every requested output must be below the gate ----
    let scale_of = |v: &[Float]| v.iter().map(log10abs).fold(-1e9f64, f64::max);
    let mut und_detail = Vec::new();
    for (u, comb) in undetermined.iter().enumerate() {
        let mut f0: Vec<Fam> = fams0.iter().map(|f| Fam { coef: Some(vec![Float::new(prec); f.num.basis.len()]), null: vec![], ..f.clone() }).collect();
        let mut fi: Vec<Fam> = famsi.iter().map(|f| Fam { coef: Some(vec![Float::new(prec); f.num.basis.len()]), null: vec![], ..f.clone() }).collect();
        for (side, k, c) in comb {
            let f = if *side == 0 { &mut f0[*k] } else { &mut fi[*k] };
            let cf = f.coef.as_mut().unwrap();
            for (x, y) in cf.iter_mut().zip(c) {
                *x += y;
            }
        }
        let (mu, _) = moments_of(&f0, &fi, &format!("und{u}-"), false)?;
        // contribution relative to the output scale of each moment (per master: relative to max(|value|, 10^-30 scale))
        let mut worst = -1e9f64;
        let mut affected = Vec::new();
        for ((m, cu), (_, v)) in mu.iter().zip(&moments) {
            let sc = scale_of(v);
            for (j, (x, y)) in cu.iter().zip(v).enumerate() {
                let r = log10abs(x) - log10abs(y).max(sc - 30.0);
                worst = worst.max(r);
                if r > gate {
                    affected.push(json!({"m": m, "master": sys.masters[j], "relative_log10": r}));
                }
            }
        }
        acc.check(&format!("undetermined direction {u}"), affected.is_empty(), format!("largest relative effect 10^{worst:.1} on {} outputs", affected.len()));
        note_err(&mut errq, "undetermined_directions", worst);
        und_detail.push(json!({"direction": u, "largest_relative_effect_log10": worst, "affected": affected}));
    }
    diag["undetermined_directions"] = json!(und_detail);
    // ---- independent checks: exactly vanishing integrals and closed-form masters (the graph evaluator, no provided values) ----
    if let Some(k1) = moments.iter().position(|(m, _)| *m == 1) {
        let v = &moments[k1].1;
        let sc = scale_of(v);
        let none = graph::Provided::default();
        let (mut nz, mut nc, mut worst_zero, mut worst_cf) = (0usize, 0usize, -1e9f64, -1e9f64);
        let mut bad = Vec::new();
        for (j, a) in sys.masters.iter().enumerate() {
            let Some(ex) = graph::massless(g, a, &h, prec, &none) else { continue };
            let sgn = if a.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
            let x = Float::with_val(prec, &v[j] * sgn);
            if ex.is_zero() {
                nz += 1;
                let r = log10abs(&x) - sc;
                worst_zero = worst_zero.max(r);
                if r > gate {
                    bad.push(format!("{a:?} should vanish: 10^{r:.1} of the scale"));
                }
            } else {
                nc += 1;
                let r = log10abs(&Float::with_val(prec, &x - &ex)) - log10abs(&ex);
                worst_cf = worst_cf.max(r);
                if r > gate {
                    bad.push(format!("{a:?}: 10^{r:.1} from its closed form"));
                }
            }
        }
        acc.check("vanishing and closed-form masters", bad.is_empty(), format!("{nz} vanishing (worst 10^{worst_zero:.1}), {nc} closed forms (worst 10^{worst_cf:.1}); {}", bad.join("; ")));
        diag["independent_checks"] = json!({"vanishing": nz, "closed_forms": nc, "worst_vanishing_log10": worst_zero, "worst_closed_form_log10": worst_cf});
        if nz > 0 {
            note_err(&mut errq, "vanishing_masters", worst_zero);
        }
        if nc > 0 {
            note_err(&mut errq, "closed_form_masters", worst_cf);
        }
    }
    // ---- the provided table: every entry the evaluator can reduce itself agrees; every irreducible graph is provided ----
    let missing = cfg.provided.missing.lock().unwrap().clone();
    acc.check("irreducible lower-loop parts provided", missing.is_empty(), format!("{} missing: {:?}", missing.len(), missing));
    let mism = cfg.provided.mismatch.lock().unwrap().clone();
    let mlist: Vec<String> = mism.iter().map(|(k, d, n)| format!("{k}: {d} ({n} uses)")).collect();
    acc.check("provided values = evaluator where it can reduce", mism.is_empty(), format!("{} mismatching entries: {:?}", mism.len(), mlist));
    diag["provided_values"] = json!(cfg.provided.values.len());
    diag["provided_digits"] = json!(cfg.provided.digits);
    diag["tail_heals"] = json!(heals.load(std::sync::atomic::Ordering::SeqCst));
    diag["tail_heal_policy"] = json!(crate::dispersive::kernel::TAIL_HEAL_POLICY);
    // M3: the worst transport diagnostics over the final transports (information; the tail is gated by the kernel and healed)
    let mut tmax = json!({});
    for key in ["max_overlap_log10", "max_tail_log10", "max_dropped_imag_log10"] {
        let w = diag.as_object().unwrap().iter().filter(|(k, _)| k.starts_with("transport_m")).filter_map(|(_, v)| v[key].as_f64()).fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x))));
        tmax[key] = json!(w);
    }
    diag["transports_worst"] = tmax;
    // the digits a run can claim: the precision cap and the minimum over every gate quantity that measures an error
    let mut est = digits_cap;
    let mut by = json!({"digits_cap": digits_cap});
    for (k, x) in &errq {
        by[*k] = json!(-x);
        est = est.min(-x);
    }
    diag["digits_estimate"] = json!(est);
    diag["digits_by_check"] = by;
    diag["seconds"] = json!(t0.elapsed().as_secs_f64());
    let accepted = acc.failures.is_empty();
    diag["acceptance"] = json!({"accepted": accepted, "failures": acc.failures, "checks": acc.checks});
    if !accepted {
        eprintln!("[pflow] NOT ACCEPTED: {}", diag["acceptance"]["failures"]);
    }
    Ok(PResult { moments, diag, accepted })
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
