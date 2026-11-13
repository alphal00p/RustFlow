//! Spectral gluing: one Kallen-Lehmann bubble glued into an auxiliary-mass line.
//!
//! A vacuum integral whose graph contains a one-loop bubble (two unit-mass
//! lines between two vertices) equals a dispersion integral over the bubble's
//! spectral density times the remaining graph with one line of mass^2 `s`.
//! The remaining graph is an auxiliary-mass flow target with `z = 1 + eta`;
//! [`crate::PreparedFlow::evaluate_spectral`] supplies its values on a contour
//! and its `eta -> infinity` expansions, and [`gluing`] assembles the glued
//! integrals.
//!
//! Conventions (Euclidean, measure `d^d l / pi^(d/2)`, unit masses,
//! `d = 4 - 2 eps`):
//!
//! * `Bub_E(Q^2) = int_4^oo ds rho(s) / (Q^2 + s)` and
//!   `B_E(z) := Bub_E(Q^2 = -z) = Gamma(eps) 2F1(1, eps; 3/2; z/4)`.
//! * `int_4^oo rho(s) V(s) ds = (1/2 pi i) int_Gamma B_E(z) V(z) dz` with a
//!   V-shaped contour `Gamma` through `z = c` (default 5/2), arms `c +- i t`.
//!   For real `eps` only the lower arm is needed: with
//!   `L' = int_c^{c - i oo} W(z) V(z) dz` the integral is `-Im(L') / pi`.
//! * The arm is integrated with Gauss-Legendre panels on `t in [0, T0]` and
//!   analytically beyond `T0` from the large-`z` expansions (terms `z^beta`
//!   integrate to `-z0^(beta+1)/(beta+1)`, continued in `d`).
//!
//! Numerics use fixed-precision MPFR/MPC through `rug` directly. Numerica's
//! `Float` tracks precision loss dynamically, which is unsuitable for long
//! quadrature sums whose terms cancel.
use rug::float::Constant;
use rug::ops::Pow;
use std::collections::BTreeMap;

use crate::{Error, Result};

#[cfg(feature = "automatic")]
pub mod gluing;

pub type R = rug::Float;
pub type C = rug::Complex;

/// Fixed binary precision of every gluing quantity.
#[derive(Clone, Copy, Debug)]
pub struct Ctx {
    pub bits: u32,
}

impl Ctx {
    /// `d` decimal digits plus 24 guard bits.
    pub fn digits(d: u32) -> Self {
        Ctx {
            bits: (d as f64 * std::f64::consts::LOG2_10).ceil() as u32 + 24,
        }
    }
    pub fn r(&self, v: i64) -> R {
        R::with_val(self.bits, v)
    }
    pub fn f(&self, v: f64) -> R {
        R::with_val(self.bits, v)
    }
    pub fn q(&self, n: i64, d: i64) -> R {
        self.r(n) / self.r(d)
    }
    /// Exact rational `NUM[/DEN]` rounded to the working precision.
    pub fn rational(&self, text: &str) -> Result<R> {
        let parse = |s: &str| {
            rug::Integer::from_str_radix(s.trim(), 10)
                .map_err(|e| Error::InvalidInput(format!("rational {text:?}: {e}")))
        };
        let value = match text.split_once('/') {
            Some((n, d)) => rug::Rational::from((parse(n)?, parse(d)?)),
            None => rug::Rational::from(parse(text)?),
        };
        Ok(R::with_val(self.bits, value))
    }
    pub fn cc(&self, re: R, im: R) -> C {
        C::with_val(self.bits, (re, im))
    }
    pub fn cr(&self, x: &R) -> C {
        C::with_val(self.bits, (x, 0))
    }
    pub fn ci(&self, re: i64, im: i64) -> C {
        C::with_val(self.bits, (re, im))
    }
    pub fn pi(&self) -> R {
        R::with_val(self.bits, Constant::Pi)
    }
    pub fn euler(&self) -> R {
        R::with_val(self.bits, Constant::Euler)
    }
    pub fn gamma(&self, x: &R) -> R {
        R::with_val(self.bits, x).gamma()
    }
    /// `2^-(bits + 8)`: relative series truncation threshold.
    pub fn tiny(&self) -> R {
        R::with_val(self.bits, 1) >> (self.bits + 8)
    }
}

pub fn abs(z: &C) -> R {
    R::with_val(z.prec().0, z.abs_ref())
}
pub fn re(z: &C) -> R {
    z.real().clone()
}
pub fn im(z: &C) -> R {
    z.imag().clone()
}
pub fn cpow(z: &C, e: &C) -> C {
    z.clone().pow(e)
}
pub fn creal_pow(ctx: &Ctx, z: &C, e: &R) -> C {
    cpow(z, &ctx.cr(e))
}
/// Full-precision decimal string.
pub fn fmt(x: &R) -> String {
    x.to_string_radix(10, None)
}
pub fn fmt_digits(x: &R, d: usize) -> String {
    x.to_string_radix(10, Some(d))
}

/// Significant digits of agreement `-log10(|a - b| / |b|)` (999 when equal).
pub fn agree(a: &C, b: &C) -> f64 {
    let d = abs(&(a.clone() - b));
    if d.is_zero() {
        return 999.0;
    }
    log10(&abs(b)) - log10(&d)
}
/// Real counterpart of [`agree`].
pub fn agree_r(a: &R, b: &R) -> f64 {
    let d = (a.clone() - b).abs();
    if d.is_zero() {
        return 999.0;
    }
    log10(&b.clone().abs()) - log10(&d)
}
/// `log10 |x|` without overflow of `f64` (`-1e9` for zero).
pub fn log10(x: &R) -> f64 {
    if x.is_zero() {
        return -1e9;
    }
    let (m, e) = x.to_f64_exp();
    m.abs().log10() + e as f64 * std::f64::consts::LOG10_2
}

/// `sum_{(k, n)} c z^(n + k eps)` with integers `k`, `n` and a fixed real `eps`.
#[derive(Clone, Debug, Default)]
pub struct Asym {
    pub terms: BTreeMap<(i64, i64), C>,
}

impl Asym {
    pub fn new() -> Self {
        Asym {
            terms: BTreeMap::new(),
        }
    }
    pub fn add_term(&mut self, k: i64, n: i64, c: C) {
        match self.terms.get_mut(&(k, n)) {
            Some(e) => *e += &c,
            None => {
                self.terms.insert((k, n), c);
            }
        }
    }
    /// Product, truncated below `z^nmin` in the integer power.
    pub fn mul(&self, o: &Asym, nmin: i64) -> Asym {
        let mut out = Asym::new();
        for ((k1, n1), c1) in &self.terms {
            for ((k2, n2), c2) in &o.terms {
                if n1 + n2 >= nmin {
                    out.add_term(k1 + k2, n1 + n2, c1.clone() * c2);
                }
            }
        }
        out
    }
    pub fn scale(&self, s: &C) -> Asym {
        Asym {
            terms: self
                .terms
                .iter()
                .map(|(k, c)| (*k, c.clone() * s))
                .collect(),
        }
    }
    pub fn add(&self, o: &Asym) -> Asym {
        let mut out = self.clone();
        for ((k, n), c) in &o.terms {
            out.add_term(*k, *n, c.clone());
        }
        out
    }
    pub fn eval(&self, ctx: &Ctx, eps: &R, z: &C) -> C {
        let lz = z.clone().ln();
        let mut s = ctx.ci(0, 0);
        for ((k, n), c) in &self.terms {
            let e = ctx.cr(&(ctx.r(*n) + eps.clone() * *k));
            s += c.clone() * (lz.clone() * &e).exp();
        }
        s
    }
    /// `int_{z0}^{z0 + oo dir} z^beta dz = -z0^(beta+1)/(beta+1)` (analytic
    /// continuation), summed termwise; for `beta = -1` exactly (`k = 0`,
    /// `n = -1`) the finite part `-ln(z0)`.
    ///
    /// The integer-power terms `P(z)` (`k = 0`, `n >= -1`) are analytic right
    /// of the contour, so their contour plus arc-at-infinity integral vanishes;
    /// the eps-regulated terms have no arc contribution. With real coefficients
    /// the finite parts make `P`'s lower-arm contribution to `Im(L')` vanish
    /// (`n >= 0`: `-c^(n+1)/(n+1)` real; `n = -1`: `-c ln c` real), i.e. the
    /// arc is included.
    pub fn tail_integral(&self, ctx: &Ctx, eps: &R, z0: &C) -> C {
        let lz = z0.clone().ln();
        let mut s = ctx.ci(0, 0);
        let mut zk: BTreeMap<i64, C> = BTreeMap::new();
        for ((k, n), c) in &self.terms {
            if *k == 0 && *n == -1 {
                s -= c.clone() * &lz;
                continue;
            }
            let b1 = ctx.r(*n + 1) + eps.clone() * *k;
            let pk = zk
                .entry(*k)
                .or_insert_with(|| (lz.clone() * ctx.cr(&(eps.clone() * *k))).exp())
                .clone();
            let pn = z0.clone().pow(*n as i32 + 1);
            s -= c.clone() * pk * pn / ctx.cr(&b1);
        }
        s
    }
}

struct Chart {
    w0: C,
    rad: R,
    coef: Vec<C>,
}

/// `B_E(z) = Gamma(eps) f(z/4)` with `f = 2F1(1, eps; 3/2; w)`: direct series
/// for `|w| <= 0.6`, Taylor charts of the hypergeometric equation along the
/// lower arm, the connection formula for `|w| >= 2`.
pub struct Bubble {
    pub ctx: Ctx,
    pub eps: R,
    pub gamma_eps: R,
    charts: Vec<Chart>,
    order: usize,
    /// Agreement of the last chart with the connection formula at
    /// `w = c/4 - 5i/2` (an internal consistency check).
    pub self_check_digits: f64,
}

impl Bubble {
    pub fn new(ctx: Ctx, eps: &R, c: &R) -> Result<Self> {
        let gamma_eps = ctx.gamma(eps);
        let order = ctx.bits as usize + 40;
        let mut b = Bubble {
            ctx,
            eps: eps.clone(),
            gamma_eps,
            charts: vec![],
            order,
            self_check_digits: 999.0,
        };
        // Charts along the lower arm w = c/4 - i y until |w| > 3 (the
        // connection formula is used beyond |w| = 2).
        let mut w = ctx.cr(&(c.clone() / 4));
        let (mut f0, mut f1) = b.direct(&w)?;
        loop {
            let rad = b.dist_sing(&w) / 2i32;
            let coef = b.taylor(&w, &f0, &f1);
            let next = ctx.cc(re(&w), im(&w) - &rad);
            let h = next.clone() - &w;
            let (g0, g1) = eval_taylor(&ctx, &coef, &h);
            b.charts.push(Chart {
                w0: w.clone(),
                rad: rad.clone(),
                coef,
            });
            if abs(&w) > 3 {
                break;
            }
            w = next;
            f0 = g0;
            f1 = g1;
        }
        let wt = ctx.cc(c.clone() / 4, -ctx.q(5, 2));
        let a = b
            .f_chart(&wt)
            .ok_or_else(|| Error::Numerical("bubble: check point not on a chart".into()))?;
        let bb = b.f_connection(&wt)?;
        b.self_check_digits = agree(&a, &bb);
        Ok(b)
    }
    fn dist_sing(&self, w: &C) -> R {
        let d0 = abs(w);
        let d1 = abs(&(w.clone() - 1));
        if d0 < d1 { d0 } else { d1 }
    }
    /// Direct series at `|w| < 1`: `(f, f')`.
    pub fn direct(&self, w: &C) -> Result<(C, C)> {
        let ctx = &self.ctx;
        let mut t = ctx.ci(1, 0);
        let mut f = ctx.ci(1, 0);
        let mut fp = ctx.ci(0, 0);
        let tiny = ctx.tiny();
        let mut n = 0i64;
        loop {
            let ratio = (self.eps.clone() + n) / (ctx.q(3, 2) + n);
            let tn1 = t.clone() * &ratio * w;
            n += 1;
            fp += tn1.clone() * n / w;
            f += &tn1;
            t = tn1;
            if abs(&t) < tiny.clone() * abs(&f) && n > 10 {
                break;
            }
            if n >= 400_000 {
                return Err(Error::Numerical(
                    "bubble: direct hypergeometric series did not converge".into(),
                ));
            }
        }
        Ok((f, fp))
    }
    fn taylor(&self, w0: &C, f0: &C, f1: &C) -> Vec<C> {
        let ctx = &self.ctx;
        let p0 = w0.clone() * (ctx.ci(1, 0) - w0);
        let p1 = ctx.ci(1, 0) - w0.clone() * 2i32;
        let ab1 = self.eps.clone() + 2i32;
        let q0: C = ctx.cr(&ctx.q(3, 2)) - C::with_val(ctx.bits, w0 * &ab1);
        let ab = self.eps.clone();
        let mut a = vec![f0.clone(), f1.clone()];
        for n in 0..self.order {
            let nn = n as i64;
            let t1 = (p1.clone() * (nn * (nn + 1)) + q0.clone() * (nn + 1)) * &a[n + 1];
            let t0 = a[n].clone() * (ctx.r(-(nn * (nn - 1))) - ab1.clone() * nn - &ab);
            let den = p0.clone() * ((nn + 1) * (nn + 2));
            a.push(-(t1 + t0) / den);
        }
        a
    }
    fn f_chart(&self, w: &C) -> Option<C> {
        let mut best: Option<(R, usize)> = None;
        for (i, ch) in self.charts.iter().enumerate() {
            let d = abs(&(w.clone() - &ch.w0));
            if d <= ch.rad.clone() * 1.01f64 {
                let rel = d / &ch.rad;
                if best.as_ref().is_none_or(|(b, _)| rel < *b) {
                    best = Some((rel, i));
                }
            }
        }
        let (_, i) = best?;
        let ch = &self.charts[i];
        Some(eval_taylor(&self.ctx, &ch.coef, &(w.clone() - &ch.w0)).0)
    }
    /// Connection formula about `w = infinity` (valid for `|w| > 1` off the cut).
    pub fn f_connection(&self, w: &C) -> Result<C> {
        let ctx = &self.ctx;
        let e = &self.eps;
        let u = ctx.ci(1, 0) / w;
        let mut t = ctx.ci(1, 0);
        let mut s = ctx.ci(1, 0);
        let tiny = ctx.tiny();
        let mut n = 0i64;
        loop {
            let ratio = (ctx.q(1, 2) + n) / (ctx.r(2) - e + n);
            t = t * &ratio * &u;
            s += &t;
            n += 1;
            if abs(&t) < tiny && n > 10 {
                break;
            }
            if n >= 200_000 {
                return Err(Error::Numerical(
                    "bubble: connection series did not converge".into(),
                ));
            }
        }
        let mw = -w.clone();
        let soft = (ctx.ci(1, 0) / &mw) * &s / (ctx.r(2) * (e.clone() - 1));
        let kappa =
            ctx.gamma(&ctx.q(3, 2)) * ctx.gamma(&(ctx.r(1) - e)) / ctx.gamma(&(ctx.q(3, 2) - e));
        let hard = cpow(&mw, &ctx.cr(&(-e.clone())))
            * cpow(&(ctx.ci(1, 0) - &u), &ctx.cr(&(ctx.q(1, 2) - e)));
        Ok(soft + hard * &kappa)
    }
    pub fn f(&self, w: &C) -> Result<C> {
        let aw = abs(w);
        if aw >= 2 && !(w.imag().is_zero() && *w.real() > 0) {
            return self.f_connection(w);
        }
        if aw <= 0.6f64 {
            return Ok(self.direct(w)?.0);
        }
        self.f_chart(w).ok_or_else(|| {
            Error::Numerical(format!(
                "bubble: point w = {} {:+}i is not covered by a chart",
                w.real().to_f64(),
                w.imag().to_f64()
            ))
        })
    }
    /// `B_E(z) = Bub_E(Q^2 = -z)`.
    pub fn b(&self, z: &C) -> Result<C> {
        Ok(self.f(&(z.clone() / 4))? * &self.gamma_eps)
    }
    /// `B21_E(z) = [Gamma(eps) - (1 - 2 eps) B_E(z)] / (4 - z)` (one line squared).
    pub fn b21(&self, z: &C) -> Result<C> {
        let one_m = self.ctx.r(1) - self.eps.clone() * 2;
        Ok((self.ctx.cr(&self.gamma_eps) - self.b(z)? * &one_m) / (self.ctx.ci(4, 0) - z))
    }
    /// Taylor coefficients of `B_E` around `z = 1`: `B_E(1 + h) = sum_m beta_m h^m`.
    pub fn taylor_at_one(&self, m: usize) -> Result<Vec<C>> {
        let ctx = &self.ctx;
        let w0 = ctx.cr(&ctx.q(1, 4));
        let (f0, f1) = self.direct(&w0)?;
        let a = self.taylor(&w0, &f0, &f1);
        let mut out = vec![];
        let mut s = ctx.r(1);
        for k in 0..m {
            out.push(a[k].clone() * (s.clone() * &self.gamma_eps));
            s /= 4;
        }
        Ok(out)
    }
    /// Large-`z` expansion of `B_E` on the LOWER half plane: soft terms
    /// `z^(-1-n)` and hard terms `z^(-eps-j)` for `n, j <= nmax`.
    pub fn asym(&self, nmax: i64) -> Asym {
        let ctx = &self.ctx;
        let e = &self.eps;
        let mut out = Asym::new();
        let mut s = ctx.r(1);
        let mut p4 = ctx.r(1);
        let pref: R = self.gamma_eps.clone() * -4i32 / (ctx.r(2) * (e.clone() - 1i32));
        for n in 0..=nmax {
            out.add_term(0, -1 - n, ctx.cr(&(pref.clone() * &p4 * &s)));
            s = s * (ctx.q(1, 2) + n) / (ctx.r(2) - e + n);
            p4 *= 4;
        }
        let kappa =
            ctx.gamma(&ctx.q(3, 2)) * ctx.gamma(&(ctx.r(1) - e)) / ctx.gamma(&(ctx.q(3, 2) - e));
        let four_eps = (ctx.r(4).ln() * e).exp();
        let pe = ctx.pi() * e;
        let ph = ctx.cc(pe.clone().cos(), -pe.sin());
        let base = ph * (self.gamma_eps.clone() * kappa * four_eps);
        let a = ctx.q(1, 2) - e;
        let mut bin = ctx.r(1);
        let mut m4 = ctx.r(1);
        for j in 0..=nmax {
            out.add_term(-1, -j, base.clone() * (bin.clone() * &m4));
            bin = bin * (a.clone() - j) / (j + 1);
            m4 *= -4;
        }
        out
    }
    /// Expansion of `B21_E = -(1/z) sum_k (4/z)^k [Gamma(eps) - (1 - 2 eps) B_E]`.
    pub fn asym21(&self, nmax: i64) -> Asym {
        let ctx = &self.ctx;
        let mut inner = self
            .asym(nmax)
            .scale(&ctx.cr(&(-(ctx.r(1) - self.eps.clone() * 2i32))));
        inner.add_term(0, 0, ctx.cr(&self.gamma_eps));
        let mut geo = Asym::new();
        let mut p = ctx.r(-1);
        for k in 0..=nmax {
            geo.add_term(0, -1 - k, ctx.cr(&p));
            p *= 4;
        }
        inner.mul(&geo, -nmax - 2)
    }
}

/// Value and derivative of `sum_n a_n h^n`.
pub fn eval_taylor(ctx: &Ctx, a: &[C], h: &C) -> (C, C) {
    let mut f = ctx.ci(0, 0);
    let mut fp = ctx.ci(0, 0);
    for (n, c) in a.iter().enumerate().rev() {
        if n >= 1 {
            fp = fp * h + c.clone() * n as i64;
        }
        f = f * h + c;
    }
    (f, fp)
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`.
pub fn gauss_legendre(ctx: &Ctx, n: usize) -> (Vec<R>, Vec<R>) {
    let mut xs = vec![];
    let mut ws = vec![];
    let tol = ctx.tiny() << 24u32;
    let leg = |x: &R| -> (R, R) {
        let (mut p0, mut p1) = (ctx.r(1), x.clone());
        for k in 1..n {
            let k_ = k as i64;
            let p2 = (x.clone() * &p1 * (2 * k_ + 1) - p0.clone() * k_) / (k_ + 1);
            p0 = p1;
            p1 = p2;
        }
        let dp = (x.clone() * &p1 - &p0) * n as i64 / (x.clone() * x - 1);
        (p1, dp)
    };
    for i in 1..=n {
        let x0 = (std::f64::consts::PI * (i as f64 - 0.25) / (n as f64 + 0.5)).cos();
        let mut x = ctx.f(x0);
        for _ in 0..100 {
            let (p, dp) = leg(&x);
            let dx = p / dp;
            x -= &dx;
            if dx.abs() < tol {
                break;
            }
        }
        let (_, dp) = leg(&x);
        ws.push(ctx.r(2) / ((ctx.r(1) - x.clone() * &x) * &dp * &dp));
        xs.push(x);
    }
    (xs, ws)
}

/// Lower-arm nodes `z_j = c - i t_j`, `t in [0, T0]`: geometric Gauss-Legendre
/// panels whose orders follow from the Bernstein ellipses of the singular
/// points `sing` (on the real `z` axis) and the requested `digits`.
pub struct Arm {
    pub c: R,
    pub t0: R,
    pub z: Vec<C>,
    pub w: Vec<R>,
    pub panels: Vec<(f64, f64, usize)>,
}

impl Arm {
    pub fn new(ctx: &Ctx, c: &R, t0: &R, sing: &[f64], digits: u32) -> Self {
        let cf = c.to_f64();
        let t0f = t0.to_f64();
        let mut edges = vec![0.0, 1.0];
        while *edges.last().unwrap() * 2.0 < t0f {
            let l = *edges.last().unwrap();
            edges.push(l * 2.0);
        }
        edges.push(t0f);
        let mut cache: BTreeMap<usize, (Vec<R>, Vec<R>)> = BTreeMap::new();
        let (mut z, mut w, mut panels) = (vec![], vec![], vec![]);
        for p in 0..edges.len() - 1 {
            let (a, b) = (edges[p], edges[p + 1]);
            let mut rho_min = f64::INFINITY;
            for &zs in sing {
                // t_s = i (z_s - c)
                let ts = NC(0.0, zs - cf);
                let u = ts.scale(2.0 / (b - a)).sub(NC((a + b) / (b - a), 0.0));
                let s = u.mul(u).sub(NC(1.0, 0.0)).sqrt();
                let r1 = u.add(s).abs();
                let r2 = u.sub(s).abs();
                rho_min = rho_min.min(r1.max(r2));
            }
            let n = ((digits as f64 + 6.0) * std::f64::consts::LN_10 / (2.0 * rho_min.ln())).ceil()
                as usize
                + 12;
            let (xs, ws) = cache
                .entry(n)
                .or_insert_with(|| gauss_legendre(ctx, n))
                .clone();
            let ar = ctx.f(a);
            let br = if p + 2 == edges.len() {
                t0.clone()
            } else {
                ctx.f(b)
            };
            let half = (br.clone() - &ar) / 2i32;
            let mid = (br + &ar) / 2i32;
            for (x, wt) in xs.iter().zip(ws.iter()) {
                let t = mid.clone() + half.clone() * x;
                z.push(ctx.cc(c.clone(), -t));
                w.push(wt.clone() * &half);
            }
            panels.push((a, b, n));
        }
        Arm {
            c: c.clone(),
            t0: t0.clone(),
            z,
            w,
            panels,
        }
    }
    pub fn z0(&self, ctx: &Ctx) -> C {
        ctx.cc(self.c.clone(), -self.t0.clone())
    }
    /// `L' = int_c^{c - i oo} F dz`: node sum plus the analytic tail of
    /// `asym`. Returns `(L', -Im(L') / pi)`.
    pub fn integrate(&self, ctx: &Ctx, eps: &R, f: &[C], asym: &Asym) -> (C, R) {
        let mut s = ctx.ci(0, 0);
        for (fj, wj) in f.iter().zip(&self.w) {
            s += fj.clone() * wj;
        }
        let num = s * ctx.ci(0, -1);
        let tail = asym.tail_integral(ctx, eps, &self.z0(ctx));
        let l = num + &tail;
        let i = -im(&l) / ctx.pi();
        (l, i)
    }
    /// As [`Arm::integrate`], with only the nodes at `t <= t_cut` (a panel
    /// edge) and the analytic tail from `c - i t_cut`: an independent split of
    /// numeric and asymptotic parts for the same integral.
    pub fn integrate_upto(&self, ctx: &Ctx, eps: &R, f: &[C], asym: &Asym, t_cut: &R) -> (C, R) {
        let mut s = ctx.ci(0, 0);
        for ((fj, wj), zj) in f.iter().zip(&self.w).zip(&self.z) {
            if -zj.imag().clone() <= *t_cut {
                s += fj.clone() * wj;
            }
        }
        let num = s * ctx.ci(0, -1);
        let z0 = ctx.cc(self.c.clone(), -t_cut.clone());
        let l = num + asym.tail_integral(ctx, eps, &z0);
        let i = -im(&l) / ctx.pi();
        (l, i)
    }
}

/// Double-precision complex arithmetic for panel sizing only.
#[derive(Clone, Copy)]
struct NC(f64, f64);
impl NC {
    fn add(self, o: NC) -> NC {
        NC(self.0 + o.0, self.1 + o.1)
    }
    fn sub(self, o: NC) -> NC {
        NC(self.0 - o.0, self.1 - o.1)
    }
    fn mul(self, o: NC) -> NC {
        NC(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn scale(self, s: f64) -> NC {
        NC(self.0 * s, self.1 * s)
    }
    fn abs(self) -> f64 {
        self.0.hypot(self.1)
    }
    fn sqrt(self) -> NC {
        let r = self.abs();
        let a = ((r + self.0) / 2.0).sqrt();
        let b = ((r - self.0) / 2.0).sqrt();
        NC(a, if self.1 < 0.0 { -b } else { b })
    }
}

/// `(z - 1)^p` as a large-`z` expansion, terms down to `z^nmin`.
pub fn zm1_pow(ctx: &Ctx, p: i64, nmin: i64) -> Asym {
    let mut a = Asym::new();
    // (z-1)^p = z^p sum_m binom(p, m) (-1)^m z^-m
    let mut c = ctx.r(1);
    let mut m = 0i64;
    while p - m >= nmin {
        if c.is_zero() && p >= 0 {
            break;
        }
        a.add_term(0, p - m, ctx.cr(&c));
        // binom(p, m+1)(-1)^(m+1) = binom(p, m)(-1)^m * (m - p)/(m+1)
        c = c * (m - p) / (m + 1);
        m += 1;
    }
    a
}

/// `eta -> infinity` expansion classes `sum_N s_N x^(base + slope eps + N)`,
/// `x = 1/eta = 1/(z - 1)`, re-expanded in `z` down to `z^nmin`.
pub fn eta_classes_to_z(ctx: &Ctx, eps: &R, classes: &[(i64, i64, Vec<C>)], nmin: i64) -> Asym {
    let mut out = Asym::new();
    for (base, slope, s) in classes {
        for (n, sn) in s.iter().enumerate() {
            if sn.is_zero() {
                continue;
            }
            // x^(beta+N) = (z-1)^-(beta+N) = z^-(beta+N) sum_m (beta+N)_m/m! z^-m
            let a = eps.clone() * *slope + (*base + n as i64);
            let lead = -(*base + n as i64);
            let mut c = ctx.cr(&ctx.r(1)) * sn;
            let mut m = 0i64;
            while lead - m >= nmin {
                out.add_term(-*slope, lead - m, c.clone());
                c = c * (a.clone() + m) / (m + 1);
                m += 1;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One-loop two-mass vacuum bubble with masses^2 1 and z:
    /// `V1(z) = Gamma(eps - 1) (1 - z^(1 - eps)) / (z - 1)`.
    fn v1_closed(ctx: &Ctx, eps: &R, z: &C) -> C {
        let g = ctx.gamma(&(eps.clone() - 1));
        (ctx.ci(1, 0) - creal_pow(ctx, z, &(ctx.r(1) - eps))) / (z.clone() - 1) * g
    }
    fn v1_asym(ctx: &Ctx, eps: &R, nmax: i64) -> Asym {
        let g = ctx.cr(&ctx.gamma(&(eps.clone() - 1)));
        let mut a = Asym::new();
        for n in 0..=nmax {
            a.add_term(0, -1 - n, g.clone());
            a.add_term(-1, -n, -g.clone());
        }
        a
    }

    /// The two-loop equal-mass sunrise as a spectral integral of the bubble
    /// against the two-mass one-loop vacuum bubble, compared with independent
    /// mpmath Feynman-parameter values (90 digits; campaign gate 1).
    #[test]
    fn sunrise_from_spectral_contour_matches_feynman_parameters() -> Result<()> {
        let references = [
            (
                "1/7",
                "-100.986582190664374837449567568300878086748535238417713322968197269610510430894078075745343",
            ),
            (
                "-1/5",
                "-25.7934906552508689394367929092488371549333951391814549331747094954233216432212502357227812",
            ),
            (
                "1/10000",
                "-150027688.787777318298684455766403471829309700668213179272765434262258907700876699362154374",
            ),
        ];
        let ctx = Ctx::digits(70);
        let c = ctx.q(5, 2);
        for (e, reference) in references {
            let eps = ctx.rational(e)?;
            let bubble = Bubble::new(ctx, &eps, &c)?;
            assert!(
                bubble.self_check_digits > 60.0,
                "{}",
                bubble.self_check_digits
            );
            let asym = bubble.asym(90).mul(&v1_asym(&ctx, &eps, 90), -92);
            let arm = Arm::new(&ctx, &c, &ctx.r(120), &[0.0, 4.0], 66);
            let f = arm
                .z
                .iter()
                .map(|z| Ok(bubble.b(z)? * v1_closed(&ctx, &eps, z)))
                .collect::<Result<Vec<_>>>()?;
            let (_, value) = arm.integrate(&ctx, &eps, &f, &asym);
            let (_, split) = arm.integrate_upto(&ctx, &eps, &f, &asym, &ctx.r(64));
            let reference = R::with_val(ctx.bits, R::parse(reference).unwrap());
            assert!(
                agree_r(&value, &reference) > 60.0,
                "eps {e}: {}",
                agree_r(&value, &reference)
            );
            assert!(agree_r(&split, &value) > 60.0, "eps {e}: tail split");
        }
        Ok(())
    }

    #[test]
    fn bubble_expansion_matches_connection_formula_far_out() -> Result<()> {
        let ctx = Ctx::digits(50);
        let eps = ctx.rational("3/10")?;
        let bubble = Bubble::new(ctx, &eps, &ctx.q(5, 2))?;
        let z = ctx.cc(ctx.q(5, 2), ctx.r(-4000));
        let direct = bubble.b(&z)?;
        let expanded = bubble.asym(40).eval(&ctx, &eps, &z);
        assert!(
            agree(&expanded, &direct) > 45.0,
            "{}",
            agree(&expanded, &direct)
        );
        let direct21 = bubble.b21(&z)?;
        let expanded21 = bubble.asym21(40).eval(&ctx, &eps, &z);
        assert!(agree(&expanded21, &direct21) > 45.0);
        Ok(())
    }

    #[test]
    fn binomial_and_eta_class_expansions_are_consistent() -> Result<()> {
        let ctx = Ctx::digits(40);
        let eps = ctx.rational("1/7")?;
        let z = ctx.cc(ctx.r(50), ctx.r(-30));
        let zm1: C = z.clone() - 1i32;
        let closed: C = zm1.pow(-3i32);
        let expanded = zm1_pow(&ctx, -3, -80).eval(&ctx, &eps, &z);
        assert!(agree(&expanded, &closed) > 35.0);
        // x^(1 + 2 eps) * (1 + 2 x), x = 1/(z - 1): classes base 1, slope 2.
        let classes = vec![(1, 2, vec![ctx.ci(1, 0), ctx.ci(2, 0)])];
        let expanded = eta_classes_to_z(&ctx, &eps, &classes, -80).eval(&ctx, &eps, &z);
        let x = ctx.ci(1, 0) / (z.clone() - 1);
        let power = ctx.r(1) + eps.clone() * 2;
        let closed = creal_pow(&ctx, &x, &power) * (ctx.ci(1, 0) + x * 2i32);
        assert!(
            agree(&expanded, &closed) > 35.0,
            "{}",
            agree(&expanded, &closed)
        );
        Ok(())
    }
}
