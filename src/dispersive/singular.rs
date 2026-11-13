//! Singular points of the e-block in z = 1 + eta (oracle: values/d5l/r11/tools/fastsing.py).
//!
//! Every distinct denominator d(eta) is deflated exactly at eta = -1 (z = 0, always a singular point),
//! made primitive, and divided EXACTLY by the members of a growing factor base (mod-p prefilter, then exact
//! division over Q, repeated); a remainder of degree >= 1 joins the base. The union of the roots of the base
//! is the union of the roots of all denominators. Roots: Aberth-Ehrlich in MPFR complex arithmetic (precision
//! doubled until every root has a Newton correction below 2^-(prec/2) relative), then z = eta + 1.
use super::poly;
use rug::{Complex, Float, Integer};

const PREFILTER_PRIME: u64 = 2147483647; // 2^31 - 1

#[derive(Clone, Debug)]
pub struct SingularSet {
    /// singular points in z (including z = 0), deduplicated with tolerance 1e-9 (1 + |s|) as the oracle
    pub points: Vec<(f64, f64)>,
    /// the same points at high precision
    pub points_hp: Vec<Complex>,
    pub base: Vec<Vec<Integer>>,
    pub n_distinct_denominators: usize,
}

impl SingularSet {
    pub fn r0(&self) -> f64 {
        self.points.iter().map(|&(a, b)| a.hypot(b)).filter(|&r| r > 1e-9).fold(f64::INFINITY, f64::min)
    }
    pub fn rmax(&self) -> f64 {
        self.points.iter().map(|&(a, b)| a.hypot(b)).fold(0.0, f64::max)
    }
    /// singular points on the negative real axis (z < 0), descending (closest to 0 first)
    pub fn axis(&self) -> Vec<f64> {
        let mut v: Vec<f64> = self
            .points
            .iter()
            .filter(|&&(a, b)| a < -1e-9 && b.abs() < 1e-12 * (1.0 + a.hypot(b)))
            .map(|&(a, _)| a)
            .collect();
        v.sort_by(|a, b| b.partial_cmp(a).unwrap());
        v
    }
}

/// exact deflation of the roots at x = -1 (synthetic division by (x + 1))
fn deflate_minus_one(p: &[Integer]) -> Vec<Integer> {
    let mut p = p.to_vec();
    while p.len() > 1 {
        let n = p.len();
        let mut q = vec![Integer::new(); n - 1];
        let mut acc = Integer::new();
        for i in (1..n).rev() {
            acc = -acc + &p[i];
            q[i - 1] = acc.clone();
        }
        // remainder = p(-1) = -acc + p[0]
        let rem = Integer::from(&p[0] - &acc);
        if rem != 0 {
            break;
        }
        p = q;
    }
    p
}

fn primitive_signed(p: Vec<Integer>) -> Vec<Integer> {
    let q: Vec<rug::Rational> = p.into_iter().map(rug::Rational::from).collect();
    poly::ZPoly::from_q(&q).prim_signed()
}

pub fn factor_base(dens: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    let mut dens: Vec<Vec<Integer>> = dens.to_vec();
    dens.sort_by(|a, b| a.len().cmp(&b.len()).then(a.cmp(b)));
    let mut base: Vec<Vec<Integer>> = Vec::new();
    for d in dens {
        let p = deflate_minus_one(&d);
        if p.len() <= 1 {
            continue;
        }
        let mut p = primitive_signed(p);
        let mut moved = true;
        while p.len() > 1 && moved {
            moved = false;
            for f in &base {
                if f.len() <= p.len() && poly::divides_mod_p(&p, f, PREFILTER_PRIME) {
                    if let Some(q) = poly::exact_div(&p, f) {
                        p = q;
                        moved = true;
                        break;
                    }
                }
            }
        }
        if p.len() > 1 {
            base.push(p);
        }
    }
    base
}

/// all complex roots of an integer polynomial (ascending), Aberth-Ehrlich at increasing precision
pub fn roots(p: &[Integer], prec0: u32) -> Vec<Complex> {
    let deg = p.len() - 1;
    if deg == 1 {
        let r = Float::with_val(prec0, -rug::Rational::from((p[0].clone(), p[1].clone())));
        return vec![Complex::with_val(prec0, (r, 0))];
    }
    let mut prec = prec0;
    loop {
        if let Some(r) = aberth(p, prec, 400 + 40 * deg) {
            return r;
        }
        prec *= 2;
        assert!(prec < 1 << 16, "root finding does not converge");
    }
}

fn aberth(p: &[Integer], prec: u32, maxit: usize) -> Option<Vec<Complex>> {
    let deg = p.len() - 1;
    let c: Vec<Complex> = p.iter().map(|x| Complex::with_val(prec, (x, 0))).collect();
    let dc: Vec<Complex> = (1..=deg).map(|k| Complex::with_val(prec, &c[k] * k as u32)).collect();
    // Fujiwara bound for the initial radius
    let lead = Float::with_val(prec, c[deg].real());
    let mut rad = 0f64;
    for k in 0..deg {
        let v = Float::with_val(prec, c[k].real() / &lead).abs();
        if v.is_zero() {
            continue;
        }
        let r = (v.to_f64().ln() / (deg - k) as f64).exp();
        rad = rad.max(2.0 * r);
    }
    let rad = if rad.is_finite() && rad > 0.0 { rad } else { 1.0 };
    let mut z: Vec<Complex> = (0..deg)
        .map(|k| {
            let th = 2.0 * std::f64::consts::PI * (k as f64 + 0.25) / deg as f64 + 0.4;
            Complex::with_val(prec, (rad * th.cos(), rad * th.sin()))
        })
        .collect();
    let tol = Float::with_val(prec, Float::i_exp(1, -(prec as i32) / 2));
    let ev = |pc: &[Complex], x: &Complex| -> Complex {
        let mut acc = Complex::new(prec);
        for a in pc.iter().rev() {
            acc *= x;
            acc += a;
        }
        acc
    };
    for _ in 0..maxit {
        let mut done = true;
        let mut nz = z.clone();
        for i in 0..deg {
            let pv = ev(&c, &z[i]);
            let dv = ev(&dc, &z[i]);
            if dv.is_zero() {
                return None;
            }
            let w = Complex::with_val(prec, &pv / &dv);
            let mut s = Complex::new(prec);
            for j in 0..deg {
                if j != i {
                    let diff = Complex::with_val(prec, &z[i] - &z[j]);
                    if diff.is_zero() {
                        return None;
                    }
                    s += Complex::with_val(prec, 1 / diff);
                }
            }
            let den = Complex::with_val(prec, 1 - Complex::with_val(prec, &w * &s));
            let corr = Complex::with_val(prec, &w / &den);
            let rel = Float::with_val(prec, corr.clone().abs().real()) / Float::with_val(prec, z[i].clone().abs().real()).max(&Float::with_val(prec, 1e-300));
            if rel > tol {
                done = false;
            }
            nz[i] -= corr;
        }
        z = nz;
        if done {
            return Some(z);
        }
    }
    None
}

pub fn singular_points(dens: &[Vec<Integer>], prec: u32) -> SingularSet {
    let base = factor_base(dens);
    let mut sing: Vec<Complex> = vec![Complex::new(prec)];
    for f in &base {
        for r in roots(f, prec) {
            sing.push(r + 1u32);
        }
    }
    let mut pts: Vec<(f64, f64)> = Vec::new();
    let mut hp: Vec<Complex> = Vec::new();
    for s in sing {
        let (a, b) = (s.real().to_f64(), s.imag().to_f64());
        if pts.iter().all(|&(x, y)| (a - x).hypot(b - y) > 1e-9 * (1.0 + a.hypot(b))) {
            pts.push((a, b));
            hp.push(s);
        }
    }
    SingularSet { points: pts, points_hp: hp, base, n_distinct_denominators: dens.len() }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn zp(v: &[i64]) -> Vec<Integer> {
        v.iter().map(|&x| Integer::from(x)).collect()
    }
    #[test]
    fn factor_base_and_roots() {
        // eta-denominators: (1 + eta) (eta - 2), (eta - 2)^2 (eta^2 + 1), 3 (eta + 1/3)
        let d1 = zp(&[-2, -1, 1]);
        let d2 = zp(&[4, -4, 5, -4, 1]);
        let d3 = zp(&[1, 3]);
        let s = singular_points(&[d1, d2, d3], 128);
        // z = eta + 1: 0, 3, 1 +- i, 2/3
        let mut want = vec![(0.0, 0.0), (3.0, 0.0), (1.0, 1.0), (1.0, -1.0), (2.0 / 3.0, 0.0)];
        assert_eq!(s.points.len(), want.len(), "{:?}", s.points);
        for p in &s.points {
            let k = want.iter().position(|w| (w.0 - p.0).hypot(w.1 - p.1) < 1e-12).expect("unexpected point");
            want.remove(k);
        }
        assert_eq!(s.base.len(), 3);
    }
}
