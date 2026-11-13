//! Path planning on z = -s, s in [s0, S] (oracle: prep5.py for rho0 / rho_inf / s0 / S, v9c zplan.py for the patches).
//!
//! Steps h = R / ratio (R = distance to the nearest singular point, z = 0 included) rounded to 24 significant
//! bits (exact dyadic), polygon (12-gon) detour over the top of every singular point on the negative axis with
//! radius 0.4 x (distance to the nearest other singular point or to 0), Taylor order nord = ceil((W + 8) /
//! log10(R / |h|)) + 4. Float operations mirror the oracle so that the same dyadic path results.
use super::poly;
use rug::Rational;

#[derive(Clone, Debug)]
pub struct Patch {
    pub centre: (Rational, Rational),
    pub step: (Rational, Rational),
    pub order: usize,
}

impl Patch {
    pub fn is_real(&self) -> bool {
        self.centre.1 == 0 && self.step.1 == 0
    }
}

/// x rounded to `den` grid (python round = ties to even)
pub fn ratround(x: f64, den: u64) -> Rational {
    Rational::from(((x * den as f64).round_ties_even() as i64, den))
}

/// ratround on the den grid, refined by powers of 2 while x * den < minq (>= 6 significant bits)
pub fn ratround_sig(x: f64, den: u64, minq: f64) -> Rational {
    assert!(x > 0.0);
    let mut den = den;
    while x * (den as f64) < minq {
        den <<= 1;
    }
    ratround(x, den)
}

/// x rounded to 24 significant bits, exact dyadic (zplan rr)
pub fn rr(x: f64) -> Rational {
    if x == 0.0 {
        return Rational::new();
    }
    let (_, ex) = frexp(x);
    let k = 24 - ex;
    let m = (x * 2f64.powi(k)).round_ties_even() as i64;
    if k >= 0 {
        Rational::from((rug::Integer::from(m), rug::Integer::from(1) << k as u32))
    } else {
        Rational::from(rug::Integer::from(m) << (-k) as u32)
    }
}

fn frexp(x: f64) -> (f64, i32) {
    if x == 0.0 || !x.is_finite() {
        return (x, 0);
    }
    let bits = x.abs().to_bits();
    let e = ((bits >> 52) & 0x7ff) as i32;
    if e == 0 {
        let (m, e2) = frexp(x * 2f64.powi(64));
        return (m, e2 - 64);
    }
    let ex = e - 1022;
    (x / 2f64.powi(ex), ex)
}

fn qf(x: &Rational) -> f64 {
    x.to_f64()
}

pub struct PathPlan {
    pub patches: Vec<Patch>,
    pub s0: Rational,
    pub s_end: Rational,
    pub rho0: Rational,
    pub rho_inf: Rational,
    pub working_digits: u32,
}

impl PathPlan {
    pub fn n_complex(&self) -> usize {
        self.patches.iter().filter(|p| !p.is_real()).count()
    }
    pub fn sum_orders(&self) -> usize {
        self.patches.iter().map(|p| p.order).sum()
    }
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "patches": self.patches.iter().map(|p| serde_json::json!([p.centre.0.to_string(), p.centre.1.to_string(), p.step.0.to_string(), p.step.1.to_string(), p.order])).collect::<Vec<_>>(),
            "n": self.patches.len(), "n_complex": self.n_complex(), "sum_nord": self.sum_orders(),
            "s0": self.s0.to_string(), "S": self.s_end.to_string(), "rho0": self.rho0.to_string(), "rho_inf": self.rho_inf.to_string(),
            "working_digits": self.working_digits,
        })
    }
}

/// rho0 = 0.95 R0 (>= 6 significant bits on the 2^-10 grid), rho_inf = 1 / ratround(1.05 Rmax, 2^-4)
pub fn radii(r0: f64, rmax: f64) -> (Rational, Rational) {
    let rho0 = ratround_sig(r0 * 0.95, 1 << 10, 64.0);
    let rinf = ratround(rmax * 1.05, 1 << 4);
    (rho0, Rational::from(1) / rinf)
}

/// s0 from the soft stable range (n_st0 orders), S from the infinity stable ranges (n_inf orders)
pub fn endpoints(rho0: &Rational, rmax: f64, w: u32, n_st0: i64, n_inf: i64) -> (Rational, Rational) {
    let w = w as f64;
    let s0 = ratround_sig(qf(rho0) * 10f64.powf(-(w + 5.0) / n_st0 as f64) * 0.9, 1 << 24, 64.0);
    let s = Rational::from(rug::Integer::from_f64((rmax * 1.05 * 10f64.powf((w + 5.0) / n_inf as f64)).ceil()).unwrap());
    (s0, s)
}

/// zplan.py: patches z = -s0 -> -S
pub fn plan(points: &[(f64, f64)], axis: &[f64], s0: &Rational, s_end: &Rational, w: u32, ratio: f64) -> Vec<Patch> {
    let radius = |zr: f64, zi: f64| -> f64 { points.iter().map(|&(a, b)| (zr - a).hypot(zi - b)).fold(f64::INFINITY, f64::min) };
    let nord = |r: f64, h: f64| -> usize { ((w as f64 + 8.0) / (r / h.abs()).log10()).ceil() as usize + 4 };
    let mut patches = Vec::new();
    let step_real = |zq: Rational, target: &Rational, patches: &mut Vec<Patch>| -> Rational {
        let mut zq = zq;
        while zq > *target {
            let zc = qf(&zq);
            let r = radius(zc, 0.0);
            let h = r / ratio;
            let hq = if zc - h * 1.5 <= qf(target) { Rational::from(target - &zq) } else { -rr(h) };
            patches.push(Patch { centre: (zq.clone(), Rational::new()), step: (hq.clone(), Rational::new()), order: nord(r, qf(&hq)) });
            zq += hq;
        }
        zq
    };
    let mut zq = -s0.clone();
    let mend = -s_end.clone();
    for &xp in axis {
        if !(qf(&mend) < xp && xp < qf(&zq)) {
            continue;
        }
        let gap = points
            .iter()
            .filter(|&&(a, b)| (a - xp).hypot(b) > 1e-9 * (1.0 + xp.abs()))
            .map(|&(a, b)| (a - xp).hypot(b))
            .fold(xp.abs(), f64::min);
        let rd = 0.4 * gap;
        let start = rr(xp + rd);
        let end = rr(xp - rd);
        if start < zq {
            zq = step_real(zq, &start, &mut patches);
        }
        let nv = 12;
        let mut verts: Vec<(Rational, Rational)> = (1..nv)
            .map(|k| {
                let th = std::f64::consts::PI * k as f64 / nv as f64;
                (rr(xp + rd * th.cos()), rr(rd * th.sin()))
            })
            .collect();
        verts.push((end, Rational::new()));
        let (mut cr, mut ci) = (zq.clone(), Rational::new());
        for (vr, vi) in verts {
            while (&cr, &ci) != (&vr, &vi) {
                let (fr, fi) = (qf(&cr), qf(&ci));
                let r = radius(fr, fi);
                let (dr, di) = (qf(&Rational::from(&vr - &cr)), qf(&Rational::from(&vi - &ci)));
                let dn = dr.hypot(di);
                let hq = if dn <= 1.5 * r / ratio {
                    (Rational::from(&vr - &cr), Rational::from(&vi - &ci))
                } else {
                    (rr(dr / dn * r / ratio), rr(di / dn * r / ratio))
                };
                let hn = qf(&hq.0).hypot(qf(&hq.1));
                patches.push(Patch { centre: (cr.clone(), ci.clone()), step: hq.clone(), order: nord(r, hn) });
                cr += &hq.0;
                ci += &hq.1;
            }
        }
        assert!(ci == 0);
        zq = cr;
    }
    let zq = step_real(zq, &mend, &mut patches);
    assert!(zq == mend);
    let _ = poly::trim;
    patches
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rr_is_24_bits() {
        let x = rr(0.1);
        assert_eq!(x.denom().significant_bits() - 1, 27);
        assert!((x.to_f64() - 0.1).abs() < 1e-8);
        assert_eq!(rr(-3.0), Rational::from(-3));
    }
    #[test]
    fn small_path_detours_and_ends() {
        // singular points 0, -1, -3 (+ a complex pair)
        let pts = vec![(0.0, 0.0), (-1.0, 0.0), (-3.0, 0.0), (-2.0, 2.0), (-2.0, -2.0)];
        let axis = vec![-1.0, -3.0];
        let s0 = Rational::from((1, 1 << 20));
        let s = Rational::from(10);
        let p = plan(&pts, &axis, &s0, &s, 30, 3.0);
        let mut z = (-s0.clone(), Rational::new());
        for q in &p {
            assert_eq!((&z.0, &z.1), (&q.centre.0, &q.centre.1));
            z.0 += &q.step.0;
            z.1 += &q.step.1;
            // every patch stays inside its convergence disk with margin
            let r = pts.iter().map(|&(a, b)| (q.centre.0.to_f64() - a).hypot(q.centre.1.to_f64() - b)).fold(f64::INFINITY, f64::min);
            assert!(q.step.0.to_f64().hypot(q.step.1.to_f64()) <= 0.51 * r, "patch too long");
        }
        assert_eq!(z, (Rational::from(-10), Rational::new()));
        assert!(p.iter().filter(|q| !q.is_real()).count() >= 2 * 11);
    }
}
