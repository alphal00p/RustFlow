//! Laurent matrices of dY/dz = A(z) Y at an arbitrary point of the z-plane, in the native variable of the system
//! (`EBlockSystem` coefficients are polynomials in x with z = offset + x; offset 0 for the propagator families of PFLOW,
//! 1 for the vacuum classes of the dispersive mode):
//!   finite c:  x_loc = z - c,  theta Y = B(x_loc) Y,  B = x_loc A(c + x_loc)
//!   infinity:  w = 1/z,        theta Y = C(w) Y,      C = -(1/w) A(1/w)
//! Exact rational shift / reversal, then the series of each entry in the scalar type S (mod p or MPFR).
use crate::dispersive::field::{Scalar, series_div};
use crate::dispersive::laurent::Lmat;
use crate::dispersive::poly::{self, QPoly};
use crate::dispersive::system::EBlockSystem;
use rug::Rational;
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub enum Pt {
    Finite(Rational),
    Infinity,
}

/// p(x + c) exactly
pub fn shift(p: &[Rational], c: &Rational) -> QPoly {
    let mut a: QPoly = p.to_vec();
    let n = a.len();
    if n < 2 || *c == 0 {
        return a;
    }
    for i in 0..n - 1 {
        for j in (i..n - 1).rev() {
            let t = Rational::from(c * &a[j + 1]);
            a[j] += t;
        }
    }
    a
}

fn local_entry(num: &QPoly, den: &QPoly, offset: &Rational, pt: &Pt) -> Option<(QPoly, QPoly, i64)> {
    // native variable: z = offset + x  ->  polynomials in z: p_z(z) = p(z - offset)
    let mo = Rational::from(-offset);
    let mut nz = shift(num, &mo);
    let mut dz = shift(den, &mo);
    poly::trim(&mut nz);
    poly::trim(&mut dz);
    poly::valuation(&nz)?;
    match pt {
        Pt::Finite(c) => {
            let mut n = shift(&nz, c);
            let mut d = shift(&dz, c);
            poly::trim(&mut n);
            poly::trim(&mut d);
            let vn = poly::valuation(&n)?;
            let vd = poly::valuation(&d).unwrap();
            Some((n[vn..].to_vec(), d[vd..].to_vec(), 1 + vn as i64 - vd as i64))
        }
        Pt::Infinity => {
            let nn = nz.len() as i64 - 1;
            let nd = dz.len() as i64 - 1;
            let rn: QPoly = nz.iter().rev().cloned().collect();
            let rd: QPoly = dz.iter().rev().cloned().collect();
            let vn = poly::valuation(&rn).unwrap();
            let vd = poly::valuation(&rd).unwrap();
            let n2: QPoly = rn[vn..].iter().map(|x| Rational::from(-x)).collect();
            Some((n2, rd[vd..].to_vec(), (nd - nn) + (vn as i64 - vd as i64) - 1))
        }
    }
}

/// Laurent matrices B_k, k <= kmax, of the system at `pt` (offset: z = offset + x of the stored polynomials)
pub fn lmat<S: Scalar>(sys: &EBlockSystem, offset: &Rational, pt: &Pt, kmax: i64, prec: u32) -> Lmat<S> {
    let mut acc: BTreeMap<(i64, usize, usize), S> = BTreeMap::new();
    for t in &sys.terms {
        let Some((n, d, shift_k)) = local_entry(&t.num.to_q(), &t.den.to_q(), offset, pt) else { continue };
        let kk = (kmax - shift_k).max(1) as usize;
        let nf: Vec<S> = n.iter().map(|x| S::from_rational(x, prec)).collect();
        let df: Vec<S> = d.iter().map(|x| S::from_rational(x, prec)).collect();
        for (k, v) in series_div(&nf, &df, kk, prec).into_iter().enumerate() {
            let kk = k as i64 + shift_k;
            if v.is_zero() || kk > kmax {
                continue;
            }
            match acc.get_mut(&(kk, t.i, t.j)) {
                Some(x) => x.add_assign_s(&v),
                None => {
                    acc.insert((kk, t.i, t.j), v);
                }
            }
        }
    }
    let mut out: Lmat<S> = BTreeMap::new();
    for ((k, i, j), v) in acc {
        if !v.is_zero() {
            out.entry(k).or_default().push((i, j, v));
        }
    }
    out
}

/// singular points of the system in z (numerical, for path planning): roots of the distinct denominators
pub fn denominators_z(sys: &EBlockSystem, offset: &Rational) -> Vec<QPoly> {
    let mo = Rational::from(-offset);
    let mut out: Vec<QPoly> = Vec::new();
    for t in &sys.terms {
        let mut d = shift(&t.den.to_q(), &mo);
        poly::trim(&mut d);
        if d.len() > 1 && !out.contains(&d) {
            out.push(d);
        }
    }
    out
}
