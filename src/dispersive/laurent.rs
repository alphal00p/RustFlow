//! Laurent matrices of the e-block at z = 0 and z = infinity (oracle: disp `laurent`, local_solutions.laurent_matrices).
//!   point 0:   x = z,     theta_x D = B(x) D,  B = z A(z)
//!   point inf: x = w = 1/z, theta_x D = C(x) D, C = -(1/w) A(1/w)
//! Exact rational shifts (eta = z - 1) and reversals, then the series of every entry in the scalar type S
//! (mod p or MPFR) to kmax orders. Scaling x = rho x' multiplies B_k by rho^k (applied per family).
use super::field::{Scalar, series_div};
use super::poly::{self, QPoly};
use super::system::EBlockSystem;
use rug::Rational;
use std::collections::BTreeMap;

pub type Lmat<S> = BTreeMap<i64, Vec<(usize, usize, S)>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Point {
    Zero,
    Infinity,
}

/// exact (numerator, denominator, shift) of one entry's local expansion: entry = x^shift n(x)/d(x), n(0), d(0) != 0
fn local_entry(num: &poly::ZPoly, den: &poly::ZPoly, point: Point) -> Option<(QPoly, QPoly, i64)> {
    let mut n = poly::shift_minus_one(num);
    let mut d = poly::shift_minus_one(den);
    poly::trim(&mut n);
    poly::trim(&mut d);
    let vn = poly::valuation(&n)?;
    match point {
        Point::Zero => {
            let vd = poly::valuation(&d).unwrap();
            let n2 = n[vn..].to_vec();
            let d2 = d[vd..].to_vec();
            Some((n2, d2, 1 + vn as i64 - vd as i64))
        }
        Point::Infinity => {
            // c(1/w - 1) = g(1/w), g = c(. - 1): w^N c(1/w - 1) = reversed g
            let nn = n.len() as i64 - 1;
            let nd = d.len() as i64 - 1;
            let rn: QPoly = n.iter().rev().cloned().collect();
            let rd: QPoly = d.iter().rev().cloned().collect();
            let vn = poly::valuation(&rn).unwrap();
            let vd = poly::valuation(&rd).unwrap();
            let n2: QPoly = rn[vn..].iter().map(|x| Rational::from(-x)).collect();
            let d2 = rd[vd..].to_vec();
            Some((n2, d2, (nd - nn) + (vn as i64 - vd as i64) - 1))
        }
    }
}

/// unscaled Laurent matrices B_k, k <= kmax, in S
pub fn laurent<S: Scalar>(
    sys: &EBlockSystem,
    point: Point,
    kmax: i64,
    prec: u32,
    nthreads: usize,
) -> Lmat<S> {
    let nt = nthreads.max(1);
    let chunk = sys.terms.len().div_ceil(nt).max(1);
    let parts: Vec<Vec<(i64, usize, usize, S)>> = std::thread::scope(|sc| {
        let hs: Vec<_> = sys
            .terms
            .chunks(chunk)
            .map(|tch| {
                sc.spawn(move || {
                    let mut out = Vec::new();
                    for t in tch {
                        let Some((n, d, shift)) = local_entry(&t.num, &t.den, point) else {
                            continue;
                        };
                        let kk = (kmax - shift).max(1) as usize;
                        let nf: Vec<S> = n.iter().map(|x| S::from_rational(x, prec)).collect();
                        let df: Vec<S> = d.iter().map(|x| S::from_rational(x, prec)).collect();
                        for (k, v) in series_div(&nf, &df, kk, prec).into_iter().enumerate() {
                            if !v.is_zero() {
                                out.push((k as i64 + shift, t.i, t.j, v));
                            }
                        }
                    }
                    out
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut acc: BTreeMap<(i64, usize, usize), S> = BTreeMap::new();
    for part in parts {
        for (k, i, j, v) in part {
            if k > kmax {
                continue;
            }
            match acc.get_mut(&(k, i, j)) {
                Some(x) => x.add_assign_s(&v),
                None => {
                    acc.insert((k, i, j), v);
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

/// B_k -> B_k rho^k
pub fn scaled<S: Scalar>(b: &Lmat<S>, rho: &Rational, prec: u32) -> Lmat<S> {
    b.iter()
        .map(|(&k, v)| {
            let rk = S::from_rational(&pow_q(rho, k), prec);
            (
                k,
                v.iter()
                    .map(|(i, j, x)| (*i, *j, x.mul_s(&rk, prec)))
                    .collect(),
            )
        })
        .collect()
}

pub fn pow_q(x: &Rational, k: i64) -> Rational {
    let mut r = Rational::from(1);
    let b = if k >= 0 {
        x.clone()
    } else {
        Rational::from(1) / x
    };
    for _ in 0..k.unsigned_abs() {
        r *= &b;
    }
    r
}
