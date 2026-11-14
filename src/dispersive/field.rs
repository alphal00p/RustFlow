//! Scalars of the linear-algebra Frobenius: exact residues mod p (structure: ranks, pivot coordinates) and
//! fixed-precision MPFR floats (values). One generic implementation runs on both.
use rug::{Assign, Float, Rational};
use std::cmp::Ordering;

pub const PRIME: u64 = 1048573; // 2^20 - 3 (oracle: local_solutions.py P_DEFAULT)

pub trait Scalar: Clone + Send + Sync + std::fmt::Debug {
    /// working precision in bits (ignored by ModP)
    fn zero(prec: u32) -> Self;
    fn one(prec: u32) -> Self;
    fn from_rational(x: &Rational, prec: u32) -> Self;
    fn is_zero(&self) -> bool;
    fn set_zero(&mut self);
    fn add_assign_s(&mut self, o: &Self);
    /// self -= a * b; returns the binary exponent of a * b (i32::MIN for exact types)
    fn sub_mul(&mut self, a: &Self, b: &Self, tmp: &mut Self) -> i32;
    fn mul_s(&self, o: &Self, prec: u32) -> Self;
    fn div_s(&self, o: &Self, prec: u32) -> Self;
    fn neg_s(&self) -> Self;
    /// |self| compared with |o| (ModP: nonzero > zero)
    fn cmp_mag(&self, o: &Self) -> Ordering;
    fn exponent(&self) -> i32;
    /// snap to zero when the value lies below 2^(emax - keep) (cancellation noise); no-op when exact
    fn snap(&mut self, emax: i32, keep: i32);
    fn log10_abs(&self) -> f64;
    /// exact zero pattern and first-nonzero pivoting (ModP, Shadowed)
    const EXACT: bool;
    /// carries an MPFR value that `snap` cleans (Float, Shadowed)
    const SNAP: bool;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModP(pub u32);

impl ModP {
    pub fn inv(self) -> ModP {
        ModP(super::poly::pow_mod(self.0 as u64, PRIME - 2, PRIME) as u32)
    }
}

impl Scalar for ModP {
    fn zero(_: u32) -> Self {
        ModP(0)
    }
    fn one(_: u32) -> Self {
        ModP(1)
    }
    fn from_rational(x: &Rational, _: u32) -> Self {
        ModP(super::poly::rat_mod_p(x, PRIME) as u32)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    fn set_zero(&mut self) {
        self.0 = 0;
    }
    #[inline]
    fn add_assign_s(&mut self, o: &Self) {
        self.0 = ((self.0 as u64 + o.0 as u64) % PRIME) as u32;
    }
    #[inline]
    fn sub_mul(&mut self, a: &Self, b: &Self, _tmp: &mut Self) -> i32 {
        let p = (a.0 as u64 * b.0 as u64) % PRIME;
        self.0 = ((self.0 as u64 + PRIME - p) % PRIME) as u32;
        i32::MIN
    }
    fn mul_s(&self, o: &Self, _: u32) -> Self {
        ModP(((self.0 as u64 * o.0 as u64) % PRIME) as u32)
    }
    fn div_s(&self, o: &Self, _: u32) -> Self {
        assert!(o.0 != 0, "division by zero mod p");
        self.mul_s(&o.inv(), 0)
    }
    fn neg_s(&self) -> Self {
        ModP(((PRIME - self.0 as u64) % PRIME) as u32)
    }
    fn cmp_mag(&self, o: &Self) -> Ordering {
        (self.0 != 0).cmp(&(o.0 != 0))
    }
    fn exponent(&self) -> i32 {
        i32::MIN
    }
    fn snap(&mut self, _: i32, _: i32) {}
    fn log10_abs(&self) -> f64 {
        if self.0 == 0 { -1e9 } else { 0.0 }
    }
    const EXACT: bool = true;
    const SNAP: bool = false;
}

impl Scalar for Float {
    fn zero(prec: u32) -> Self {
        Float::new(prec)
    }
    fn one(prec: u32) -> Self {
        Float::with_val(prec, 1)
    }
    fn from_rational(x: &Rational, prec: u32) -> Self {
        Float::with_val(prec, x)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        Float::is_zero(self)
    }
    fn set_zero(&mut self) {
        self.assign(rug::float::Special::Zero);
    }
    #[inline]
    fn add_assign_s(&mut self, o: &Self) {
        *self += o;
    }
    #[inline]
    fn sub_mul(&mut self, a: &Self, b: &Self, tmp: &mut Self) -> i32 {
        tmp.assign(a * b);
        *self -= &*tmp;
        tmp.get_exp().unwrap_or(i32::MIN)
    }
    fn mul_s(&self, o: &Self, prec: u32) -> Self {
        Float::with_val(prec, self * o)
    }
    fn div_s(&self, o: &Self, prec: u32) -> Self {
        Float::with_val(prec, self / o)
    }
    fn neg_s(&self) -> Self {
        Float::with_val(self.prec(), -self)
    }
    fn cmp_mag(&self, o: &Self) -> Ordering {
        self.cmp_abs(o).unwrap_or(Ordering::Equal)
    }
    fn exponent(&self) -> i32 {
        self.get_exp().unwrap_or(i32::MIN)
    }
    #[inline]
    fn snap(&mut self, emax: i32, keep: i32) {
        if let Some(ex) = self.get_exp() {
            if emax != i32::MIN && ex < emax.saturating_sub(keep) {
                self.assign(rug::float::Special::Zero);
            }
        }
    }
    fn log10_abs(&self) -> f64 {
        if Float::is_zero(self) { -1e9 } else { self.clone().abs().log10().to_f64() }
    }
    const EXACT: bool = false;
    const SNAP: bool = true;
}

/// MPFR value with an exact shadow residue mod 2^61 - 1: the shadow decides which entries are exactly zero (no rounding noise and no
/// spurious fill survive) and the pivots are the exact run's (first nonzero); the float carries the value. A nonzero rational has a zero
/// residue with probability ~ 2^-61 per entry.
#[derive(Clone, Debug)]
pub struct Shadowed {
    pub m: u64,
    pub f: Float,
}

pub const P61: u64 = (1u64 << 61) - 1;

#[inline]
fn mul61(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P61 as u128) as u64
}

fn pow61(mut b: u64, mut e: u64) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = mul61(r, b);
        }
        b = mul61(b, b);
        e >>= 1;
    }
    r
}

pub fn rat_mod61(x: &Rational) -> u64 {
    let p = rug::Integer::from(P61);
    let red = |v: &rug::Integer| -> u64 {
        let mut r = rug::Integer::from(v % &p);
        if r < 0 {
            r += &p;
        }
        r.to_u64().unwrap()
    };
    let n = red(x.numer());
    let d = red(x.denom());
    assert!(d != 0, "denominator divisible by 2^61 - 1");
    mul61(n, pow61(d, P61 - 2))
}

impl Scalar for Shadowed {
    fn zero(prec: u32) -> Self {
        Shadowed { m: 0, f: Float::new(prec) }
    }
    fn one(prec: u32) -> Self {
        Shadowed { m: 1, f: Float::with_val(prec, 1) }
    }
    fn from_rational(x: &Rational, prec: u32) -> Self {
        Shadowed { m: rat_mod61(x), f: Float::with_val(prec, x) }
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.m == 0
    }
    fn set_zero(&mut self) {
        self.m = 0;
        self.f.assign(rug::float::Special::Zero);
    }
    #[inline]
    fn add_assign_s(&mut self, o: &Self) {
        self.m = (self.m + o.m) % P61;
        self.f += &o.f;
    }
    #[inline]
    fn sub_mul(&mut self, a: &Self, b: &Self, tmp: &mut Self) -> i32 {
        self.m = (self.m + P61 - mul61(a.m, b.m)) % P61;
        tmp.f.assign(&a.f * &b.f);
        self.f -= &tmp.f;
        tmp.f.get_exp().unwrap_or(i32::MIN)
    }
    fn mul_s(&self, o: &Self, prec: u32) -> Self {
        Shadowed { m: mul61(self.m, o.m), f: Float::with_val(prec, &self.f * &o.f) }
    }
    fn div_s(&self, o: &Self, prec: u32) -> Self {
        assert!(o.m != 0, "division by an exact zero");
        Shadowed { m: mul61(self.m, pow61(o.m, P61 - 2)), f: Float::with_val(prec, &self.f / &o.f) }
    }
    fn neg_s(&self) -> Self {
        Shadowed { m: (P61 - self.m) % P61, f: Float::with_val(self.f.prec(), -&self.f) }
    }
    fn cmp_mag(&self, o: &Self) -> Ordering {
        match (self.m == 0, o.m == 0) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => self.f.cmp_abs(&o.f).unwrap_or(Ordering::Equal),
        }
    }
    fn exponent(&self) -> i32 {
        self.f.get_exp().unwrap_or(i32::MIN)
    }
    #[inline]
    fn snap(&mut self, _emax: i32, _keep: i32) {
        if self.m == 0 {
            self.f.assign(rug::float::Special::Zero);
        }
    }
    fn log10_abs(&self) -> f64 {
        if self.m == 0 || self.f.is_zero() { -1e9 } else { self.f.clone().abs().log10().to_f64() }
    }
    const EXACT: bool = true;
    const SNAP: bool = true;
}

/// power series n/d to k terms (d[0] != 0)
pub fn series_div<S: Scalar>(n: &[S], d: &[S], k: usize, prec: u32) -> Vec<S> {
    let mut out: Vec<S> = Vec::with_capacity(k);
    let mut tmp = S::zero(prec);
    for i in 0..k {
        let mut acc = if i < n.len() { n[i].clone() } else { S::zero(prec) };
        for m in 1..=i.min(d.len() - 1) {
            acc.sub_mul(&d[m], &out[i - m], &mut tmp);
        }
        out.push(acc.div_s(&d[0], prec));
    }
    out
}
