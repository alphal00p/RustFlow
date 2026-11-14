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
    /// exact zero pattern and first-nonzero pivoting (ModP)
    const EXACT: bool;
    /// carries an MPFR value that `snap` cleans (Float)
    const SNAP: bool;
    /// exact serialisation for the in-family checkpoints (W3b, ckpt.rs); false / None: this type is not checkpointed
    fn ckpt_put(&self, _out: &mut Vec<u8>) -> bool {
        false
    }
    fn ckpt_get(_b: &[u8], _p: &mut usize, _prec: u32) -> Option<Self> {
        None
    }
}

fn take<'a>(b: &'a [u8], p: &mut usize, n: usize) -> Option<&'a [u8]> {
    let s = b.get(*p..*p + n)?;
    *p += n;
    Some(s)
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
        if let Some(ex) = self.get_exp()
            && emax != i32::MIN
            && ex < emax.saturating_sub(keep)
        {
            self.assign(rug::float::Special::Zero);
        }
    }
    fn log10_abs(&self) -> f64 {
        if Float::is_zero(self) {
            -1e9
        } else {
            self.clone().abs().log10().to_f64()
        }
    }
    const EXACT: bool = false;
    const SNAP: bool = true;
    /// tag (0 +0, 4 -0, 1 positive, 2 negative), precision, then exponent and mantissa digits: x = m 2^e exactly
    fn ckpt_put(&self, out: &mut Vec<u8>) -> bool {
        if !self.is_finite() {
            return false;
        }
        if Float::is_zero(self) {
            out.push(if self.is_sign_negative() { 4 } else { 0 });
            out.extend_from_slice(&self.prec().to_le_bytes());
            return true;
        }
        let (m, e) = self.to_integer_exp().expect("finite");
        out.push(if m < 0 { 2 } else { 1 });
        out.extend_from_slice(&self.prec().to_le_bytes());
        out.extend_from_slice(&e.to_le_bytes());
        let d = m.as_abs().to_digits::<u8>(rug::integer::Order::Lsf);
        out.extend_from_slice(&(d.len() as u32).to_le_bytes());
        out.extend_from_slice(&d);
        true
    }
    fn ckpt_get(b: &[u8], p: &mut usize, _prec: u32) -> Option<Self> {
        let tag = take(b, p, 1)?[0];
        let pr = u32::from_le_bytes(take(b, p, 4)?.try_into().ok()?);
        match tag {
            0 => Some(Float::new(pr)),
            4 => Some(-Float::new(pr)),
            1 | 2 => {
                let e = i32::from_le_bytes(take(b, p, 4)?.try_into().ok()?);
                let n = u32::from_le_bytes(take(b, p, 4)?.try_into().ok()?) as usize;
                let m = rug::Integer::from_digits(take(b, p, n)?, rug::integer::Order::Lsf);
                let mut x = Float::with_val(pr, m);
                x <<= e;
                Some(if tag == 2 { -x } else { x })
            }
            _ => None,
        }
    }
}

/// power series n/d to k terms (d[0] != 0)
pub fn series_div<S: Scalar>(n: &[S], d: &[S], k: usize, prec: u32) -> Vec<S> {
    let mut out: Vec<S> = Vec::with_capacity(k);
    let mut tmp = S::zero(prec);
    for i in 0..k {
        let mut acc = if i < n.len() {
            n[i].clone()
        } else {
            S::zero(prec)
        };
        for m in 1..=i.min(d.len() - 1) {
            acc.sub_mul(&d[m], &out[i - m], &mut tmp);
        }
        out.push(acc.div_s(&d[0], prec));
    }
    out
}
