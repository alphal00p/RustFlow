//! Number types of the patch kernel.
//!
//! [`KernelScalar`] is the small set of in-place operations the recurrence
//! needs. Fixed precision (default): `rug::Float` / `rug::Complex`, every
//! operation correctly rounded to the working precision, no allocation in the
//! inner loop. Tracked (option): Symbolica's `Float` / `Complex<Float>` with
//! their own precision-tracking operators; the precision left in a result is
//! the tracking accuracy indicator.
//!
//! The exact state between patches is always stored as `rug::Complex`
//! (component precisions preserved), so checkpoints and resumes are exact in
//! both modes.
use rug::{Assign, Complex, Float, Rational};
use serde::{Deserialize, Serialize};
use symbolica::domains::float::{Complex as SymComplex, Float as SymFloat};

/// Arithmetic of the transport.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArithmeticMode {
    /// MPFR at one working precision (correctly rounded operations).
    Fixed,
    /// Symbolica's precision-tracking floats (slower; diagnostic option).
    Tracked,
}

impl std::str::FromStr for ArithmeticMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "fixed" => Ok(Self::Fixed),
            "tracked" => Ok(Self::Tracked),
            _ => Err(format!("unknown arithmetic mode {s:?} (fixed|tracked)")),
        }
    }
}

/// In-place arithmetic used by the patch recurrence.
pub trait KernelScalar: Clone + Send + Sync {
    /// Complex arithmetic (false: real values, imaginary parts never formed).
    const COMPLEX: bool;
    fn zero(bits: u32) -> Self;
    /// `a / b` rounded to `bits`, from guarded (higher-precision) values.
    /// Real scalars use the real parts only.
    fn quotient(a: &Complex, b: &Complex, bits: u32) -> Self;
    /// An exact rational complex number rounded to `bits`.
    fn from_rational(re: &Rational, im: &Rational, bits: u32) -> Self;
    /// From the exact stored state (real scalars take the real part).
    fn from_state(v: &Complex, bits: u32) -> Self;
    /// To the stored state, exactly.
    fn to_state(&self) -> Complex;
    fn set_zero(&mut self);
    /// `self += x * y` (tmp is scratch space).
    fn add_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self);
    fn add_ref(&mut self, o: &Self);
    fn sub_ref(&mut self, o: &Self);
    fn mul_ref(&mut self, o: &Self);
    fn div_ref(&mut self, o: &Self);
    fn div_u32(&mut self, u: u32);
}

impl KernelScalar for Float {
    const COMPLEX: bool = false;
    fn zero(bits: u32) -> Self {
        Float::new(bits)
    }
    fn quotient(a: &Complex, b: &Complex, bits: u32) -> Self {
        Float::with_val(bits, a.real() / b.real())
    }
    fn from_rational(re: &Rational, _im: &Rational, bits: u32) -> Self {
        Float::with_val(bits, re)
    }
    fn from_state(v: &Complex, bits: u32) -> Self {
        let mut x = Float::new(bits);
        x.assign(v.real());
        x
    }
    fn to_state(&self) -> Complex {
        Complex::with_val((self.prec(), self.prec()), (self, 0))
    }
    #[inline]
    fn set_zero(&mut self) {
        self.assign(rug::float::Special::Zero);
    }
    #[inline]
    fn add_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self) {
        tmp.assign(x * y);
        *self += &*tmp;
    }
    #[inline]
    fn add_ref(&mut self, o: &Self) {
        *self += o;
    }
    #[inline]
    fn sub_ref(&mut self, o: &Self) {
        *self -= o;
    }
    #[inline]
    fn mul_ref(&mut self, o: &Self) {
        *self *= o;
    }
    #[inline]
    fn div_ref(&mut self, o: &Self) {
        *self /= o;
    }
    #[inline]
    fn div_u32(&mut self, u: u32) {
        *self /= u;
    }
}

impl KernelScalar for Complex {
    const COMPLEX: bool = true;
    fn zero(bits: u32) -> Self {
        Complex::new(bits)
    }
    fn quotient(a: &Complex, b: &Complex, bits: u32) -> Self {
        Complex::with_val(bits, a / b)
    }
    fn from_rational(re: &Rational, im: &Rational, bits: u32) -> Self {
        Complex::with_val(bits, (re, im))
    }
    fn from_state(v: &Complex, bits: u32) -> Self {
        let mut x = Complex::new(bits);
        x.assign(v);
        x
    }
    fn to_state(&self) -> Complex {
        self.clone()
    }
    #[inline]
    fn set_zero(&mut self) {
        self.mut_real().assign(rug::float::Special::Zero);
        self.mut_imag().assign(rug::float::Special::Zero);
    }
    #[inline]
    fn add_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self) {
        tmp.assign(x * y);
        *self += &*tmp;
    }
    #[inline]
    fn add_ref(&mut self, o: &Self) {
        *self += o;
    }
    #[inline]
    fn sub_ref(&mut self, o: &Self) {
        *self -= o;
    }
    #[inline]
    fn mul_ref(&mut self, o: &Self) {
        *self *= o;
    }
    #[inline]
    fn div_ref(&mut self, o: &Self) {
        *self /= o;
    }
    #[inline]
    fn div_u32(&mut self, u: u32) {
        *self /= u;
    }
}

/// Real precision-tracking scalar (Symbolica `Float`).
#[derive(Clone)]
pub struct TrackedReal(pub SymFloat);
/// Complex precision-tracking scalar (Symbolica `Complex<Float>`).
#[derive(Clone)]
pub struct TrackedComplex(pub SymComplex<SymFloat>);

fn sym(x: &Float) -> SymFloat {
    SymFloat::from(x.clone())
}

impl KernelScalar for TrackedReal {
    const COMPLEX: bool = false;
    fn zero(bits: u32) -> Self {
        Self(SymFloat::new(bits))
    }
    fn quotient(a: &Complex, b: &Complex, bits: u32) -> Self {
        Self(sym(&<Float as KernelScalar>::quotient(a, b, bits)))
    }
    fn from_rational(re: &Rational, _im: &Rational, bits: u32) -> Self {
        Self(sym(&Float::with_val(bits, re)))
    }
    fn from_state(v: &Complex, _bits: u32) -> Self {
        Self(sym(v.real()))
    }
    fn to_state(&self) -> Complex {
        let x = self.0.to_raw();
        let p = x.prec();
        Complex::with_val((p, p), (x, 0))
    }
    fn set_zero(&mut self) {
        let p = self.0.prec();
        self.0 = SymFloat::new(p);
    }
    fn add_prod(&mut self, x: &Self, y: &Self, _tmp: &mut Self) {
        self.0 += &(x.0.clone() * &y.0);
    }
    fn add_ref(&mut self, o: &Self) {
        self.0 += &o.0;
    }
    fn sub_ref(&mut self, o: &Self) {
        self.0 -= &o.0;
    }
    fn mul_ref(&mut self, o: &Self) {
        self.0 *= &o.0;
    }
    fn div_ref(&mut self, o: &Self) {
        self.0 = self.0.clone() / &o.0;
    }
    fn div_u32(&mut self, u: u32) {
        // An exact integer does not reduce the tracked precision.
        let d = SymFloat::with_val(self.0.prec().max(32), u);
        self.0 = self.0.clone() / &d;
    }
}

impl KernelScalar for TrackedComplex {
    const COMPLEX: bool = true;
    fn zero(bits: u32) -> Self {
        Self(SymComplex::new(SymFloat::new(bits), SymFloat::new(bits)))
    }
    fn quotient(a: &Complex, b: &Complex, bits: u32) -> Self {
        let q = <Complex as KernelScalar>::quotient(a, b, bits);
        Self::from_state(&q, bits)
    }
    fn from_rational(re: &Rational, im: &Rational, bits: u32) -> Self {
        Self::from_state(&Complex::with_val(bits, (re, im)), bits)
    }
    fn from_state(v: &Complex, _bits: u32) -> Self {
        Self(SymComplex::new(sym(v.real()), sym(v.imag())))
    }
    fn to_state(&self) -> Complex {
        let re = self.0.re.to_raw();
        let im = self.0.im.to_raw();
        Complex::with_val((re.prec(), im.prec()), (re, im))
    }
    fn set_zero(&mut self) {
        let (a, b) = (self.0.re.prec(), self.0.im.prec());
        self.0 = SymComplex::new(SymFloat::new(a), SymFloat::new(b));
    }
    fn add_prod(&mut self, x: &Self, y: &Self, _tmp: &mut Self) {
        self.0 += &(x.0.clone() * &y.0);
    }
    fn add_ref(&mut self, o: &Self) {
        self.0 += &o.0;
    }
    fn sub_ref(&mut self, o: &Self) {
        self.0 = self.0.clone() - &o.0;
    }
    fn mul_ref(&mut self, o: &Self) {
        self.0 *= &o.0;
    }
    fn div_ref(&mut self, o: &Self) {
        self.0 = self.0.clone() / &o.0;
    }
    fn div_u32(&mut self, u: u32) {
        let d = SymFloat::with_val(self.0.re.prec().max(32), u);
        self.0 = SymComplex::new(self.0.re.clone() / &d, self.0.im.clone() / &d);
    }
}
