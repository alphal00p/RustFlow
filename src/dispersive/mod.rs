//! DISPERSIVE MODE (W3): the value at z = 1 of a master integral from the discontinuity of its eta-DE across the
//! tachyonic cut z = 1 + eta < 0 (one-line dispersive gluing; oracle = the campaign's in-house tools soft.py /
//! soft_gen.py, prep5.py, zplan.py, disp3y-pt, post.py). See W3A-DESIGN.md for the layout.
//!
//! The core modules depend only on rug + serde_json (fixed-precision MPFR, exact rug rationals); `harness`
//! connects them to RustFlow's RustRed harness reader.
// Index loops in the numeric kernels are deliberate: the order of operations is part of the bit-identical (thread-count independent)
// results. Negated comparisons (`!(x <= bound)`) are deliberate: they treat NaN as a failure.
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]
pub mod boundary;
pub mod ckpt;
pub mod cli;
pub mod dispersion;
pub mod field;
pub mod frobenius;
#[cfg(all(feature = "native", feature = "automatic"))]
pub mod harness;
pub mod kernel;
pub mod lane;
pub mod laurent;
pub mod lower;
pub mod path;
pub mod poly;
pub mod run;
pub mod singular;
pub mod softplan;
pub mod system;
#[cfg(test)]
mod tests;
pub mod transport;

/// Decimal digits per bit as the reference tool writes it. Kept as this literal (not `std::f64::consts::LOG10_2`) so that
/// precision thresholds and printed digit counts stay identical to the reference tool and to existing checkpoints.
#[allow(clippy::approx_constant)]
pub const DIGITS_PER_BIT: f64 = 0.30103;

/// Bits per decimal digit as the reference tool writes it (see [`DIGITS_PER_BIT`]): working precisions `bits(W)` depend on it.
#[allow(clippy::approx_constant)]
pub const BITS_PER_DIGIT: f64 = 3.3219281;
