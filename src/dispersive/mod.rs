//! DISPERSIVE MODE (W3): the value at z = 1 of a master integral from the discontinuity of its eta-DE across the
//! tachyonic cut z = 1 + eta < 0 (one-line dispersive gluing; oracle = the campaign's in-house tools soft.py /
//! soft_gen.py, prep5.py, zplan.py, disp3y-pt, post.py). See W3A-DESIGN.md for the layout.
//!
//! The core modules depend only on rug + serde_json (fixed-precision MPFR, exact rug rationals); `harness`
//! connects them to RustFlow's RustRed harness reader.
pub mod path;
pub mod poly;
pub mod singular;
pub mod system;
pub mod field;
pub mod frobenius;
pub mod laurent;
pub mod boundary;
pub mod cli;
pub mod dispersion;
#[cfg(all(feature = "native", feature = "automatic"))]
pub mod harness;
pub mod run;
pub mod transport;
#[cfg(test)]
mod tests;
pub mod softplan;
