//! PFLOW: massless propagator integrals (p-integrals) at a fixed rational d0 from a ONE-MASS flow, dispersive route (B)
//! (PFLOW-S1.md / PFLOW-S1B.md). The engine reuses the dispersive mode's exact machinery (`dispersive::{system, frobenius,
//! field, dispersion}`: e-block of a RustRed harness DE, Frobenius families with exact rank decisions over a prime field and
//! MPFR values, relative least squares) and the fixed-precision patch kernel (`patch_kernel`, Cauchy weight (a - z)^-m).
//! Configuration-driven: the family graph (lines, external vertices) and the massive line come from the input; nothing is
//! topology specific.
//!  * `local`  : Laurent matrices at any point of the z-plane (finite or infinity) in the system's own variable
//!  * `graph`  : region enumeration and leading coefficients (soft sets at 0, hard subgraphs at infinity), massless evaluator
//!  * `engine` : families per exponent class at 0 / infinity / threshold, boundary fits, two half-line transports, threshold
//!               join, Mellin moments (FP int rho s^-m ds)
pub mod engine;
pub mod graph;
pub mod local;
#[cfg(test)]
mod tests;
