//! Fixed-precision Taylor-patch transport of a linear system `dy/dx = A(x) y`
//! with exact rational-function coefficients, along a planned path.
//!
//! This is the transport kernel of RustFlow's dispersive mode and of any
//! long, planned transport where the adaptive certified stepper of
//! [`crate::ode`] is too expensive:
//!
//! * **Arithmetic.** MPFR at one working precision (`rug::Float` on real
//!   patches, `rug::Complex` on complex ones) is the default. Symbolica's
//!   precision-tracking `Float` stays available as an option
//!   ([`ArithmeticMode::Tracked`]); its tracked precision is then an extra
//!   accuracy indicator. Accuracy of the fixed-precision transport is assessed
//!   by recomputation instead: the same path at two working precisions (with
//!   the truncation orders of each precision), patch-overlap consistency, and
//!   series-tail estimates ([`accuracy`]).
//! * **Exact coefficients.** Every matrix entry is an exact quotient of integer
//!   polynomials. Each patch shifts them EXACTLY to its centre (integer
//!   synthetic division, [`shift`]) and rounds every normalised coefficient
//!   once, so no digits are lost to cancellation in the shift.
//! * **Row-parallel recurrence.** One patch computes its Taylor coefficients
//!   by the polynomial recurrence `D_g u_g = sum_t N_t y_j(t)` per group
//!   (row, denominator). Each row is owned by one worker and receives the same
//!   operations in the same order whatever the thread count, so the output is
//!   bit-identical for any number of threads ([`patch`]).
//! * **Checkpoints.** The state after every patch is exact; it is written
//!   atomically with fingerprints of system, path, precision, mode and binary
//!   and resumed with any thread count ([`transport`], [`checkpoint`]).
//!
//! Nothing here uses global mutable state; threads are scoped per patch.
pub mod accuracy;
pub mod checkpoint;
pub mod path;
pub mod patch;
pub mod scalar;
pub mod shift;
pub mod system;
pub mod transport;
#[cfg(test)]
mod tests;

pub use accuracy::{DigitsReport, compare_digits};
pub use path::{PatchSpec, PlannedPath, PlanInput};
pub use scalar::ArithmeticMode;
pub use system::RationalSystem;
pub use transport::{PatchRecord, Transport, TransportConfig, TransportOutcome, TransportState};

/// Version of the kernel algorithm. It enters every checkpoint fingerprint:
/// change it whenever the arithmetic or the order of operations changes.
pub const KERNEL_VERSION: &str = "rustflow-patch-kernel-1";

/// Binary digits for `digits` decimal digits plus ten guard digits, as the
/// campaign's in-house tool (`floor((W + 10) log2 10)`; 581 bits at W = 165).
pub fn bits_for_digits(digits: u32) -> u32 {
    ((f64::from(digits) + 10.0) * std::f64::consts::LOG2_10).floor() as u32
}
