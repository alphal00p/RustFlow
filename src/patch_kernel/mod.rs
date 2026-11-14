//! Fixed-precision Taylor-patch transport of a linear system `dy/dx = A(x) y`
//! with exact rational-function coefficients, along a planned path.
//!
//! RustFlow has two transports of a rational linear system:
//!
//! * [`crate::ode::CompiledSystem::transport`]: ADAPTIVE and CERTIFIED (step
//!   control with residual bounds on Symbolica floats, precision refinement);
//!   the default of the ordinary auxiliary-mass flow.
//! * this module: PLANNED-PATH FIXED-PRECISION transport (patches given in
//!   advance, MPFR at one working precision, accuracy assessed by
//!   recomputation). The default for dispersive transport and for any long
//!   planned path. Build its system from the solver's own types with
//!   [`RationalSystem::from_differential_system`] or
//!   [`crate::ode::CompiledSystem::patch_kernel_system`].
//!
//! The library takes a cancellation token only; process signal handling is
//! the caller's (see `examples/patch_transport.rs`).
//!
//! Properties:
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
//!   (row, denominator), scheduled dynamically over fixed chunks; every value
//!   receives the same operations in the same order whichever thread computes
//!   it, so the output is bit-identical for any number of threads ([`patch`]).
//! * **Checkpoints.** The state after every patch is exact; it is written
//!   atomically with fingerprints of system, path, precision, mode and
//!   algorithm, never over another run's checkpoint, under an advisory lock,
//!   and resumed with any thread count ([`transport`], [`checkpoint`]).
//! * **Fail fast.** Non-finite values and divergent patches (series tail above
//!   the bound) stop the run with an error, the state left at the last good
//!   patch.
//!
//! Nothing here uses global mutable state; threads are scoped per patch.
// Index loops in the numeric kernels are deliberate: the order of operations is part of the bit-identical (thread-count independent)
// results. Negated comparisons (`!(x <= bound)`) are deliberate: they treat NaN as a failure.
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]
pub mod accuracy;
pub mod checkpoint;
pub mod patch;
pub mod path;
pub mod scalar;
pub mod shift;
pub mod system;
#[cfg(test)]
mod tests;
pub mod transport;

pub use accuracy::{DigitsReport, compare_digits};
pub use path::{PatchSpec, PlanInput, PlannedPath};
pub use scalar::ArithmeticMode;
pub use system::RationalSystem;
pub use transport::{PatchRecord, Transport, TransportConfig, TransportOutcome, TransportState};

/// Version of the kernel algorithm. It enters every checkpoint fingerprint:
/// change it whenever the arithmetic or the order of operations changes.
pub const KERNEL_VERSION: &str = "rustflow-patch-kernel-2";

/// Automatic algorithm identity (build.rs): digest of `src/patch_kernel` and
/// of the resolved dependency graph. Always part of the checkpoint check.
pub const KERNEL_ALGORITHM: &str = env!("PATCH_KERNEL_SOURCE_DIGEST");

/// Binary digits for `digits` decimal digits plus ten guard digits, as the
/// campaign's in-house tool (`floor((W + 10) log2 10)`; 581 bits at W = 165).
pub fn bits_for_digits(digits: u32) -> u32 {
    ((f64::from(digits) + 10.0) * std::f64::consts::LOG2_10).floor() as u32
}
