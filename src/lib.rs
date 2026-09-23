//! Symbolica-native auxiliary mass flow.
#![doc = include_str!("../docs/usage.md")]
pub use symbolica;
mod algebra;
pub mod benchmarks;
pub mod boundary;
mod bubble;
pub mod cache;
pub mod engine;
pub mod epsilon;
pub mod error;
pub mod family;
pub mod frobenius;
pub mod ft;
pub mod gaussian;
pub mod integrand;
mod native;
pub mod normalize;
pub mod numeric;
pub mod ode;
pub mod options;
pub mod recursive;
pub mod reduction;
pub mod refine;
pub mod regions;
pub mod tensor;
pub mod vacuum;
pub use engine::{PreparedFlow, evaluate_samples, solve_integrals, solve_prepared};
pub use epsilon::{LaurentExpansion, fit_epsilon};
pub use error::{Error, Result};
pub use family::{Integral, IntegralFamily, KinematicPoint, Propagator};
pub use numeric::{ComplexFloat, Precision};
pub use ode::{BoundaryData, DifferentialSystem, FlowDiagnostics, FlowResult};
pub use options::{
    CancellationToken, FlowOptions, MassMode, Prescription, Progress, RecursionMode, RunContext,
};
pub use reduction::{ReductionBackend, RustRedBackend, TableBackend};
