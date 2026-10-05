//! Symbolica-native auxiliary mass flow.
#![doc = include_str!("../docs/usage.md")]
pub use symbolica;
mod algebra;
pub mod algebraic;
pub mod asymptotic;
pub mod benchmarks;
pub mod boundary;
mod bubble;
pub mod cache;
pub mod contour;
pub mod cuts;
pub mod diffexp;
pub mod engine;
pub mod epsilon;
pub mod epsilon_flow;
pub mod epsilon_shearing;
pub mod error;
pub mod family;
mod fixed_series;
pub mod frobenius;
pub mod ft;
pub mod gaussian;
pub mod hepkit;
mod integer_shearing;
pub mod integrand;
pub mod kinematic_derivative;
pub mod kinematics;
pub mod linear;
pub mod local_coordinates;
mod native;
pub mod normalize;
pub mod numeric;
pub mod ode;
pub mod options;
pub mod phase_space;
mod physical_conditions;
pub mod physical_family;
pub mod physical_targets;
pub mod physical_transport;
pub mod recursive;
pub mod reduction;
pub mod refine;
pub mod regions;
mod root_path;
pub mod tensor;
pub mod transport_cache;
pub mod vacuum;
mod vacuum_peel;
pub use asymptotic::AsymptoticConstraint;
pub use engine::{
    PreparedFlow, evaluate_samples, solve_integral_combinations, solve_integrals, solve_prepared,
};
pub use epsilon::{LaurentExpansion, fit_epsilon};
pub use epsilon_flow::EpsilonShearedFlow;
pub use epsilon_shearing::EpsilonShearing;
pub use error::{Error, Result};
pub use family::{Integral, IntegralFamily, KinematicPoint, Propagator};
pub use local_coordinates::LocalCoordinate;
pub use numeric::{ComplexFloat, Precision};
pub use ode::{BoundaryData, DifferentialSystem, FlowDiagnostics, FlowResult};
pub use options::{
    CancellationToken, FlowOptions, MassMode, Prescription, Progress, RecursionMode, RunContext,
    StepSizeStrategy,
};
pub use physical_family::PreparedPhysicalFamily;
pub use physical_targets::ProjectedLaurentExpansion;
pub use physical_transport::RustFlow;
pub use reduction::{ReductionBackend, RustRedBackend, ScopedTableBackend, TableBackend};
pub use transport_cache::RustFlowCache;
