//! Finite-density algebra and admission owners.
//!
//! Occupation coordinates are distinct from required cuts. These owners do not
//! turn a formal cut decomposition into a numerical evaluation certificate.
//! See `docs/finite-density.md` for the current implementation and acceptance status.
pub mod assembly;
mod basis;
pub mod boundary;
pub mod compact;
mod contour;
pub mod flow;
pub mod flow_boundary;
pub mod geometry;
pub mod guarded;
mod input;
pub mod interface;
pub mod massless_contour;
pub mod massless_endpoint;
pub mod measure;
pub mod normalization;
pub mod parametric_endpoint;
pub mod preparation;
pub mod reduction;
pub mod terminal;
pub mod vacuum;

pub use basis::InversePropagatorBasis;
pub use input::{
    DensityEdge, DensityInput, DensityTarget, NumeratorConvention, PreparedDensityInput,
};
