//! Published gg → Hg integral definitions, certified basis maps and amplitudes.
pub mod amplitude;
pub mod data;
mod form_factors;
mod plugin_basis;
mod published;
mod workflow;

pub use form_factors::{HiggsJetFormFactorResult, HiggsJetFormFactors};
pub use plugin_basis::{PluginBasisMap, PluginFamilyKind, PluginMapEvidence};
pub use workflow::{HiggsJetConfiguration, HiggsJetIntegralSystem};

pub use published::{
    CanonicalIntegral, CanonicalIntegralHead, CanonicalTerm, PublishedCanonicalBasis,
    PublishedFamily, PublishedFamilyKind, PublishedMeasure, PublishedNormalization, SourceEvidence,
};
