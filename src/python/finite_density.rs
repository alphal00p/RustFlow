//! Algebraic finite-density input preparation in the embedding Symbolica host.
use pyo3::prelude::*;

/// Validate graph/charge evidence and convert polynomial targets, returning the
/// same JSON report as `rustflow finite-density-prepare`. No numerical result is
/// returned by this preparation entry point.
#[pyfunction]
#[pyo3(signature = (input_json, cut_budget=65536))]
pub fn prepare_finite_density(input_json: &str, cut_budget: usize) -> PyResult<String> {
    let input: crate::finite_density::DensityInput = serde_json::from_str(input_json)
        .map_err(|e| super::InvalidInputError::new_err(e.to_string()))?;
    let report = input
        .prepare()
        .and_then(|p| p.preparation_report(cut_budget))
        .map_err(super::error)?;
    serde_json::to_string(&report).map_err(|e| super::InvalidInputError::new_err(e.to_string()))
}
