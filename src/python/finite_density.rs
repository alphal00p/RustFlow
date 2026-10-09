//! Native finite-density requests in the embedding Symbolica host.
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

/// Return a JSON string containing a complete native finite-density evaluation.
/// `epsilon` is an exact nonzero rational string, or None for Laurent fitting.
/// Values use the unscaled Euclidean measure and preserve target/cut ordering.
/// Every occupied sector must pass native admission and closure; any unresolved
/// sector raises an exception instead of returning a partial amplitude.
#[pyfunction]
#[pyo3(signature = (input_json, epsilon=None, *, options=None, control=None, start_scale=8))]
pub fn evaluate_finite_density(
    py: Python<'_>,
    input_json: &str,
    epsilon: Option<&str>,
    options: Option<&super::PyEvaluationOptions>,
    control: Option<&super::PyComputationControl>,
    start_scale: u32,
) -> PyResult<String> {
    let input: crate::finite_density::DensityInput = serde_json::from_str(input_json)
        .map_err(|e| super::InvalidInputError::new_err(e.to_string()))?;
    let options = options
        .map(|options| options.inner.clone())
        .unwrap_or_else(|| crate::finite_density::interface::default_options(&input));
    let context = super::context(control);
    let report = py
        .detach(|| {
            crate::finite_density::interface::evaluate_request(
                &input,
                epsilon,
                &options,
                crate::finite_density::interface::default_closure_options(),
                &context,
                start_scale,
            )
        })
        .map_err(super::error)?;
    serde_json::to_string(&report).map_err(|e| super::InvalidInputError::new_err(e.to_string()))
}
