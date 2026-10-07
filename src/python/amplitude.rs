//! Application binding over the native HEPKit tensor and scalar assembly pipeline.
use super::*;
use crate::gg_hg::amplitude::{HiggsJetAmplitude, HiggsJetFormFactors, HiggsJetObservables};
use feynkit_py::{PyFeynmanDiagram, PyModel};
use std::sync::Arc;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Coherent W/Z and infinite-top Higgs-plus-jet observables in the effective model's basis.
/// Construction reuses the supplied model for generation, display, and contraction.
/// Cancellation is checked between native stages; an active HEPKit or Idenso
/// operation completes before acknowledging a cancellation request.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "HiggsJetAmplitude",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyHiggsJetAmplitude {
    inner: Arc<HiggsJetAmplitude>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyHiggsJetAmplitude {
    /// Add the symbolic W, Z and HEFT gggH vertices to an existing Standard Model.
    /// Returns a native HEPKit model retaining its existing particles, parameters
    /// and interactions. Numerical form factors still come from loop evaluation;
    /// conflicting declaration names are rejected rather than overwritten.
    #[staticmethod]
    fn with_form_factor_vertices(py: Python<'_>, model: &PyModel) -> PyResult<PyModel> {
        let model = model.as_model().clone();
        py.detach(|| HiggsJetAmplitude::with_form_factor_vertices(&model))
            .inspect(|_| super::citations::higgs_jet())
            .map(Into::into)
            .map_err(error)
    }

    #[new]
    #[pyo3(signature=(model, *, control=None))]
    fn new(
        py: Python<'_>,
        model: &PyModel,
        control: Option<&PyComputationControl>,
    ) -> PyResult<Self> {
        let model = model.as_model().clone();
        let context = context(control);
        py.detach(|| HiggsJetAmplitude::new(model, &context))
            .inspect(|_| super::citations::higgs_jet())
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(error)
    }
    /// Native effective diagrams sharing the original model instance.
    #[getter]
    fn diagrams(&self) -> Vec<PyFeynmanDiagram> {
        self.inner
            .diagrams()
            .iter()
            .map(|diagram| diagram.as_ref().clone().into())
            .collect()
    }
    /// Exact EW square, interference, and effective square kernels, in that order.
    #[getter]
    fn expressions(&self) -> Vec<PythonExpression> {
        self.inner
            .expressions()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
    #[getter]
    fn coordinate_symbols(&self) -> Vec<PythonExpression> {
        self.inner
            .coordinate_symbols()
            .into_iter()
            .map(|s| Atom::var(s).into())
            .collect()
    }
    /// Evaluate normalized W1..W4 and Z1..Z4 form factors and their absolute allowances.
    /// Model parameters are exact Symbolica expressions; model float defaults are unused.
    /// Coefficients must use the same physical kinematics and W/Z masses as the model inputs.
    /// Includes incoming spin/color averages and final-state sums, with no phase-space factor.
    #[pyo3(signature=(s, t, higgs_mass_squared, w_factors, z_factors, w_errors, z_errors, *, parameters, provenance, digits=20, guard_digits=40, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        py: Python<'_>,
        s: PythonExpression,
        t: PythonExpression,
        higgs_mass_squared: PythonExpression,
        w_factors: Vec<PythonMultiPrecisionComplex>,
        z_factors: Vec<PythonMultiPrecisionComplex>,
        w_errors: Vec<PythonMultiPrecisionFloat>,
        z_errors: Vec<PythonMultiPrecisionFloat>,
        parameters: Expressions,
        provenance: &str,
        digits: u32,
        guard_digits: u32,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyAmplitudeResult> {
        if [
            w_factors.len(),
            z_factors.len(),
            w_errors.len(),
            z_errors.len(),
        ] != [4; 4]
        {
            return Err(InvalidInputError::new_err(
                "four W and four Z coefficients and allowances required",
            ));
        }
        let factors = HiggsJetFormFactors {
            values: w_factors
                .into_iter()
                .chain(z_factors)
                .map(|v| v.0)
                .collect::<Vec<_>>()
                .try_into()
                .map_err(|_| InvalidInputError::new_err("form-factor count"))?,
            absolute_errors: w_errors
                .into_iter()
                .chain(z_errors)
                .map(|v| v.0)
                .collect::<Vec<_>>()
                .try_into()
                .map_err(|_| InvalidInputError::new_err("form-factor error count"))?,
            provenance: provenance.to_owned(),
        };
        let parameters = super::point(parameters).0;
        let context = context(control);
        let started = std::time::Instant::now();
        let inner = py
            .detach(|| {
                self.inner.evaluate(
                    &s.expr,
                    &t.expr,
                    &higgs_mass_squared.expr,
                    &factors,
                    &parameters,
                    digits,
                    guard_digits,
                    &context,
                )
            })
            .map_err(error)?;
        Ok(PyAmplitudeResult {
            inner,
            elapsed: started.elapsed(),
        })
    }
}

/// Physical observables with conditional propagated uncertainties and arithmetic refinement.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AmplitudeResult",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAmplitudeResult {
    inner: HiggsJetObservables,
    elapsed: std::time::Duration,
}

impl PyAmplitudeResult {
    fn results(&self) -> [(&str, &crate::gg_hg::amplitude::HiggsJetObservable); 3] {
        [
            ("electroweak_squared", &self.inner.electroweak_squared),
            ("interference", &self.inner.interference),
            ("effective_squared", &self.inner.effective_squared),
        ]
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAmplitudeResult {
    #[getter]
    fn values(&self) -> BTreeMap<String, PythonMultiPrecisionFloat> {
        self.results()
            .into_iter()
            .map(|(name, value)| {
                (
                    name.to_owned(),
                    PythonMultiPrecisionFloat(value.value.clone()),
                )
            })
            .collect()
    }
    #[getter]
    fn absolute_errors(&self) -> BTreeMap<String, PythonMultiPrecisionFloat> {
        self.results()
            .into_iter()
            .map(|(name, value)| {
                (
                    name.to_owned(),
                    PythonMultiPrecisionFloat(value.absolute_error.clone()),
                )
            })
            .collect()
    }
    #[getter]
    fn verified_relative_digits(&self) -> BTreeMap<String, Option<u32>> {
        self.results()
            .into_iter()
            .map(|(name, value)| (name.to_owned(), value.verified_relative_digits))
            .collect()
    }
    #[getter]
    fn arithmetic_changes(&self) -> BTreeMap<String, PythonMultiPrecisionFloat> {
        self.results()
            .into_iter()
            .map(|(name, value)| {
                (
                    name.to_owned(),
                    PythonMultiPrecisionFloat(value.arithmetic_change.clone()),
                )
            })
            .collect()
    }
    #[getter]
    fn working_bits(&self) -> u32 {
        self.inner.working_bits
    }
    #[getter]
    fn provenance(&self) -> &str {
        &self.inner.input_provenance
    }
    #[getter]
    fn elapsed_nanoseconds(&self) -> u128 {
        self.elapsed.as_nanos()
    }
}
