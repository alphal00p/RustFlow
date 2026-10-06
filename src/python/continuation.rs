//! Exact continuation declarations bound to the existing native route planner.
use super::*;
use crate::contour::PolynomialPrescription;
use crate::transport_cache::PhysicalContinuation;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

fn prescription(value: &str) -> PyResult<crate::Prescription> {
    match value {
        "+i0" => Ok(crate::Prescription::PlusI0),
        "-i0" => Ok(crate::Prescription::MinusI0),
        _ => Err(InvalidInputError::new_err(
            "prescription must be +i0 or -i0",
        )),
    }
}
fn name(value: crate::Prescription) -> &'static str {
    match value {
        crate::Prescription::PlusI0 => "+i0",
        crate::Prescription::MinusI0 => "-i0",
    }
}

/// Exact polynomial i0 declarations for physical threshold continuation.
/// Declarations retain their order. The native connection validates their
/// variables and original domains when this object is bound to a transport.
/// A declaration fixes local sides; it does not establish global monodromy.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ContinuationPrescription",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyContinuationPrescription {
    inner: PhysicalContinuation,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyContinuationPrescription {
    #[new]
    #[pyo3(signature=(prescriptions, *, domain, unprescribed_side="+i0"))]
    fn new(
        prescriptions: Vec<(PythonExpression, String)>,
        domain: String,
        unprescribed_side: &str,
    ) -> PyResult<Self> {
        if domain.trim().is_empty() {
            return Err(InvalidInputError::new_err(
                "physical continuation needs an explicit homotopy domain",
            ));
        }
        Ok(Self {
            inner: PhysicalContinuation {
                prescriptions: prescriptions
                    .into_iter()
                    .map(|(polynomial, side)| {
                        Ok(PolynomialPrescription {
                            polynomial: polynomial.expr,
                            prescription: prescription(&side)?,
                        })
                    })
                    .collect::<PyResult<_>>()?,
                unprescribed_side: prescription(unprescribed_side)?,
                domain,
            },
        })
    }

    #[getter]
    fn prescriptions(&self) -> Vec<(PythonExpression, &'static str)> {
        self.inner
            .prescriptions
            .iter()
            .map(|p| (p.polynomial.clone().into(), name(p.prescription)))
            .collect()
    }
    #[getter]
    fn domain(&self) -> &str {
        &self.inner.domain
    }
    #[getter]
    fn unprescribed_side(&self) -> &'static str {
        name(self.inner.unprescribed_side)
    }
}

pub(super) fn bind<S>(
    flow: crate::RustFlow<S>,
    continuation: Option<&PyContinuationPrescription>,
) -> crate::Result<crate::RustFlow<S>> {
    match continuation {
        Some(value) => flow.with_prescribed_continuation(value.inner.clone()),
        None => Ok(flow),
    }
}
