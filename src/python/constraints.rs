use super::*;
use crate::asymptotic::{AsymptoticSelector, ExactAsymptoticConstraints, ExactAsymptoticRelation};
use crate::frobenius::ExactFrobeniusLimits;
use crate::singular_endpoint::EndpointConstraints;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Select a power/log coefficient in a declared endpoint chart.
/// For epsilon hierarchies, component is epsilon_offset * basis_size + component.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AsymptoticCoefficient",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAsymptoticCoefficient {
    selector: AsymptoticSelector,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAsymptoticCoefficient {
    #[new]
    #[pyo3(signature=(component, power, log_power=0))]
    fn new(component: usize, power: PythonExpression, log_power: usize) -> Self {
        Self {
            selector: AsymptoticSelector {
                component,
                power: power.expr,
                log_power,
            },
        }
    }
    #[getter]
    fn component(&self) -> usize {
        self.selector.component
    }
    #[getter]
    fn power(&self) -> PythonExpression {
        self.selector.power.clone().into()
    }
    #[getter]
    fn log_power(&self) -> usize {
        self.selector.log_power
    }
}

/// Assert sum(weight * selected coefficient) = value using exact expressions.
/// This is supplied physical information, never inferred from numerical errors.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AsymptoticRelation",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAsymptoticRelation {
    relation: ExactAsymptoticRelation,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAsymptoticRelation {
    #[new]
    fn new(
        terms: Vec<(PyAsymptoticCoefficient, PythonExpression)>,
        value: PythonExpression,
    ) -> Self {
        Self {
            relation: ExactAsymptoticRelation {
                terms: terms
                    .into_iter()
                    .map(|(s, w)| (s.selector, w.expr))
                    .collect(),
                value: value.expr,
            },
        }
    }
    #[getter]
    fn terms(&self) -> Vec<(PyAsymptoticCoefficient, PythonExpression)> {
        self.relation
            .terms
            .iter()
            .map(|(s, w)| {
                (
                    PyAsymptoticCoefficient {
                        selector: s.clone(),
                    },
                    w.clone().into(),
                )
            })
            .collect()
    }
    #[getter]
    fn value(&self) -> PythonExpression {
        self.relation.value.clone().into()
    }
}

/// Explicit exact endpoint relations, their provenance and bounded work settings.
/// Weights and values are Gaussian-rational; selected powers are real rational.
/// Numerical matching and uncertainty remain separate from these assumptions.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EndpointConstraints",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEndpointConstraints {
    pub(crate) constraints: EndpointConstraints,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEndpointConstraints {
    #[new]
    #[pyo3(signature=(relations, provenance, *, max_dimension=64, max_order=256, max_coefficient_bits=65536, max_scalar_cells=4000000))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        relations: Vec<PyAsymptoticRelation>,
        provenance: &str,
        max_dimension: usize,
        max_order: usize,
        max_coefficient_bits: u64,
        max_scalar_cells: usize,
    ) -> PyResult<Self> {
        let constraints = EndpointConstraints {
            asymptotic: ExactAsymptoticConstraints {
                relations: relations.into_iter().map(|r| r.relation).collect(),
                provenance: provenance.into(),
            },
            limits: ExactFrobeniusLimits {
                max_dimension,
                max_order,
                max_coefficient_bits,
                max_scalar_cells,
            },
        };
        py.detach(|| {
            // This lower bound validates the declarations. The endpoint owner's
            // preflight checks the actual complete epsilon/basis dimension again.
            let dimension = constraints
                .asymptotic
                .relations
                .iter()
                .flat_map(|r| &r.terms)
                .map(|(s, _)| s.component)
                .max()
                .unwrap_or(0)
                .checked_add(1)
                .ok_or_else(|| crate::Error::Limit("asymptotic component index overflow".into()))?;
            constraints.asymptotic.validate(
                dimension,
                &constraints.limits,
                &crate::RunContext::default(),
            )
        })
        .map_err(error)?;
        Ok(Self { constraints })
    }
    #[getter]
    fn relations(&self) -> Vec<PyAsymptoticRelation> {
        self.constraints
            .asymptotic
            .relations
            .iter()
            .cloned()
            .map(|relation| PyAsymptoticRelation { relation })
            .collect()
    }
    #[getter]
    fn provenance(&self) -> &str {
        &self.constraints.asymptotic.provenance
    }
    #[getter]
    fn max_dimension(&self) -> usize {
        self.constraints.limits.max_dimension
    }
    #[getter]
    fn max_order(&self) -> usize {
        self.constraints.limits.max_order
    }
    #[getter]
    fn max_coefficient_bits(&self) -> u64 {
        self.constraints.limits.max_coefficient_bits
    }
    #[getter]
    fn max_scalar_cells(&self) -> usize {
        self.constraints.limits.max_scalar_cells
    }
}
