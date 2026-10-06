use super::*;
use crate::algebraic::{CanonicalAlgebraicSystem, SquareRoot};
use crate::transport_cache::{
    BoundaryAccuracy, BoundaryIdentity, CachedBoundary, CachedPoint, EpsilonRange, PointKind,
    RootGerm, RootSheet, ScaledDistance,
};
use std::sync::Arc;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

enum Connection {
    Rational(crate::RustFlow),
    Canonical(crate::RustFlow<CanonicalAlgebraicSystem>),
}
impl Connection {
    fn identity(&self) -> &BoundaryIdentity {
        match self {
            Self::Rational(f) => f.identity(),
            Self::Canonical(f) => f.identity(),
        }
    }
}
pub(super) fn germ(sheets: Option<HashMap<PythonExpression, i8>>) -> PyResult<Option<RootGerm>> {
    sheets
        .map(|sheets| {
            sheets
                .into_iter()
                .map(|(key, sign)| {
                    Ok((
                        symbol(&key)?,
                        match sign {
                            1 => RootSheet::Principal,
                            -1 => RootSheet::Opposite,
                            _ => {
                                return Err(InvalidInputError::new_err(
                                    "root sheet signs must be +1 or -1",
                                ));
                            }
                        },
                    ))
                })
                .collect::<PyResult<BTreeMap<_, _>>>()
                .map(|sheets| RootGerm { sheets })
        })
        .transpose()
}

/// Transport a common master basis in physical invariants, reusing verified cached points.
/// The declared branch domain and explicit straight-path admission are part of the
/// physical continuation contract; distance alone never establishes a sheet.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "KinematicTransport",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyKinematicTransport {
    connection: Arc<Connection>,
    options: crate::FlowOptions,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyKinematicTransport {
    /// Evaluate finite epsilon coefficient limits in an exact endpoint chart.
    /// Admit the regular matching path and final endpoint approach separately.
    /// Recorded input errors and independent precision/order profiles determine
    /// reusable evidence; this is not symbolic dimensional-sector projection.
    #[pyo3(signature=(cache, route, leading, last, *, admit_matching_path=false, admit_endpoint=false, scales=None, max_lift_dimension=256, series_order=64, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate_endpoint(
        &self,
        py: Python<'_>,
        cache: &PyBoundaryCache,
        route: &PyEndpointRoute,
        leading: i32,
        last: i32,
        admit_matching_path: bool,
        admit_endpoint: bool,
        scales: Option<Expressions>,
        max_lift_dimension: usize,
        series_order: usize,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyEndpointResult> {
        let request = crate::singular_endpoint::EndpointRequest {
            chart: route.chart.clone(),
            range: EpsilonRange::new(leading, last).map_err(error)?,
            options: crate::singular_endpoint::EndpointOptions {
                max_lift_dimension,
                series_order,
            },
        };
        let context = context(control);
        let policy = ScaledDistance {
            scales: coordinates(scales.unwrap_or_default())?,
            admissible: |source: &CachedBoundary, target: &CachedPoint| {
                Ok(admit_matching_path
                    || (source.point.restart_coordinates()? == target.restart_coordinates()?
                        && source.point.root_germ() == target.root_germ()))
            },
        };
        let admission =
            |_: &CachedBoundary, _: &crate::singular_endpoint::EndpointChart| Ok(admit_endpoint);
        let started = std::time::Instant::now();
        py.detach(|| {
            cache
                .access(|cache| match self.connection.as_ref() {
                    Connection::Rational(flow) => flow.evaluate_endpoint(
                        cache,
                        &request,
                        &self.options,
                        &context,
                        &policy,
                        &admission,
                    ),
                    Connection::Canonical(flow) => flow.evaluate_endpoint(
                        cache,
                        &request,
                        &self.options,
                        &context,
                        &policy,
                        &admission,
                    ),
                })
                .and_then(PyEndpointResult::from_result)
        })
        .map(|result| result.timed(started))
        .map_err(error)
    }
    #[new]
    #[pyo3(signature=(epsilon, derivatives, basis, normalization, *, branch_domain, options=None, nonzero_conditions=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        epsilon: PythonExpression,
        derivatives: HashMap<PythonExpression, Vec<Vec<PythonExpression>>>,
        basis: Vec<PythonExpression>,
        normalization: PythonExpression,
        branch_domain: &str,
        options: Option<&PyEvaluationOptions>,
        nonzero_conditions: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        let options = super::options(options);
        options.validate().map_err(error)?;
        let system = crate::kinematics::KinematicSystem {
            epsilon: symbol(&epsilon)?,
            derivatives: derivatives
                .into_iter()
                .map(|(key, matrix)| {
                    Ok((
                        symbol(&key)?,
                        matrix
                            .into_iter()
                            .map(|r| r.into_iter().map(|x| x.expr).collect())
                            .collect(),
                    ))
                })
                .collect::<PyResult<_>>()?,
        };
        let basis = basis.into_iter().map(|x| x.expr).collect::<Vec<_>>();
        let conditions = nonzero_conditions
            .unwrap_or_default()
            .into_iter()
            .map(|x| x.expr)
            .collect::<Vec<_>>();
        let flow = py
            .detach(|| {
                crate::RustFlow::with_conditions(
                    system,
                    &basis,
                    &normalization.expr,
                    options.prescription,
                    branch_domain,
                    &conditions,
                )
            })
            .map_err(error)?;
        Ok(Self {
            connection: Arc::new(Connection::Rational(flow)),
            options,
        })
    }
    /// Construct dY = epsilon sum(M_a dlog(letter_a))Y without dense physical assembly.
    /// Named square roots retain separate generators and explicit endpoint sheet signs.
    #[staticmethod]
    #[pyo3(signature=(epsilon, variables, letters, matrices, basis, normalization, *, branch_domain, roots=None, options=None, nonzero_conditions=None))]
    #[allow(clippy::too_many_arguments)]
    fn canonical(
        py: Python<'_>,
        epsilon: PythonExpression,
        variables: Vec<PythonExpression>,
        letters: Vec<PythonExpression>,
        matrices: Vec<Vec<Vec<PythonExpression>>>,
        basis: Vec<PythonExpression>,
        normalization: PythonExpression,
        branch_domain: &str,
        roots: Option<Expressions>,
        options: Option<&PyEvaluationOptions>,
        nonzero_conditions: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        let options = super::options(options);
        options.validate().map_err(error)?;
        let epsilon = symbol(&epsilon)?;
        let variables = variables.iter().map(symbol).collect::<PyResult<Vec<_>>>()?;
        let letters = letters.into_iter().map(|x| x.expr).collect::<Vec<_>>();
        let matrices = matrices
            .into_iter()
            .map(|m| {
                m.into_iter()
                    .map(|r| r.into_iter().map(|x| x.expr).collect())
                    .collect()
            })
            .collect::<Vec<_>>();
        let roots = coordinates(roots.unwrap_or_default())?
            .into_iter()
            .map(|(symbol, radicand)| SquareRoot { symbol, radicand })
            .collect();
        let basis = basis.into_iter().map(|x| x.expr).collect::<Vec<_>>();
        let conditions = nonzero_conditions
            .unwrap_or_default()
            .into_iter()
            .map(|x| x.expr)
            .collect::<Vec<_>>();
        let flow = py
            .detach(|| {
                let system =
                    CanonicalAlgebraicSystem::new(epsilon, &variables, &letters, &matrices, roots)?;
                crate::RustFlow::with_canonical_conditions(
                    system,
                    &basis,
                    &normalization.expr,
                    options.prescription,
                    branch_domain,
                    &conditions,
                )
            })
            .map_err(error)?;
        Ok(Self {
            connection: Arc::new(Connection::Canonical(flow)),
            options,
        })
    }
    #[getter]
    fn identity(&self) -> &str {
        self.connection.identity().key()
    }
    #[getter]
    fn dimension(&self) -> usize {
        self.connection.identity().dimension()
    }
    #[getter]
    fn nonzero_conditions(&self) -> Vec<PythonExpression> {
        self.connection
            .identity()
            .nonzero_conditions()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
    /// Insert caller-established boundary evidence. Declared accuracy is not certified here.
    /// Coefficients and absolute error estimates are epsilon-major arrays.
    #[pyo3(signature=(cache, point, coefficients, leading, *, verified_digits, comparison_errors, provenance, root_sheets=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_boundary(
        &self,
        py: Python<'_>,
        cache: &PyBoundaryCache,
        point: Expressions,
        coefficients: Vec<Vec<PythonMultiPrecisionComplex>>,
        leading: i32,
        verified_digits: u32,
        comparison_errors: Vec<Vec<PythonMultiPrecisionFloat>>,
        provenance: &str,
        root_sheets: Option<HashMap<PythonExpression, i8>>,
    ) -> PyResult<PyTransportResult> {
        let mut point = CachedPoint::Exact(coordinates(point)?);
        if let Some(germ) = germ(root_sheets)? {
            point = point.with_root_germ(germ).map_err(error)?;
        }
        let coefficients = coefficients
            .into_iter()
            .map(|r| r.into_iter().map(|v| v.0).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let bits = coefficients
            .iter()
            .flatten()
            .map(|v| v.re.prec().min(v.im.prec()))
            .min()
            .ok_or_else(|| InvalidInputError::new_err("nonempty boundary required"))?;
        let count = i32::try_from(coefficients.len())
            .map_err(|_| ResourceLimitError::new_err("epsilon range overflow"))?;
        let last = leading
            .checked_add(count - 1)
            .ok_or_else(|| ResourceLimitError::new_err("epsilon range overflow"))?;
        let accuracy = BoundaryAccuracy::supplied(
            verified_digits,
            bits,
            comparison_errors
                .into_iter()
                .map(|r| r.into_iter().map(|v| v.0).collect())
                .collect(),
            provenance,
        )
        .map_err(error)?;
        let boundary = CachedBoundary {
            identity: self.connection.identity().clone(),
            point,
            kind: PointKind::Physical,
            range: EpsilonRange::new(leading, last).map_err(error)?,
            coefficients,
            accuracy,
        };
        py.detach(|| cache.access(|cache| cache.insert(boundary.clone())))
            .map_err(error)?;
        PyTransportResult::cached(boundary).map_err(error)
    }
    /// A true cache hit needs no admission. Other points require an explicitly
    /// admitted regular affine path within the constructor's branch domain.
    #[pyo3(signature=(cache, destination, leading, last, *, root_sheets=None, admit_straight_path=false, scales=None, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        py: Python<'_>,
        cache: &PyBoundaryCache,
        destination: Expressions,
        leading: i32,
        last: i32,
        root_sheets: Option<HashMap<PythonExpression, i8>>,
        admit_straight_path: bool,
        scales: Option<Expressions>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyTransportResult> {
        let destination = coordinates(destination)?;
        let germ = germ(root_sheets)?;
        let scales = coordinates(scales.unwrap_or_default())?;
        let context = context(control);
        let range = EpsilonRange::new(leading, last).map_err(error)?;
        let p = crate::Precision::decimal(self.options.digits + self.options.guard_digits)
            .map_err(error)?;
        let identity = self.connection.identity();
        let policy = ScaledDistance {
            scales,
            admissible: |source: &CachedBoundary, target: &CachedPoint| {
                if source.point.restart_coordinates()? == target.restart_coordinates()?
                    && source.point.root_germ() == target.root_germ()
                {
                    return Ok(true);
                }
                if !admit_straight_path {
                    return Ok(false);
                }
                identity.conditions_admit_straight_path(
                    &source.point,
                    target,
                    p,
                    self.options.digits,
                )
            },
        };
        let started = std::time::Instant::now();
        py.detach(|| {
            cache
                .access(|cache| match self.connection.as_ref() {
                    Connection::Rational(flow) => {
                        if germ.is_some() {
                            return Err(crate::Error::InvalidInput(
                                "rational system has no root sheets".into(),
                            ));
                        }
                        flow.evaluate_to(
                            cache,
                            &destination,
                            range,
                            &self.options,
                            &context,
                            &policy,
                        )
                    }
                    Connection::Canonical(flow) => flow.evaluate_to(
                        cache,
                        &destination,
                        germ.as_ref(),
                        range,
                        &self.options,
                        &context,
                        &policy,
                    ),
                })
                .and_then(PyTransportResult::physical)
        })
        .map(|result| result.timed(started))
        .map_err(error)
    }
}
