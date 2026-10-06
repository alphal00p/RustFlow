use super::*;
use crate::algebraic::{AlgebraicKinematicSystem, CanonicalAlgebraicSystem, SquareRoot};
use crate::transport_cache::{
    BoundaryAccuracy, BoundaryIdentity, CachedBoundary, CachedPoint, EpsilonRange, PointKind,
    RootGerm, RootSheet, ScaledDistance,
};
use std::sync::Arc;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

enum Connection {
    Rational(crate::RustFlow),
    Algebraic(crate::RustFlow<AlgebraicKinematicSystem>),
    Canonical(Arc<crate::RustFlow<CanonicalAlgebraicSystem>>),
}
impl Connection {
    fn identity(&self) -> &BoundaryIdentity {
        match self {
            Self::Rational(f) => f.identity(),
            Self::Algebraic(f) => f.identity(),
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
impl PyKinematicTransport {
    pub(super) fn from_canonical(
        flow: Arc<crate::RustFlow<CanonicalAlgebraicSystem>>,
        options: crate::FlowOptions,
    ) -> Self {
        Self {
            connection: Arc::new(Connection::Canonical(flow)),
            options,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyKinematicTransport {
    /// Evaluate finite epsilon coefficient limits in an exact endpoint chart.
    /// Admit the regular matching path and final endpoint approach separately.
    /// With a continuation declaration, matching uses the native prescribed
    /// contour and admit_matching_path admits its declared global homotopy.
    /// Recorded input errors and independent precision/order profiles determine
    /// reusable evidence; this is not symbolic dimensional-sector projection.
    #[pyo3(signature=(cache, route, leading, last, *, admit_matching_path=false, admit_endpoint=false, scales=None, max_lift_dimension=256, series_order=64, constraints=None, control=None))]
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
        constraints: Option<&PyEndpointConstraints>,
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
                    || (source.point.same_coordinates(target)?
                        && source.point.root_germ() == target.root_germ()))
            },
        };
        let admission =
            |_: &CachedBoundary, _: &crate::singular_endpoint::EndpointChart| Ok(admit_endpoint);
        let route_admission =
            |_: &CachedBoundary, _: &CachedPoint, _: &crate::physical_transport::PhysicalRoute| {
                Ok(admit_matching_path)
            };
        let prescribed = self.connection.identity().physical_continuation().is_some();
        let started = std::time::Instant::now();
        py.detach(|| {
            cache
                .access(|cache| {
                    macro_rules! evaluate {
                        ($flow:expr) => {
                            match (prescribed, constraints) {
                                (true, Some(asserted)) => $flow
                                    .evaluate_prescribed_constrained_endpoint(
                                        cache,
                                        &request,
                                        &asserted.constraints,
                                        &self.options,
                                        &context,
                                        &policy,
                                        &route_admission,
                                        &admission,
                                    ),
                                (false, Some(asserted)) => $flow.evaluate_constrained_endpoint(
                                    cache,
                                    &request,
                                    &asserted.constraints,
                                    &self.options,
                                    &context,
                                    &policy,
                                    &admission,
                                ),
                                (true, None) => $flow.evaluate_prescribed_endpoint(
                                    cache,
                                    &request,
                                    &self.options,
                                    &context,
                                    &policy,
                                    &route_admission,
                                    &admission,
                                ),
                                (false, None) => $flow.evaluate_endpoint(
                                    cache,
                                    &request,
                                    &self.options,
                                    &context,
                                    &policy,
                                    &admission,
                                ),
                            }
                        };
                    }
                    match self.connection.as_ref() {
                        Connection::Rational(flow) => evaluate!(flow),
                        Connection::Canonical(flow) => evaluate!(flow),
                        Connection::Algebraic(flow) => evaluate!(flow),
                    }
                })
                .and_then(PyEndpointResult::from_result)
        })
        .map(|result| result.timed(started))
        .map_err(error)
    }
    /// Construct a general physical connection, optionally with named square roots.
    /// Registered roots require explicit source and destination sheet signs.
    /// Matrix entries may contain epsilon-independent and higher epsilon terms;
    /// their admission and expansion use the native connection implementation.
    #[new]
    #[pyo3(signature=(epsilon, derivatives, basis, normalization, *, branch_domain, roots=None, options=None, nonzero_conditions=None, continuation=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        epsilon: PythonExpression,
        derivatives: HashMap<PythonExpression, Vec<Vec<PythonExpression>>>,
        basis: Vec<PythonExpression>,
        normalization: PythonExpression,
        branch_domain: &str,
        roots: Option<Expressions>,
        options: Option<&PyEvaluationOptions>,
        nonzero_conditions: Option<Vec<PythonExpression>>,
        continuation: Option<&PyContinuationPrescription>,
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
        let roots = coordinates(roots.unwrap_or_default())?
            .into_iter()
            .map(|(symbol, radicand)| SquareRoot { symbol, radicand })
            .collect::<Vec<_>>();
        let connection = py
            .detach(|| {
                if roots.is_empty() {
                    let flow = crate::RustFlow::with_conditions(
                        system,
                        &basis,
                        &normalization.expr,
                        options.prescription,
                        branch_domain,
                        &conditions,
                    )?;
                    Ok(Connection::Rational(super::continuation::bind(
                        flow,
                        continuation,
                    )?))
                } else {
                    let flow = crate::RustFlow::with_algebraic_conditions(
                        AlgebraicKinematicSystem {
                            epsilon: system.epsilon,
                            derivatives: system.derivatives,
                            roots,
                        },
                        &basis,
                        &normalization.expr,
                        options.prescription,
                        branch_domain,
                        &conditions,
                    )?;
                    Ok(Connection::Algebraic(super::continuation::bind(
                        flow,
                        continuation,
                    )?))
                }
            })
            .map_err(error)?;
        Ok(Self {
            connection: Arc::new(connection),
            options,
        })
    }
    /// Construct dY = epsilon sum(M_a dlog(letter_a))Y without dense physical assembly.
    /// Named square roots retain separate generators and explicit endpoint sheet signs.
    #[staticmethod]
    #[pyo3(signature=(epsilon, variables, letters, matrices, basis, normalization, *, branch_domain, roots=None, options=None, nonzero_conditions=None, continuation=None))]
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
        continuation: Option<&PyContinuationPrescription>,
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
                let flow = crate::RustFlow::with_canonical_conditions(
                    system,
                    &basis,
                    &normalization.expr,
                    options.prescription,
                    branch_domain,
                    &conditions,
                )?;
                super::continuation::bind(flow, continuation)
            })
            .map_err(error)?;
        Ok(Self::from_canonical(Arc::new(flow), options))
    }
    #[getter]
    fn identity(&self) -> &str {
        self.connection.identity().key()
    }
    #[getter]
    fn dimension(&self) -> usize {
        self.connection.identity().dimension()
    }
    /// Native root declarations; discrete sheet signs live on each boundary/query.
    #[getter]
    fn roots(&self) -> Expressions {
        self.connection
            .identity()
            .roots()
            .iter()
            .map(|root| (Atom::var(root.symbol).into(), root.radicand.clone().into()))
            .collect()
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
    /// admitted regular affine path or, with continuation declarations, an
    /// explicitly admitted prescribed homotopy. These admissions are separate.
    #[pyo3(signature=(cache, destination, leading, last, *, root_sheets=None, admit_straight_path=false, admit_prescribed_path=false, scales=None, control=None))]
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
        admit_prescribed_path: bool,
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
        let prescribed = identity.physical_continuation().is_some();
        let policy = ScaledDistance {
            scales,
            admissible: |source: &CachedBoundary, target: &CachedPoint| {
                if source.point.same_coordinates(target)?
                    && source.point.root_germ() == target.root_germ()
                {
                    return Ok(true);
                }
                if prescribed {
                    return Ok(admit_prescribed_path);
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
        let route_admission =
            |_: &CachedBoundary, _: &CachedPoint, _: &crate::physical_transport::PhysicalRoute| {
                Ok(admit_prescribed_path)
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
                        if prescribed {
                            return flow.evaluate_prescribed_to(
                                cache,
                                &destination,
                                range,
                                &self.options,
                                &context,
                                &policy,
                                &route_admission,
                            );
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
                    Connection::Canonical(flow) if prescribed => flow.evaluate_prescribed_to(
                        cache,
                        &destination,
                        germ.as_ref(),
                        range,
                        &self.options,
                        &context,
                        &policy,
                        &route_admission,
                    ),
                    Connection::Algebraic(flow) => {
                        let germ = germ.as_ref().ok_or_else(|| {
                            crate::Error::InvalidInput(
                                "registered-root transport requires destination root_sheets".into(),
                            )
                        })?;
                        if prescribed {
                            flow.evaluate_prescribed_to(
                                cache,
                                &destination,
                                germ,
                                range,
                                &self.options,
                                &context,
                                &policy,
                                &route_admission,
                            )
                        } else {
                            flow.evaluate_to(
                                cache,
                                &destination,
                                germ,
                                range,
                                &self.options,
                                &context,
                                &policy,
                            )
                        }
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

#[cfg(test)]
mod sharing_tests {
    use super::*;

    #[test]
    fn higgs_adapter_shares_the_existing_canonical_owner() -> crate::Result<()> {
        let system = crate::gg_hg::HiggsJetIntegralSystem::load(
            crate::gg_hg::PluginFamilyKind::Planar,
            "portable_owner_test",
        )?;
        let connection = system.transport.clone();
        let adapter =
            PyKinematicTransport::from_canonical(connection.clone(), crate::FlowOptions::default());
        match adapter.connection.as_ref() {
            Connection::Canonical(shared) => assert!(Arc::ptr_eq(shared, &connection)),
            _ => panic!("canonical owner changed"),
        }
        drop(system);
        assert_eq!(
            adapter.connection.identity().key(),
            connection.identity().key()
        );
        Ok(())
    }
}
