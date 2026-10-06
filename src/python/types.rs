use super::*;
use std::collections::VecDeque;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Requested accuracy and resource limits. Working precision is not achieved accuracy.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EvaluationOptions",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEvaluationOptions {
    pub(crate) inner: crate::FlowOptions,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEvaluationOptions {
    #[new]
    #[pyo3(signature = (*, digits=20, guard_digits=40, series_order=80, max_steps=1000, workers=1, dimension=4, recursion="auxiliary_mass", prescription="+i0", mass_mode="automatic", refine_basis=false, skip_reduction=false, sampled_reduction=true, max_precision_attempts=3, max_boundary_attempts=8, cache_directory=None, local_coordinate="identity", sample_cache_directory=None, reuse_samples=true, pade_degree=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        digits: u32,
        guard_digits: u32,
        series_order: usize,
        max_steps: usize,
        workers: usize,
        dimension: i64,
        recursion: &str,
        prescription: &str,
        mass_mode: &str,
        refine_basis: bool,
        skip_reduction: bool,
        sampled_reduction: bool,
        max_precision_attempts: usize,
        max_boundary_attempts: usize,
        cache_directory: Option<PathBuf>,
        local_coordinate: &str,
        sample_cache_directory: Option<PathBuf>,
        reuse_samples: bool,
        pade_degree: Option<usize>,
    ) -> PyResult<Self> {
        let inner = crate::FlowOptions {
            digits,
            guard_digits,
            series_order,
            pade: pade_degree.map(|degree| crate::PadeOptions {
                degree,
                ..Default::default()
            }),
            max_steps,
            workers,
            dimension,
            recursion: match recursion {
                "auxiliary_mass" => crate::RecursionMode::Amf,
                "propagator_combination" => crate::RecursionMode::Ft,
                _ => {
                    return Err(InvalidInputError::new_err(
                        "recursion must be auxiliary_mass or propagator_combination",
                    ));
                }
            },
            prescription: match prescription {
                "+i0" => crate::Prescription::PlusI0,
                "-i0" => crate::Prescription::MinusI0,
                _ => {
                    return Err(InvalidInputError::new_err(
                        "prescription must be +i0 or -i0",
                    ));
                }
            },
            mass_mode: match mass_mode {
                "automatic" => crate::MassMode::Auto,
                "all" => crate::MassMode::All,
                "mass" => crate::MassMode::Mass,
                "propagator" => crate::MassMode::Propagator,
                "branch" => crate::MassMode::Branch,
                "loop" => crate::MassMode::Loop,
                _ => {
                    return Err(InvalidInputError::new_err(
                        "unknown auxiliary mass placement",
                    ));
                }
            },
            local_coordinate: match local_coordinate {
                "identity" => crate::LocalCoordinate::Identity,
                "balanced_mobius" => crate::LocalCoordinate::BalancedMobius,
                _ => {
                    return Err(InvalidInputError::new_err(
                        "local_coordinate must be identity or balanced_mobius",
                    ));
                }
            },
            refine_basis,
            skip_reduction,
            sampled_reduction,
            max_precision_attempts,
            max_boundary_attempts,
            cache_directory,
            sample_cache_directory,
            reuse_samples,
            ..Default::default()
        };
        inner.validate().map_err(error)?;
        Ok(Self { inner })
    }
    #[getter]
    fn digits(&self) -> u32 {
        self.inner.digits
    }
    #[getter]
    fn guard_digits(&self) -> u32 {
        self.inner.guard_digits
    }
    #[getter]
    fn series_order(&self) -> usize {
        self.inner.series_order
    }
    /// Optional rational approximation degree; absent means Taylor transport.
    #[getter]
    fn pade_degree(&self) -> Option<usize> {
        self.inner.pade.as_ref().map(|p| p.degree)
    }
    #[getter]
    fn workers(&self) -> usize {
        self.inner.workers
    }
    #[getter]
    fn dimension(&self) -> i64 {
        self.inner.dimension
    }
    fn __repr__(&self) -> String {
        format!(
            "EvaluationOptions(digits={}, guard_digits={}, series_order={}, workers={})",
            self.inner.digits, self.inner.guard_digits, self.inner.series_order, self.inner.workers
        )
    }
}

/// Thread-safe cancellation and progress polling for a running calculation.
/// Run an evaluator in a worker thread and poll this object from the UI thread.
/// At most the 4096 most recent progress events are retained.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ComputationControl",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyComputationControl {
    cancellation: crate::CancellationToken,
    events: Arc<Mutex<VecDeque<String>>>,
}
impl PyComputationControl {
    pub(crate) fn context(&self) -> crate::RunContext {
        let events = self.events.clone();
        crate::RunContext {
            cancellation: self.cancellation.clone(),
            progress: Some(Arc::new(move |event| {
                if let Ok(mut queue) = events.lock() {
                    if queue.len() == 4096 {
                        queue.pop_front();
                    }
                    queue.push_back(format!("{event:?}"));
                }
            })),
        }
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyComputationControl {
    #[new]
    fn new() -> Self {
        Self::default()
    }
    fn cancel(&self) {
        self.cancellation.cancel();
    }
    #[getter]
    fn cancelled(&self) -> bool {
        self.cancellation.check().is_err()
    }
    /// Drain queued progress events. Cancellation remains set for this control.
    fn poll(&self) -> PyResult<Vec<String>> {
        Ok(self
            .events
            .lock()
            .map_err(|_| EvaluationError::new_err("progress lock poisoned"))?
            .drain(..)
            .collect())
    }
}

/// Epsilon coefficients with independent-fit evidence and absolute comparison errors.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LaurentExpansion",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLaurentExpansion {
    pub(crate) inner: crate::LaurentExpansion,
    pub(crate) absolute_errors: Option<BTreeMap<i32, Float>>,
}
impl From<crate::LaurentExpansion> for PyLaurentExpansion {
    fn from(inner: crate::LaurentExpansion) -> Self {
        Self {
            inner,
            absolute_errors: None,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyLaurentExpansion {
    /// Absolute propagated uncertainty for projected targets; absent for standalone fits.
    #[getter]
    fn absolute_errors(&self) -> Option<BTreeMap<i32, PythonMultiPrecisionFloat>> {
        self.absolute_errors.as_ref().map(|errors| {
            errors
                .iter()
                .map(|(&k, v)| (k, PythonMultiPrecisionFloat(v.clone())))
                .collect()
        })
    }
    #[getter]
    fn evidence_kind(&self) -> &'static str {
        if self.absolute_errors.is_some() {
            "propagated_boundary"
        } else if self.inner.verified_digits.is_some() {
            "independent_fits"
        } else {
            "unverified_fit"
        }
    }
    #[getter]
    fn coefficients(&self) -> BTreeMap<i32, PythonMultiPrecisionComplex> {
        self.inner
            .coefficients
            .iter()
            .map(|(&k, v)| (k, PythonMultiPrecisionComplex(v.clone())))
            .collect()
    }
    #[getter]
    fn verified_digits(&self) -> Option<u32> {
        self.inner.verified_digits
    }
    #[getter]
    fn working_bits(&self) -> u32 {
        self.inner.working_bits
    }
    #[getter]
    fn samples(&self) -> usize {
        self.inner.samples
    }
    #[getter]
    fn validation_samples(&self) -> usize {
        self.inner.validation_samples
    }
    #[getter]
    fn refinements(&self) -> usize {
        self.inner.refinements
    }
    #[getter]
    fn comparison_errors(&self) -> BTreeMap<i32, PythonMultiPrecisionFloat> {
        self.inner
            .comparison_errors
            .iter()
            .map(|(&k, v)| (k, PythonMultiPrecisionFloat(v.clone())))
            .collect()
    }
    /// Polynomial interpolation alone does not certify the fitted coefficients.
    #[staticmethod]
    #[pyo3(signature = (samples, values, leading, last, *, working_digits=80))]
    fn fit(
        py: Python<'_>,
        samples: Vec<PythonExpression>,
        values: Vec<PythonMultiPrecisionComplex>,
        leading: i32,
        last: i32,
        working_digits: u32,
    ) -> PyResult<Self> {
        let samples = samples.iter().map(rational).collect::<PyResult<Vec<_>>>()?;
        let values = values.into_iter().map(|v| v.0).collect::<Vec<_>>();
        let p = crate::Precision::decimal(working_digits).map_err(error)?;
        py.detach(|| crate::fit_epsilon(&samples, &values, leading, last, p))
            .map(Self::from)
            .map_err(error)
    }
}

/// Supplied regular-point initial values. Accuracy is an explicit caller claim.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "BoundaryData",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyBoundaryData {
    pub(crate) inner: crate::BoundaryData,
    pub(crate) verified_digits: Option<u32>,
    pub(crate) provenance: String,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBoundaryData {
    #[new]
    #[pyo3(signature = (point, values, *, verified_digits=None, provenance="caller supplied"))]
    fn new(
        point: PythonMultiPrecisionComplex,
        values: Vec<PythonMultiPrecisionComplex>,
        verified_digits: Option<u32>,
        provenance: &str,
    ) -> PyResult<Self> {
        let inner = crate::BoundaryData {
            point: point.0,
            values: values.into_iter().map(|v| v.0).collect(),
        };
        if inner.values.is_empty()
            || provenance.trim().is_empty()
            || verified_digits == Some(0)
            || std::iter::once(&inner.point)
                .chain(&inner.values)
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(InvalidInputError::new_err(
                "finite nonempty boundary and nonempty provenance required",
            ));
        }
        if let Some(digits) = verified_digits {
            let bits = crate::Precision::decimal(digits).map_err(error)?.bits;
            if inner
                .values
                .iter()
                .any(|v| v.re.prec() < bits || v.im.prec() < bits)
            {
                return Err(InvalidInputError::new_err(
                    "claimed accuracy exceeds boundary storage precision",
                ));
            }
        }
        Ok(Self {
            inner,
            verified_digits,
            provenance: provenance.to_owned(),
        })
    }
    #[getter]
    fn point(&self) -> PythonMultiPrecisionComplex {
        PythonMultiPrecisionComplex(self.inner.point.clone())
    }
    #[getter]
    fn values(&self) -> Vec<PythonMultiPrecisionComplex> {
        self.inner
            .values
            .iter()
            .cloned()
            .map(PythonMultiPrecisionComplex)
            .collect()
    }
    #[getter]
    fn verified_digits(&self) -> Option<u32> {
        self.verified_digits
    }
    #[getter]
    fn provenance(&self) -> &str {
        &self.provenance
    }
}

/// Regular-point rational differential equation dY/dx = A(x)Y.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DifferentialSystem",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDifferentialSystem {
    pub(crate) inner: crate::DifferentialSystem,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDifferentialSystem {
    #[new]
    fn new(variable: PythonExpression, matrix: Vec<Vec<PythonExpression>>) -> PyResult<Self> {
        let inner = crate::DifferentialSystem {
            variable: symbol(&variable)?,
            matrix: matrix
                .into_iter()
                .map(|row| row.into_iter().map(|x| x.expr).collect())
                .collect(),
        };
        inner.validate().map_err(error)?;
        Ok(Self { inner })
    }
    #[getter]
    fn matrix(&self) -> Vec<Vec<PythonExpression>> {
        self.inner
            .matrix
            .iter()
            .map(|r| r.iter().cloned().map(Into::into).collect())
            .collect()
    }
    #[getter]
    fn variable(&self) -> PythonExpression {
        Atom::var(self.inner.variable).into()
    }
    fn blocks(&self, py: Python<'_>) -> PyResult<Vec<Vec<usize>>> {
        py.detach(|| self.inner.blocks()).map_err(error)
    }
    fn change_basis(&self, py: Python<'_>, matrix: Vec<Vec<PythonExpression>>) -> PyResult<Self> {
        let matrix = matrix
            .into_iter()
            .map(|r| r.into_iter().map(|v| v.expr).collect())
            .collect::<Vec<_>>();
        py.detach(|| self.inner.change_basis(&matrix))
            .map(|inner| Self { inner })
            .map_err(error)
    }
    /// Fixed-precision continuation. Returned diagnostics do not certify a global error bound.
    #[pyo3(signature = (boundary, waypoints, *, parameters=None, options=None, control=None))]
    fn transport(
        &self,
        py: Python<'_>,
        boundary: &PyBoundaryData,
        waypoints: Vec<PythonMultiPrecisionComplex>,
        parameters: Option<NumericValues>,
        options: Option<&PyEvaluationOptions>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyTransportResult> {
        let options = super::options(options);
        options.validate().map_err(error)?;
        let p = crate::Precision::decimal(options.digits + options.guard_digits).map_err(error)?;
        let parameters = parameters
            .unwrap_or_default()
            .into_iter()
            .map(|(k, v)| (k.expr, v.0))
            .collect();
        let waypoints = waypoints.into_iter().map(|v| v.0).collect::<Vec<_>>();
        let context = context(control);
        let started = std::time::Instant::now();
        let result = py
            .detach(|| {
                self.inner.compile(p, &parameters)?.transport(
                    &boundary.inner,
                    &waypoints,
                    &options,
                    &context,
                )
            })
            .map_err(error)?;
        Ok(PyTransportResult::ordinary(result, boundary).timed(started))
    }
}

/// Shared, progressively filled cache of verified physical boundary values.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "BoundaryCache",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyBoundaryCache {
    pub(crate) inner: Arc<Mutex<crate::RustFlowCache>>,
}
impl PyBoundaryCache {
    pub(crate) fn access<T>(
        &self,
        f: impl FnOnce(&mut crate::RustFlowCache) -> crate::Result<T>,
    ) -> crate::Result<T> {
        let mut cache = self
            .inner
            .lock()
            .map_err(|_| crate::Error::Cache("boundary cache lock poisoned".into()))?;
        f(&mut cache)
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBoundaryCache {
    #[new]
    fn new() -> Self {
        Self::default()
    }
    #[staticmethod]
    fn load(py: Python<'_>, directory: PathBuf) -> PyResult<Self> {
        py.detach(|| crate::RustFlowCache::load(&directory))
            .map(|inner| Self {
                inner: Arc::new(Mutex::new(inner)),
            })
            .map_err(error)
    }
    /// Atomically persist versioned native binary entries and their provenance.
    fn save(&self, py: Python<'_>, directory: PathBuf) -> PyResult<()> {
        py.detach(|| self.access(|cache| cache.save(&directory)))
            .map_err(error)
    }
    /// Merge compatible entries while retaining intermediate points and stronger evidence.
    /// Repeated merges and merging a cache with itself are safe and idempotent.
    fn extend(&self, py: Python<'_>, other: &Self) -> PyResult<()> {
        py.detach(|| {
            // Release the source lock before acquiring the destination. Both
            // Python objects may refer to the same cache, or be merged from
            // different Python threads in opposite directions.
            let (entries, endpoints) = other.access(|cache| {
                Ok((cache.entries().to_vec(), cache.endpoint_entries().to_vec()))
            })?;
            self.access(|cache| cache.insert_batch(entries, endpoints))
        })
        .map_err(error)
    }
    /// Number of regular initial conditions; endpoint_entries() is separate.
    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        py.detach(|| self.access(|cache| Ok(cache.len())))
            .map_err(error)
    }
    /// Inspect retained points without exposing mutable cache entries.
    fn entries(&self, py: Python<'_>) -> PyResult<Vec<PyTransportResult>> {
        py.detach(|| {
            self.access(|cache| {
                cache
                    .entries()
                    .iter()
                    .cloned()
                    .map(PyTransportResult::cached)
                    .collect()
            })
        })
        .map_err(error)
    }
    /// Terminal endpoint records are separate from regular initial conditions.
    fn endpoint_entries(&self, py: Python<'_>) -> PyResult<Vec<PyEndpointResult>> {
        py.detach(|| {
            self.access(|cache| {
                cache
                    .endpoint_entries()
                    .iter()
                    .cloned()
                    .map(PyEndpointResult::cached)
                    .collect()
            })
        })
        .map_err(error)
    }
}

/// Numerical values and evidence for an ordinary or physical transport result.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "TransportResult",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyTransportResult {
    pub(crate) boundary: Option<crate::transport_cache::CachedBoundary>,
    pub(crate) point: Option<crate::ComplexFloat>,
    pub(crate) coordinates: BTreeMap<Symbol, Atom>,
    pub(crate) starting_coordinates: BTreeMap<Symbol, Atom>,
    pub(crate) coefficients: Vec<Vec<crate::ComplexFloat>>,
    pub(crate) leading: i32,
    pub(crate) errors: Vec<Vec<Float>>,
    pub(crate) verified_digits: Option<u32>,
    pub(crate) input_verified_digits: Option<u32>,
    pub(crate) diagnostics: crate::FlowDiagnostics,
    pub(crate) inserted_points: usize,
    pub(crate) cache_hit: bool,
    pub(crate) identity: Option<String>,
    ordinary_provenance: Option<String>,
    elapsed_nanoseconds: u128,
    attempts: Vec<crate::physical_transport::BoundaryAttempt>,
}
impl PyTransportResult {
    pub(crate) fn timed(mut self, started: std::time::Instant) -> Self {
        self.elapsed_nanoseconds = started.elapsed().as_nanos();
        self
    }
    fn ordinary(result: crate::FlowResult, boundary: &PyBoundaryData) -> Self {
        Self {
            boundary: None,
            point: Some(result.point),
            coordinates: BTreeMap::new(),
            starting_coordinates: BTreeMap::new(),
            coefficients: vec![result.values],
            leading: 0,
            errors: vec![],
            verified_digits: None,
            input_verified_digits: boundary.verified_digits,
            diagnostics: result.diagnostics,
            inserted_points: 0,
            cache_hit: false,
            identity: None,
            ordinary_provenance: Some(boundary.provenance.clone()),
            elapsed_nanoseconds: 0,
            attempts: vec![],
        }
    }
    pub(crate) fn cached(boundary: crate::transport_cache::CachedBoundary) -> crate::Result<Self> {
        Ok(Self {
            boundary: Some(boundary.clone()),
            point: None,
            coordinates: boundary.point.restart_coordinates()?,
            starting_coordinates: BTreeMap::new(),
            coefficients: boundary.coefficients,
            leading: boundary.range.leading,
            errors: boundary.accuracy.comparison_errors().to_vec(),
            verified_digits: Some(boundary.accuracy.verified_digits()),
            input_verified_digits: Some(boundary.accuracy.input_verified_digits()),
            diagnostics: crate::FlowDiagnostics {
                working_bits: boundary.accuracy.working_bits(),
                ..Default::default()
            },
            inserted_points: 0,
            cache_hit: true,
            identity: Some(boundary.identity.key().to_owned()),
            ordinary_provenance: None,
            elapsed_nanoseconds: 0,
            attempts: vec![],
        })
    }
    pub(crate) fn physical(
        result: crate::physical_transport::PhysicalResult,
    ) -> crate::Result<Self> {
        let start = result.starting_point.restart_coordinates()?;
        let mut out = Self::cached(result.boundary)?;
        out.starting_coordinates = start;
        out.inserted_points = result.inserted_points;
        out.cache_hit = result.transport.is_none();
        if let Some(attempt) = result.boundary_attempts.last() {
            out.input_verified_digits = Some(attempt.source_verified_digits);
        }
        out.attempts = result.boundary_attempts;
        if let Some(transport) = result.transport {
            out.diagnostics = transport.diagnostics;
        }
        Ok(out)
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyTransportResult {
    #[getter]
    fn provenance(&self) -> Option<&str> {
        self.boundary
            .as_ref()
            .map(|boundary| boundary.accuracy.provenance())
            .or(self.ordinary_provenance.as_deref())
    }
    /// Discrete root signs relative to the principal root at this exact point.
    #[getter]
    fn root_sheets(&self) -> HashMap<PythonExpression, i8> {
        self.boundary
            .as_ref()
            .and_then(|boundary| boundary.point.root_germ())
            .map(|germ| {
                germ.sheets
                    .iter()
                    .map(|(&s, sign)| {
                        (
                            Atom::var(s).into(),
                            match sign {
                                crate::transport_cache::RootSheet::Principal => 1,
                                crate::transport_cache::RootSheet::Opposite => -1,
                            },
                        )
                    })
                    .collect()
            })
            .unwrap_or_default()
    }
    /// Wall time in integer nanoseconds, including cache-lock wait; zero for cache inspection.
    #[getter]
    fn elapsed_nanoseconds(&self) -> u128 {
        self.elapsed_nanoseconds
    }
    /// Compatible sources actually attempted, in selection order.
    #[getter]
    fn attempted_starting_points(&self) -> PyResult<Vec<Expressions>> {
        self.attempts
            .iter()
            .map(|attempt| {
                Ok(attempt
                    .starting_point
                    .restart_coordinates()
                    .map_err(error)?
                    .into_iter()
                    .map(|(s, a)| (Atom::var(s).into(), a.into()))
                    .collect())
            })
            .collect()
    }
    #[getter]
    fn source_attempt_outcomes(&self) -> Vec<String> {
        self.attempts
            .iter()
            .map(|a| match &a.outcome {
                crate::physical_transport::BoundaryAttemptOutcome::Accepted => {
                    "accepted".to_owned()
                }
                crate::physical_transport::BoundaryAttemptOutcome::AccuracyRejected { message } => {
                    format!("accuracy rejected: {message}")
                }
            })
            .collect()
    }
    #[getter]
    fn source_attempt_costs(&self) -> Vec<PythonMultiPrecisionFloat> {
        self.attempts
            .iter()
            .map(|a| PythonMultiPrecisionFloat(a.cost.clone()))
            .collect()
    }
    #[getter]
    fn point(&self) -> Option<PythonMultiPrecisionComplex> {
        self.point.clone().map(PythonMultiPrecisionComplex)
    }
    #[getter]
    fn coordinates(&self) -> Expressions {
        self.coordinates
            .iter()
            .map(|(&k, v)| (Atom::var(k).into(), v.clone().into()))
            .collect()
    }
    #[getter]
    fn starting_coordinates(&self) -> Expressions {
        self.starting_coordinates
            .iter()
            .map(|(&k, v)| (Atom::var(k).into(), v.clone().into()))
            .collect()
    }
    /// Rows are ordered epsilon powers; columns are basis components.
    #[getter]
    fn coefficients(&self) -> Vec<Vec<PythonMultiPrecisionComplex>> {
        self.coefficients
            .iter()
            .map(|r| r.iter().cloned().map(PythonMultiPrecisionComplex).collect())
            .collect()
    }
    #[getter]
    fn leading_power(&self) -> i32 {
        self.leading
    }
    #[getter]
    fn comparison_errors(&self) -> Vec<Vec<PythonMultiPrecisionFloat>> {
        self.errors
            .iter()
            .map(|r| r.iter().cloned().map(PythonMultiPrecisionFloat).collect())
            .collect()
    }
    #[getter]
    fn verified_digits(&self) -> Option<u32> {
        self.verified_digits
    }
    #[getter]
    fn input_verified_digits(&self) -> Option<u32> {
        self.input_verified_digits
    }
    #[getter]
    fn working_bits(&self) -> u32 {
        self.diagnostics.working_bits
    }
    #[getter]
    fn steps(&self) -> usize {
        self.diagnostics.steps
    }
    #[getter]
    fn rejected_steps(&self) -> usize {
        self.diagnostics.rejected_steps
    }
    #[getter]
    fn predicate_evaluations(&self) -> usize {
        self.diagnostics.predicate_evaluations
    }
    #[getter]
    fn rational_trials(&self) -> usize {
        self.diagnostics.pade_trials
    }
    #[getter]
    fn rational_steps(&self) -> usize {
        self.diagnostics.pade_steps
    }
    #[getter]
    fn rational_fallbacks(&self) -> usize {
        self.diagnostics.pade_fallbacks
    }
    #[getter]
    fn last_rational_fallback(&self) -> Option<String> {
        self.diagnostics.last_pade_fallback.clone()
    }
    #[getter]
    fn conditioning_digits(&self) -> Option<u32> {
        self.diagnostics.conditioning_digits
    }
    #[getter]
    fn inserted_points(&self) -> usize {
        self.inserted_points
    }
    #[getter]
    fn cache_hit(&self) -> bool {
        self.cache_hit
    }
    #[getter]
    fn identity(&self) -> Option<String> {
        self.identity.clone()
    }
}
