use super::*;
use crate::reduction::{Reduction, ReductionBackend};
use feynkit_py::{PyFeynmanDiagram, PyIntegralFamily, PyKinematics};
use std::sync::Arc;

type IntegralProjection = Vec<(Vec<i16>, PythonExpression)>;
type SuppliedRules = Vec<(Vec<i16>, IntegralProjection)>;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

fn family(
    native: &PyIntegralFamily,
    epsilon: &PythonExpression,
    physical: Option<usize>,
    dimension: i64,
) -> PyResult<crate::IntegralFamily> {
    let native = native.as_family();
    crate::IntegralFamily::from_hepkit(
        native,
        physical.unwrap_or(native.denominators().len()),
        symbol(epsilon)?,
        dimension,
    )
    .map_err(error)
}

/// Exact supplied reductions bound to one ordered family and parameter domain.
/// Auxiliary and recursive families require their own explicitly admitted tables.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ReductionTables",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyReductionTables {
    inner: crate::ScopedTableBackend,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyReductionTables {
    #[new]
    #[pyo3(signature=(name="supplied reductions"))]
    fn new(name: &str) -> Self {
        Self {
            inner: crate::ScopedTableBackend::new(name),
        }
    }
    /// Return a new collection containing this table; existing scopes cannot be replaced.
    #[pyo3(signature=(family, epsilon, rules, residuals, *, nonzero_conditions=None, physical_propagators=None, dimension=4))]
    #[allow(clippy::too_many_arguments)]
    fn with_family(
        &self,
        py: Python<'_>,
        family: &PyIntegralFamily,
        epsilon: PythonExpression,
        rules: SuppliedRules,
        residuals: Vec<Vec<i16>>,
        nonzero_conditions: Option<Vec<PythonExpression>>,
        physical_propagators: Option<usize>,
        dimension: i64,
    ) -> PyResult<Self> {
        let family = self::family(family, &epsilon, physical_propagators, dimension)?;
        let mut reduction = Reduction {
            residuals: residuals.into_iter().map(crate::Integral).collect(),
            nonzero_conditions: nonzero_conditions
                .unwrap_or_default()
                .into_iter()
                .map(|x| x.expr)
                .collect(),
            ..Default::default()
        };
        for (powers, terms) in rules {
            let mut combined = crate::reduction::LinearCombination::new();
            for (powers, coefficient) in terms {
                *combined.entry(crate::Integral(powers)).or_default() += coefficient.expr;
            }
            if reduction
                .rules
                .insert(crate::Integral(powers), combined)
                .is_some()
            {
                return Err(InvalidInputError::new_err("duplicate reduction target"));
            }
        }
        let mut inner = self.inner.clone();
        py.detach(|| inner.insert(&family, reduction))
            .map_err(error)?;
        Ok(Self { inner })
    }
}

/// Evaluate native HEPKit families with exact reduction and recursive boundary construction.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "IntegralEvaluator",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyIntegralEvaluator {
    pub(crate) options: crate::FlowOptions,
    pub(crate) backend: Arc<dyn ReductionBackend>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyIntegralEvaluator {
    /// Evaluate several exact linear combinations with shared sample values and reductions.
    /// Projection weights are multiplied before Laurent fitting.
    #[pyo3(signature=(family, projections, point, epsilon, *, last=0, physical_propagators=None, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate_projections(
        &self,
        py: Python<'_>,
        family: &PyIntegralFamily,
        projections: Vec<IntegralProjection>,
        point: Expressions,
        epsilon: PythonExpression,
        last: i32,
        physical_propagators: Option<usize>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<Vec<PyLaurentExpansion>> {
        let family = self::family(
            family,
            &epsilon,
            physical_propagators,
            self.options.dimension,
        )?;
        let rows = projections
            .into_iter()
            .map(|terms| {
                let mut row = crate::reduction::LinearCombination::new();
                for (powers, coefficient) in terms {
                    *row.entry(crate::Integral(powers)).or_default() += coefficient.expr;
                }
                row
            })
            .collect();
        let point = super::point(point);
        let context = context(control);
        py.detach(|| {
            crate::engine::solve_integral_projections(
                &[(family, rows)],
                &point,
                last,
                &self.options,
                self.backend.as_ref(),
                &context,
            )
        })
        .map(|values| values.into_iter().map(PyLaurentExpansion::from).collect())
        .map_err(error)
    }
    #[new]
    #[pyo3(signature=(*, options=None, reductions=None, reduction_depth=2, reduction_targets=4096, reduction_batch_size=32))]
    fn new(
        options: Option<&PyEvaluationOptions>,
        reductions: Option<&PyReductionTables>,
        reduction_depth: u32,
        reduction_targets: usize,
        reduction_batch_size: usize,
    ) -> PyResult<Self> {
        let options = super::options(options);
        options.validate().map_err(error)?;
        if reduction_targets == 0 || reduction_batch_size == 0 {
            return Err(InvalidInputError::new_err(
                "reduction_targets and reduction_batch_size must be positive",
            ));
        }
        let backend: Arc<dyn ReductionBackend> = match reductions {
            Some(tables) => Arc::new(tables.inner.clone()),
            None => Arc::new(crate::RustRedBackend {
                max_depth: reduction_depth,
                max_targets: reduction_targets,
                max_sector_batch: reduction_batch_size,
                ..Default::default()
            }),
        };
        Ok(Self { options, backend })
    }
    #[getter]
    fn options(&self) -> PyEvaluationOptions {
        PyEvaluationOptions {
            inner: self.options.clone(),
        }
    }
    #[getter]
    fn backend_identity(&self) -> String {
        self.backend.identity()
    }
    /// Denominator powers follow the native family's exact input order.
    #[pyo3(signature=(family, powers, point, epsilon, *, last=0, physical_propagators=None, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        py: Python<'_>,
        family: &PyIntegralFamily,
        powers: Vec<Vec<i16>>,
        point: Expressions,
        epsilon: PythonExpression,
        last: i32,
        physical_propagators: Option<usize>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<Vec<PyLaurentExpansion>> {
        let family = self::family(
            family,
            &epsilon,
            physical_propagators,
            self.options.dimension,
        )?;
        let targets = powers.into_iter().map(crate::Integral).collect::<Vec<_>>();
        let point = super::point(point);
        let context = context(control);
        py.detach(|| {
            crate::solve_integrals(
                &family,
                &targets,
                &point,
                last,
                &self.options,
                self.backend.as_ref(),
                &context,
            )
        })
        .map(|values| values.into_iter().map(PyLaurentExpansion::from).collect())
        .map_err(error)
    }
    /// Evaluate prescribed exact nonzero epsilon samples with fresh recursive boundaries.
    #[pyo3(signature=(family, powers, point, epsilon, samples, *, physical_propagators=None, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate_samples(
        &self,
        py: Python<'_>,
        family: &PyIntegralFamily,
        powers: Vec<Vec<i16>>,
        point: Expressions,
        epsilon: PythonExpression,
        samples: Vec<PythonExpression>,
        physical_propagators: Option<usize>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<Vec<Vec<PythonMultiPrecisionComplex>>> {
        let family = self::family(
            family,
            &epsilon,
            physical_propagators,
            self.options.dimension,
        )?;
        let targets = powers.into_iter().map(crate::Integral).collect::<Vec<_>>();
        let samples = samples.iter().map(rational).collect::<PyResult<Vec<_>>>()?;
        if samples.is_empty() || samples.iter().any(Rational::is_zero) {
            return Err(InvalidInputError::new_err(
                "at least one nonzero epsilon sample is required",
            ));
        }
        let point = super::point(point);
        let context = context(control);
        let values = py
            .detach(|| {
                let boundary = crate::recursive::RecursiveBoundary::new(
                    self.backend.as_ref(),
                    &self.options,
                    &context,
                );
                let evaluate = |(index, epsilon): (usize, &Rational)| {
                    context.emit(crate::Progress::Sample {
                        index,
                        total: samples.len(),
                    })?;
                    crate::PreparedFlow::new_at_epsilon(
                        &family,
                        &targets,
                        &point,
                        self.backend.as_ref(),
                        &self.options,
                        &context,
                        epsilon,
                    )?
                    .evaluate(epsilon, &self.options, &boundary, &context)
                };
                if self.options.workers == 1 {
                    samples
                        .iter()
                        .enumerate()
                        .map(evaluate)
                        .collect::<crate::Result<Vec<_>>>()
                } else {
                    use rayon::prelude::*;
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(self.options.workers)
                        .build()
                        .map_err(|e| crate::Error::InvalidInput(e.to_string()))?
                        .install(|| samples.par_iter().enumerate().map(evaluate).collect())
                }
            })
            .map_err(error)?;
        Ok(values
            .into_iter()
            .map(|row| row.into_iter().map(PythonMultiPrecisionComplex).collect())
            .collect())
    }
    /// Contract the native diagram numerator, decompose it, then integrate the exact weighted sum.
    #[pyo3(signature=(diagram, kinematics, point, epsilon, *, last=0, max_partial_fraction_states=10000, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate_diagram(
        &self,
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        kinematics: &PyKinematics,
        point: Expressions,
        epsilon: PythonExpression,
        last: i32,
        max_partial_fraction_states: usize,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyLaurentExpansion> {
        let diagram = diagram.as_shared_diagram()?.clone();
        let kinematics = kinematics.as_kinematics().clone();
        let point = super::point(point);
        let epsilon = symbol(&epsilon)?;
        let context = context(control);
        let inner = py
            .detach(|| {
                let graph = crate::hepkit::GraphIntegral::new(diagram, &kinematics)?;
                let groups = graph.integral_groups(
                    &point,
                    epsilon,
                    self.options.dimension,
                    max_partial_fraction_states,
                    &context,
                )?;
                crate::solve_integral_combinations(
                    &groups,
                    // integral_groups has already applied these substitutions
                    // simultaneously. Applying a -> b, b -> 3 twice changes a.
                    &crate::KinematicPoint::default(),
                    last,
                    &self.options,
                    self.backend.as_ref(),
                    &context,
                )
            })
            .map_err(error)?;
        Ok(PyLaurentExpansion::from(inner))
    }
    /// Prepare a common derivative-closed basis in all remaining physical variables.
    #[pyo3(signature=(family, powers, variables, epsilon, *, branch_domain, physical_propagators=None, fixed_parameters=None, epsilon_shearing=false, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        &self,
        py: Python<'_>,
        family: &PyIntegralFamily,
        powers: Vec<Vec<i16>>,
        variables: Vec<PythonExpression>,
        epsilon: PythonExpression,
        branch_domain: &str,
        physical_propagators: Option<usize>,
        fixed_parameters: Option<Expressions>,
        epsilon_shearing: bool,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyPreparedIntegralFamily> {
        let family = self::family(
            family,
            &epsilon,
            physical_propagators,
            self.options.dimension,
        )?
        .at(&point(fixed_parameters.unwrap_or_default()));
        let targets = powers.into_iter().map(crate::Integral).collect::<Vec<_>>();
        let variables = variables.iter().map(symbol).collect::<PyResult<Vec<_>>>()?;
        let context = context(control);
        let prepared = py
            .detach(|| {
                let prepared = crate::PreparedPhysicalFamily::new(
                    &family,
                    &targets,
                    &variables,
                    self.backend.as_ref(),
                    &self.options,
                    branch_domain,
                    &context,
                )?;
                if epsilon_shearing {
                    prepared.with_epsilon_shearing(&context)
                } else {
                    Ok(prepared)
                }
            })
            .map_err(error)?;
        Ok(PyPreparedIntegralFamily {
            inner: Arc::new(prepared),
            backend: self.backend.clone(),
            options: self.options.clone(),
        })
    }
}

/// A derivative-closed physical family that can generate its own integration boundaries.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "PreparedIntegralFamily",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyPreparedIntegralFamily {
    inner: Arc<crate::PreparedPhysicalFamily>,
    backend: Arc<dyn ReductionBackend>,
    options: crate::FlowOptions,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyPreparedIntegralFamily {
    /// Required master epsilon orders, accounting for poles in target reduction coefficients.
    fn required_master_range(
        &self,
        py: Python<'_>,
        point: Expressions,
        leading: i32,
        last: i32,
    ) -> PyResult<(i32, i32)> {
        let point = coordinates(point)?;
        let range = crate::transport_cache::EpsilonRange::new(leading, last).map_err(error)?;
        py.detach(|| self.inner.required_master_range(&point, range))
            .map(|range| (range.leading, range.last))
            .map_err(error)
    }
    /// Apply exact target reductions to a retained master result, propagating its uncertainty.
    #[pyo3(signature=(result, leading, last, *, digits=None))]
    fn project_targets(
        &self,
        py: Python<'_>,
        result: &PyTransportResult,
        leading: i32,
        last: i32,
        digits: Option<u32>,
    ) -> PyResult<Vec<PyLaurentExpansion>> {
        let boundary = result.boundary.as_ref().ok_or_else(|| {
            InvalidInputError::new_err("target projection requires a physical master result")
        })?;
        let range = crate::transport_cache::EpsilonRange::new(leading, last).map_err(error)?;
        let projected = py
            .detach(|| {
                self.inner
                    .project_targets(boundary, range, digits.unwrap_or(self.options.digits))
            })
            .map_err(error)?;
        Ok(projected
            .into_iter()
            .map(|v| PyLaurentExpansion {
                inner: crate::LaurentExpansion {
                    coefficients: v.coefficients,
                    verified_digits: Some(v.verified_digits),
                    working_bits: v.working_bits,
                    samples: 0,
                    validation_samples: 0,
                    refinements: 0,
                    comparison_errors: BTreeMap::new(),
                },
                absolute_errors: Some(v.absolute_errors),
            })
            .collect())
    }
    #[getter]
    fn basis(&self) -> Vec<Vec<i16>> {
        self.inner.basis().iter().map(|i| i.0.clone()).collect()
    }
    #[getter]
    fn nonzero_conditions(&self) -> Vec<PythonExpression> {
        self.inner
            .nonzero_conditions()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
    #[getter]
    fn target_reductions(&self) -> Vec<Vec<(Vec<i16>, PythonExpression)>> {
        self.inner
            .target_reductions()
            .iter()
            .map(|r| {
                r.iter()
                    .map(|(i, c)| (i.0.clone(), c.clone().into()))
                    .collect()
            })
            .collect()
    }
    /// Compute and verify fresh auxiliary-mass boundary values, then insert them atomically.
    #[pyo3(signature=(cache, point, *, last=0, extra_digits=10, control=None))]
    fn generate_boundary(
        &self,
        py: Python<'_>,
        cache: &PyBoundaryCache,
        point: Expressions,
        last: i32,
        extra_digits: u32,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyTransportResult> {
        let point = coordinates(point)?;
        let context = context(control);
        let started = std::time::Instant::now();
        py.detach(|| {
            cache.access(|cache| {
                let count = cache.len();
                let boundary = self.inner.seed_cache(
                    cache,
                    &point,
                    last,
                    extra_digits,
                    &self.options,
                    self.backend.as_ref(),
                    &context,
                )?;
                let mut result = PyTransportResult::cached(boundary)?;
                result.cache_hit = false;
                result.inserted_points = cache.len().saturating_sub(count);
                Ok(result.timed(started))
            })
        })
        .map_err(error)
    }
    /// Transport the common master basis; target reductions remain available separately.
    #[pyo3(signature=(cache, destination, leading, last, *, admit_straight_path=false, scales=None, control=None))]
    #[allow(clippy::too_many_arguments)]
    fn transport(
        &self,
        py: Python<'_>,
        cache: &PyBoundaryCache,
        destination: Expressions,
        leading: i32,
        last: i32,
        admit_straight_path: bool,
        scales: Option<Expressions>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyTransportResult> {
        let destination = coordinates(destination)?;
        let scales = coordinates(scales.unwrap_or_default())?;
        let context = context(control);
        let range = crate::transport_cache::EpsilonRange::new(leading, last).map_err(error)?;
        let precision = crate::Precision::decimal(self.options.digits + self.options.guard_digits)
            .map_err(error)?;
        let identity = self.inner.flow().identity();
        let policy = crate::transport_cache::ScaledDistance {
            scales,
            admissible: |source: &crate::transport_cache::CachedBoundary,
                         target: &crate::transport_cache::CachedPoint| {
                if source.point.restart_coordinates()? == target.restart_coordinates()? {
                    return Ok(true);
                }
                if !admit_straight_path {
                    return Ok(false);
                }
                identity.conditions_admit_straight_path(
                    &source.point,
                    target,
                    precision,
                    self.options.digits,
                )
            },
        };
        let started = std::time::Instant::now();
        py.detach(|| {
            cache
                .access(|cache| {
                    self.inner.flow().evaluate_to(
                        cache,
                        &destination,
                        range,
                        &self.options,
                        &context,
                        &policy,
                    )
                })
                .and_then(PyTransportResult::physical)
        })
        .map(|result| result.timed(started))
        .map_err(error)
    }
}
