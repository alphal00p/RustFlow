use super::*;
use crate::gg_hg::{HiggsJetConfiguration, HiggsJetIntegralSystem, PluginFamilyKind};
use crate::transport_cache::{RootGerm, RootSheet};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::sync::Arc;

fn coordinates(values: Expressions) -> PyResult<BTreeMap<Symbol, Atom>> {
    values
        .into_iter()
        .map(|(s, a)| Ok((symbol(&s)?, a.expr)))
        .collect()
}
fn sheets(values: HashMap<PythonExpression, i8>) -> PyResult<RootGerm> {
    Ok(RootGerm {
        sheets: values
            .into_iter()
            .map(|(s, sign)| {
                Ok((
                    symbol(&s)?,
                    match sign {
                        1 => RootSheet::Principal,
                        -1 => RootSheet::Opposite,
                        _ => return Err(InvalidInputError::new_err("root sheet must be +1 or -1")),
                    },
                ))
            })
            .collect::<PyResult<_>>()?,
    })
}

/// Exact physical coordinates and root conventions for one mass and crossing.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "HiggsJetConfiguration",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyHiggsJetConfiguration {
    inner: HiggsJetConfiguration,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyHiggsJetConfiguration {
    #[getter]
    fn label(&self) -> &str {
        &self.inner.label
    }
    #[getter]
    fn mass(&self) -> &str {
        &self.inner.mass
    }
    #[getter]
    fn permutation(&self) -> usize {
        self.inner.permutation
    }
    #[getter]
    fn start(&self) -> Expressions {
        self.inner
            .start
            .iter()
            .map(|(&s, a)| (Atom::var(s).into(), a.clone().into()))
            .collect()
    }
    #[getter]
    fn destination(&self) -> Expressions {
        self.inner
            .destination
            .iter()
            .map(|(&s, a)| (Atom::var(s).into(), a.clone().into()))
            .collect()
    }
    #[getter]
    fn root_sheets(&self) -> HashMap<PythonExpression, i8> {
        self.inner
            .germ
            .sheets
            .iter()
            .map(|(&s, g)| {
                (
                    Atom::var(s).into(),
                    if *g == RootSheet::Principal { 1 } else { -1 },
                )
            })
            .collect()
    }
}

/// Certified planar or nonplanar Higgs-plus-jet integral system. Generates its
/// own boundary values using the supplied IntegralEvaluator's reduction backend.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "HiggsJetIntegralSystem",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyHiggsJetIntegralSystem {
    inner: Arc<HiggsJetIntegralSystem>,
    namespace: String,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyHiggsJetIntegralSystem {
    #[new]
    #[pyo3(signature=(topology, *, namespace="hep_higgs_jet"))]
    fn new(py: Python<'_>, topology: &str, namespace: &str) -> PyResult<Self> {
        let kind = match topology {
            "planar" => PluginFamilyKind::Planar,
            "nonplanar" => PluginFamilyKind::Nonplanar,
            _ => {
                return Err(InvalidInputError::new_err(
                    "topology must be planar or nonplanar",
                ));
            }
        };
        let inner = py
            .detach(|| HiggsJetIntegralSystem::load(kind, namespace))
            .map_err(error)?;
        Ok(Self {
            inner: Arc::new(inner),
            namespace: namespace.into(),
        })
    }
    /// Whether this build contains native reduction and boundary generation.
    #[getter]
    fn automatic_boundary_generation_available(&self) -> bool {
        cfg!(feature = "automatic")
    }
    #[getter]
    fn dimension(&self) -> usize {
        self.inner.map.kind.dimension()
    }
    /// Stable mathematical certificate for portable, externally supplied seeds.
    /// This does not bypass runtime-sensitive binary cache compatibility checks.
    #[getter]
    fn mathematical_fingerprint(&self) -> String {
        self.inner.mathematical_fingerprint()
    }
    /// Share this system's canonical connection and ordinary supplied-boundary
    /// validation. No reduction, boundary generation or numeric transport runs.
    #[pyo3(signature=(options=None))]
    fn kinematic_transport(
        &self,
        options: Option<&PyEvaluationOptions>,
    ) -> PyResult<PyKinematicTransport> {
        let options = super::options(options);
        options.validate().map_err(error)?;
        if options.dimension != 4 || options.prescription != crate::Prescription::PlusI0 {
            return Err(error(crate::Error::Unsupported(
                "the certified Higgs-jet canonical basis uses D=4-2 epsilon and +i0".into(),
            )));
        }
        Ok(PyKinematicTransport::from_canonical(
            self.inner.transport.clone(),
            options,
        ))
    }
    #[getter]
    fn epsilon(&self) -> PythonExpression {
        Atom::var(self.inner.map.epsilon).into()
    }
    #[getter]
    fn coordinates(&self) -> Vec<PythonExpression> {
        self.inner
            .map
            .coordinates
            .iter()
            .map(|&s| Atom::var(s).into())
            .collect()
    }
    #[getter]
    fn roots(&self) -> Expressions {
        self.inner
            .map
            .roots
            .iter()
            .map(|r| (Atom::var(r.symbol).into(), r.radicand.clone().into()))
            .collect()
    }
    #[getter]
    fn physical_powers(&self) -> Vec<Vec<i16>> {
        self.inner
            .map
            .integrals
            .iter()
            .map(|i| i.0.clone())
            .collect()
    }
    #[getter]
    fn canonical_projection(&self) -> Vec<Vec<(Vec<i16>, PythonExpression)>> {
        self.inner
            .map
            .canonical_rows()
            .into_iter()
            .map(|r| r.into_iter().map(|(i, a)| (i.0, a.into())).collect())
            .collect()
    }
    #[getter]
    fn provenance(&self) -> String {
        format!(
            "{}; sha256={}; exact inverse and differential-matrix certificates",
            self.inner.map.paper.url, self.inner.map.paper.sha256
        )
    }
    fn configurations(&self, py: Python<'_>) -> PyResult<Vec<PyHiggsJetConfiguration>> {
        py.detach(|| self.inner.configurations(&self.namespace))
            .map(|values| {
                values
                    .into_iter()
                    .map(|inner| PyHiggsJetConfiguration { inner })
                    .collect()
            })
            .map_err(error)
    }
    fn verify_basis(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.inner.map.verify_inverse()).map_err(error)
    }

    /// Complete canonical seed generation from ordinary physical integrals.
    /// Recompute bypasses verified-boundary and numerical-sample reuse while
    /// retaining the evaluator's exact reduction cache.
    #[pyo3(signature=(evaluator,cache,point,root_sheets,*,last=4,recompute=false,control=None))]
    #[allow(clippy::too_many_arguments)]
    #[cfg(feature = "automatic")]
    fn generate_boundary(
        &self,
        py: Python<'_>,
        evaluator: &PyIntegralEvaluator,
        cache: &PyBoundaryCache,
        point: Expressions,
        root_sheets: HashMap<PythonExpression, i8>,
        last: i32,
        recompute: bool,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyTransportResult> {
        let point = coordinates(point)?;
        let germ = sheets(root_sheets)?;
        let context = context(control);
        let started = std::time::Instant::now();
        py.detach(|| {
            cache.access(|cache| {
                let had_boundary = !recompute
                    && cache.entries().iter().any(|b| {
                        b.identity.key() == self.inner.transport.identity().key()
                            && b.point.restart_coordinates().ok().as_ref() == Some(&point)
                            && b.point.root_germ() == Some(&germ)
                            && b.range.leading == 0
                            && b.range.last >= last
                            && b.accuracy.verified_digits() >= evaluator.options.digits
                    });
                let boundary = self.inner.generate_boundary(
                    cache,
                    &point,
                    &germ,
                    last,
                    &evaluator.options,
                    evaluator.backend.as_ref(),
                    &context,
                    !recompute,
                )?;
                let mut result = PyTransportResult::cached(boundary)?;
                result.cache_hit = had_boundary;
                Ok(result.timed(started))
            })
        })
        .map_err(error)
    }

    /// Continue within a regular real kinematic chamber, saving accepted
    /// intermediate physical points in the shared cache.
    #[pyo3(signature=(cache,point,root_sheets,*,last=4,options=None,control=None))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        py: Python<'_>,
        cache: &PyBoundaryCache,
        point: Expressions,
        root_sheets: HashMap<PythonExpression, i8>,
        last: i32,
        options: Option<&PyEvaluationOptions>,
        control: Option<&PyComputationControl>,
    ) -> PyResult<PyTransportResult> {
        let point = coordinates(point)?;
        let germ = sheets(root_sheets)?;
        let options = super::options(options);
        let context = context(control);
        let started = std::time::Instant::now();
        py.detach(|| {
            cache.access(|cache| {
                self.inner
                    .evaluate(cache, &point, &germ, last, &options, &context)
                    .and_then(PyTransportResult::physical)
            })
        })
        .map(|r| r.timed(started))
        .map_err(error)
    }
}

/// Kinematics-dependent projection of crossed master integrals into the four
/// normalized scalar Higgs-plus-jet form factors for one vector-boson mass.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "HiggsJetFormFactorProjector",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyHiggsJetFormFactorProjector {
    inner: Arc<crate::gg_hg::HiggsJetFormFactors>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyHiggsJetFormFactorProjector {
    #[new]
    #[pyo3(signature=(*,namespace="hep_higgs_jet"))]
    fn new(py: Python<'_>, namespace: &str) -> PyResult<Self> {
        py.detach(|| crate::gg_hg::HiggsJetFormFactors::load(namespace))
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(error)
    }
    /// Planar/nonplanar results are each ordered as (s,u), (s,t), (u,t), (t,s).
    #[pyo3(signature=(s,t,higgs_mass_squared,vector_mass_squared,planar,nonplanar,*,working_digits=100))]
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        py: Python<'_>,
        s: PythonExpression,
        t: PythonExpression,
        higgs_mass_squared: PythonExpression,
        vector_mass_squared: PythonExpression,
        planar: Vec<PyTransportResult>,
        nonplanar: Vec<PyTransportResult>,
        working_digits: u32,
    ) -> PyResult<PyFormFactorResult> {
        let boundary = |results: Vec<PyTransportResult>| {
            results
                .into_iter()
                .map(|r| {
                    r.boundary.ok_or_else(|| {
                        InvalidInputError::new_err(
                            "form-factor projections require canonical physical transport results",
                        )
                    })
                })
                .collect::<PyResult<Vec<_>>>()
        };
        let planar = boundary(planar)?;
        let nonplanar = boundary(nonplanar)?;
        let p = crate::Precision::decimal(working_digits).map_err(error)?;
        py.detach(|| {
            self.inner.evaluate(
                &s.expr,
                &t.expr,
                &higgs_mass_squared.expr,
                &vector_mass_squared.expr,
                &planar,
                &nonplanar,
                p,
            )
        })
        .map(|inner| PyFormFactorResult { inner })
        .map_err(error)
    }
}

/// Four normalized form factors, with propagated master-integral uncertainties.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "FormFactorResult",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyFormFactorResult {
    inner: crate::gg_hg::HiggsJetFormFactorResult,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyFormFactorResult {
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
    fn absolute_errors(&self) -> Vec<PythonMultiPrecisionFloat> {
        self.inner
            .absolute_errors
            .iter()
            .cloned()
            .map(PythonMultiPrecisionFloat)
            .collect()
    }
    #[getter]
    fn verified_relative_digits(&self) -> Vec<Option<u32>> {
        self.inner.verified_relative_digits.to_vec()
    }
    #[getter]
    fn provenance(&self) -> &str {
        &self.inner.provenance
    }
}
