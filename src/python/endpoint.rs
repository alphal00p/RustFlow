use super::*;
use crate::singular_endpoint::{EndpointBoundary, EndpointChart, EndpointResult};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Exact endpoint chart, with endpoint at parameter zero and a regular matching point.
/// Root signs and log winding describe the caller's explicitly admitted approach.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EndpointRoute",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEndpointRoute {
    pub(crate) chart: EndpointChart,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEndpointRoute {
    #[new]
    #[pyo3(signature=(parameter, coordinates, matching_parameter, *, homotopy, root_sheets=None, winding=0))]
    fn new(
        parameter: PythonExpression,
        coordinates: Expressions,
        matching_parameter: PythonExpression,
        homotopy: &str,
        root_sheets: Option<HashMap<PythonExpression, i8>>,
        winding: i32,
    ) -> PyResult<Self> {
        let chart = EndpointChart {
            path: crate::kinematics::KinematicPath {
                parameter: symbol(&parameter)?,
                coordinates: super::coordinates(coordinates)?,
            },
            matching_parameter: matching_parameter.expr,
            matching_germ: super::transport::germ(root_sheets)?,
            winding,
            homotopy: homotopy.into(),
        };
        chart.matching_point().map_err(error)?;
        Ok(Self { chart })
    }
    #[getter]
    fn parameter(&self) -> PythonExpression {
        Atom::var(self.chart.path.parameter).into()
    }
    #[getter]
    fn coordinates(&self) -> Expressions {
        self.chart
            .path
            .coordinates
            .iter()
            .map(|(&s, a)| (Atom::var(s).into(), a.clone().into()))
            .collect()
    }
    #[getter]
    fn matching_parameter(&self) -> PythonExpression {
        self.chart.matching_parameter.clone().into()
    }
    #[getter]
    fn matching_coordinates(&self) -> PyResult<Expressions> {
        Ok(self
            .chart
            .matching_point()
            .and_then(|p| p.restart_coordinates())
            .map_err(error)?
            .into_iter()
            .map(|(s, a)| (Atom::var(s).into(), a.into()))
            .collect())
    }
    #[getter]
    fn root_sheets(&self) -> HashMap<PythonExpression, i8> {
        self.chart
            .matching_germ
            .iter()
            .flat_map(|g| g.sheets.iter())
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
    }
    #[getter]
    fn winding(&self) -> i32 {
        self.chart.winding
    }
    #[getter]
    fn homotopy(&self) -> &str {
        &self.chart.homotopy
    }
}

/// Finite epsilon coefficient limits and estimated accuracy, with their regular anchor.
/// These terminal values cannot be used as initial data for a singular ODE.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EndpointResult",
    module = "symbolica.community.hep.integration",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEndpointResult {
    boundary: EndpointBoundary,
    matching: PyTransportResult,
    cache_hit: bool,
    elapsed_nanoseconds: u128,
}
impl PyEndpointResult {
    pub(crate) fn from_result(result: EndpointResult) -> crate::Result<Self> {
        let starting_point = result
            .boundary_attempts
            .last()
            .map(|a| a.starting_point.clone())
            .unwrap_or_else(|| result.boundary.matching_boundary.point.clone());
        let matching = PyTransportResult::physical(crate::physical_transport::PhysicalResult {
            boundary: result.boundary.matching_boundary.clone(),
            starting_point,
            transport: result.matching_transport,
            inserted_points: result.inserted_regular_points,
            boundary_attempts: result.boundary_attempts,
        })?;
        Ok(Self {
            boundary: result.boundary,
            matching,
            cache_hit: result.cache_hit,
            elapsed_nanoseconds: 0,
        })
    }
    pub(crate) fn cached(boundary: EndpointBoundary) -> crate::Result<Self> {
        Self::from_result(EndpointResult {
            boundary,
            matching_transport: None,
            inserted_regular_points: 0,
            cache_hit: true,
            boundary_attempts: Vec::new(),
        })
    }
    pub(crate) fn timed(mut self, started: std::time::Instant) -> Self {
        self.elapsed_nanoseconds = started.elapsed().as_nanos();
        self
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEndpointResult {
    #[getter]
    fn route(&self) -> PyEndpointRoute {
        PyEndpointRoute {
            chart: self.boundary.chart.clone(),
        }
    }
    #[getter]
    fn matching_boundary(&self) -> PyTransportResult {
        self.matching.clone()
    }
    /// Rows are epsilon powers, columns are the original physical basis.
    #[getter]
    fn coefficients(&self) -> Vec<Vec<PythonMultiPrecisionComplex>> {
        self.boundary
            .coefficients
            .iter()
            .map(|row| {
                row.iter()
                    .cloned()
                    .map(PythonMultiPrecisionComplex)
                    .collect()
            })
            .collect()
    }
    #[getter]
    fn comparison_errors(&self) -> Vec<Vec<PythonMultiPrecisionFloat>> {
        self.boundary
            .accuracy
            .comparison_errors()
            .iter()
            .map(|row| row.iter().cloned().map(PythonMultiPrecisionFloat).collect())
            .collect()
    }
    #[getter]
    fn leading_power(&self) -> i32 {
        self.boundary.range.leading
    }
    #[getter]
    fn last_power(&self) -> i32 {
        self.boundary.range.last
    }
    #[getter]
    fn verified_digits(&self) -> u32 {
        self.boundary.accuracy.verified_digits()
    }
    #[getter]
    fn input_verified_digits(&self) -> u32 {
        self.boundary.accuracy.input_verified_digits()
    }
    #[getter]
    fn working_bits(&self) -> u32 {
        self.boundary.accuracy.working_bits()
    }
    #[getter]
    fn provenance(&self) -> &str {
        self.boundary.accuracy.provenance()
    }
    #[getter]
    fn cache_hit(&self) -> bool {
        self.cache_hit
    }
    #[getter]
    fn elapsed_nanoseconds(&self) -> u128 {
        self.elapsed_nanoseconds
    }
}
