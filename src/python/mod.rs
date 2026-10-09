//! Bindings registered into the host's existing Symbolica Python extension.
//!
//! This module deliberately has no `pymodule` entry point: loading a second
//! Symbolica runtime would invalidate native expression and family identities.
mod amplitude;
mod citations;
mod constraints;
mod continuation;
mod endpoint;
#[cfg(feature = "automatic")]
mod evaluation;
#[cfg(feature = "automatic")]
mod finite_density;
mod gg_hg;
mod transport;
mod types;

use pyo3::prelude::*;
use pyo3::types::PyModule;
use std::collections::{BTreeMap, HashMap};
use symbolica::api::python::PythonExpression;
use symbolica::domains::float::{PythonMultiPrecisionComplex, PythonMultiPrecisionFloat};
use symbolica::prelude::*;

pub use amplitude::{PyAmplitudeResult, PyHiggsJetAmplitude};
pub use citations::get_citations;
pub use constraints::{PyAsymptoticCoefficient, PyAsymptoticRelation, PyEndpointConstraints};
pub use continuation::PyContinuationPrescription;
pub use endpoint::{PyEndpointResult, PyEndpointRoute};
#[cfg(feature = "automatic")]
pub use evaluation::{PyIntegralEvaluator, PyPreparedIntegralFamily, PyReductionTables};
pub use gg_hg::{
    PyFormFactorResult, PyHiggsJetConfiguration, PyHiggsJetFormFactorProjector,
    PyHiggsJetIntegralSystem,
};
pub use transport::PyKinematicTransport;
pub use types::{
    PyBoundaryCache, PyBoundaryData, PyComputationControl, PyDifferentialSystem,
    PyEvaluationOptions, PyLaurentExpansion, PyTransportResult,
};

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::module_variable!(
    "symbolica.community.hep.integration",
    "automatic_boundary_generation_available",
    bool
);

type Expressions = HashMap<PythonExpression, PythonExpression>;
type NumericValues = HashMap<PythonExpression, PythonMultiPrecisionComplex>;

pyo3::create_exception!(
    symbolica.community.hep.integration,
    EvaluationError,
    pyo3::exceptions::PyException
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    InvalidInputError,
    EvaluationError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    UnsupportedInputError,
    EvaluationError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    ReductionError,
    EvaluationError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    IncompleteReductionError,
    ReductionError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    NumericalError,
    EvaluationError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    AccuracyError,
    NumericalError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    PrecisionError,
    AccuracyError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    ResourceLimitError,
    EvaluationError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    CalculationCancelled,
    EvaluationError
);
pyo3::create_exception!(
    symbolica.community.hep.integration,
    BoundaryCacheError,
    EvaluationError
);

#[cfg(feature = "python_stubgen")]
macro_rules! exception_stub {
    ($name:ident, $base:expr) => {
        pyo3_stub_gen::inventory::submit! {
            pyo3_stub_gen::type_info::PyClassInfo {
                pyclass_name: stringify!($name),
                struct_id: std::any::TypeId::of::<$name>,
                getters: &[], setters: &[],
                module: Some("symbolica.community.hep.integration"),
                doc: "Numerical integration failure; see the exception message for context.",
                bases: &[|| $base], has_eq: false, has_ord: false, has_hash: false,
                has_str: false, subclass: true,
            }
        }
    };
}

#[cfg(feature = "python_stubgen")]
mod exception_stubs {
    use super::*;
    use pyo3_stub_gen::TypeInfo;
    exception_stub!(EvaluationError, TypeInfo::builtin("Exception"));
    exception_stub!(InvalidInputError, TypeInfo::unqualified("EvaluationError"));
    exception_stub!(
        UnsupportedInputError,
        TypeInfo::unqualified("EvaluationError")
    );
    exception_stub!(ReductionError, TypeInfo::unqualified("EvaluationError"));
    exception_stub!(
        IncompleteReductionError,
        TypeInfo::unqualified("ReductionError")
    );
    exception_stub!(NumericalError, TypeInfo::unqualified("EvaluationError"));
    exception_stub!(AccuracyError, TypeInfo::unqualified("NumericalError"));
    exception_stub!(PrecisionError, TypeInfo::unqualified("AccuracyError"));
    exception_stub!(ResourceLimitError, TypeInfo::unqualified("EvaluationError"));
    exception_stub!(
        CalculationCancelled,
        TypeInfo::unqualified("EvaluationError")
    );
    exception_stub!(BoundaryCacheError, TypeInfo::unqualified("EvaluationError"));
}

fn error(e: crate::Error) -> PyErr {
    let message = e.to_string();
    match e {
        crate::Error::InvalidInput(_) => InvalidInputError::new_err(message),
        crate::Error::Unsupported(_) => UnsupportedInputError::new_err(message),
        crate::Error::Reduction(_) => ReductionError::new_err(message),
        crate::Error::IncompleteReduction(_) => IncompleteReductionError::new_err(message),
        crate::Error::Numerical(_) => NumericalError::new_err(message),
        crate::Error::Accuracy(_) => AccuracyError::new_err(message),
        crate::Error::InsufficientPrecision { .. } => PrecisionError::new_err(message),
        crate::Error::Limit(_) => ResourceLimitError::new_err(message),
        crate::Error::Cancelled => CalculationCancelled::new_err(message),
        crate::Error::Cache(_) => BoundaryCacheError::new_err(message),
        crate::Error::Io(e) => pyo3::exceptions::PyOSError::new_err(e.to_string()),
    }
}

fn symbol(expression: &PythonExpression) -> PyResult<Symbol> {
    match expression.expr.as_view() {
        AtomView::Var(v) => Ok(v.get_symbol()),
        _ => Err(InvalidInputError::new_err("expected a Symbolica variable")),
    }
}

fn coordinates(values: Expressions) -> PyResult<BTreeMap<Symbol, Atom>> {
    values
        .into_iter()
        .map(|(k, v)| Ok((symbol(&k)?, v.expr)))
        .collect()
}

fn point(values: Expressions) -> crate::KinematicPoint {
    crate::KinematicPoint(values.into_iter().map(|(k, v)| (k.expr, v.expr)).collect())
}

fn rational(value: &PythonExpression) -> PyResult<Rational> {
    if let AtomView::Num(n) = value.expr.as_view()
        && let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned()
        && c.im.is_zero()
    {
        return Ok(c.re);
    }
    Err(InvalidInputError::new_err(
        "epsilon samples must be exact rational Symbolica expressions",
    ))
}

fn context(control: Option<&PyComputationControl>) -> crate::RunContext {
    control
        .map(PyComputationControl::context)
        .unwrap_or_default()
}

fn options(value: Option<&PyEvaluationOptions>) -> crate::FlowOptions {
    value.map(|v| v.inner.clone()).unwrap_or_default()
}

/// Add integration classes to `symbolica.community.hep.integration`.
/// The community host owns module creation and the shared Symbolica runtime.
pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "automatic")]
    module.add_function(wrap_pyfunction!(
        finite_density::prepare_finite_density,
        module
    )?)?;
    // CPython's WASM object allocator guarantees eight-byte alignment, whereas
    // Rust's u128 (and any aggregate containing it) requires sixteen. Checking
    // every registered wrapper here also covers nested numerical payloads and
    // automatic-only classes, and runs in release as well as debug builds.
    // Keep exact durations as Duration; put any other over-aligned owned data
    // behind Box/Arc instead of narrowing it or relying on allocator luck.
    macro_rules! add_class {
        ($class:ty) => {{
            const _: () = assert!(
                std::mem::align_of::<$class>() <= 8,
                concat!(
                    stringify!($class),
                    " exceeds Python's WASM object alignment"
                ),
            );
            const _: () = assert!(
                std::mem::align_of::<<$class as pyo3::impl_::pyclass::PyClassImpl>::Layout>() <= 8,
                concat!(stringify!($class), " has an over-aligned Python layout"),
            );
            module.add_class::<$class>()?;
        }};
    }
    add_class!(PyHiggsJetIntegralSystem);
    add_class!(PyHiggsJetConfiguration);
    add_class!(PyHiggsJetFormFactorProjector);
    add_class!(PyFormFactorResult);
    add_class!(PyHiggsJetAmplitude);
    add_class!(PyAmplitudeResult);
    add_class!(PyEvaluationOptions);
    add_class!(PyComputationControl);
    #[cfg(feature = "automatic")]
    add_class!(PyIntegralEvaluator);
    #[cfg(feature = "automatic")]
    add_class!(PyPreparedIntegralFamily);
    #[cfg(feature = "automatic")]
    add_class!(PyReductionTables);
    add_class!(PyKinematicTransport);
    add_class!(PyContinuationPrescription);
    add_class!(PyEndpointRoute);
    add_class!(PyAsymptoticCoefficient);
    add_class!(PyAsymptoticRelation);
    add_class!(PyEndpointConstraints);
    add_class!(PyEndpointResult);
    add_class!(PyBoundaryCache);
    add_class!(PyBoundaryData);
    add_class!(PyDifferentialSystem);
    add_class!(PyLaurentExpansion);
    add_class!(PyTransportResult);
    macro_rules! add_error {
        ($($name:ident),* $(,)?) => { $(module.add(stringify!($name), module.py().get_type::<$name>())?;)* };
    }
    add_error!(
        EvaluationError,
        InvalidInputError,
        UnsupportedInputError,
        ReductionError,
        IncompleteReductionError,
        NumericalError,
        AccuracyError,
        PrecisionError,
        ResourceLimitError,
        CalculationCancelled,
        BoundaryCacheError
    );
    module.add(
        "__all__",
        [
            "automatic_boundary_generation_available",
            "HiggsJetIntegralSystem",
            "HiggsJetConfiguration",
            "HiggsJetFormFactorProjector",
            "FormFactorResult",
            "IntegralEvaluator",
            "HiggsJetAmplitude",
            "AmplitudeResult",
            "PreparedIntegralFamily",
            "KinematicTransport",
            "ContinuationPrescription",
            "EndpointRoute",
            "AsymptoticCoefficient",
            "AsymptoticRelation",
            "EndpointConstraints",
            "EndpointResult",
            "BoundaryCache",
            "EvaluationOptions",
            "DifferentialSystem",
            "BoundaryData",
            "LaurentExpansion",
            "TransportResult",
            "ComputationControl",
            "ReductionTables",
            "EvaluationError",
            "InvalidInputError",
            "UnsupportedInputError",
            "ReductionError",
            "IncompleteReductionError",
            "NumericalError",
            "AccuracyError",
            "PrecisionError",
            "ResourceLimitError",
            "CalculationCancelled",
            "BoundaryCacheError",
        ]
        .into_iter()
        .filter(|name| {
            cfg!(feature = "automatic")
                || !matches!(
                    *name,
                    "IntegralEvaluator" | "PreparedIntegralFamily" | "ReductionTables"
                )
        })
        .collect::<Vec<_>>(),
    )?;
    module.add(
        "automatic_boundary_generation_available",
        cfg!(feature = "automatic"),
    )?;
    Ok(())
}
