#![cfg(feature = "python")]
use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};
use symbolica::api::python::{PythonExpression, create_symbolica_module};
use symbolica::domains::float::{PythonMultiPrecisionComplex, PythonMultiPrecisionFloat};
use symbolica::prelude::*;
use symbolica_amflow::Precision;

#[test]
fn python_boundary_error_policy_runs_native_transport_and_persistent_reuse() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let locals = PyDict::new(py);
        locals.set_item("integration", integration)?;
        for key in ["x", "eps", "I1", "I2"] {
            locals.set_item(
                key,
                PythonExpression::from(
                    Atom::parse(key, "python_fundamental", Default::default()).unwrap(),
                ),
            )?;
        }
        for (key, value) in [("zero", 0), ("one", 1), ("rate", 5800)] {
            locals.set_item(key, PythonExpression::from(Atom::num(value)))?;
        }
        let p = Precision::decimal(90).unwrap();
        locals.set_item("n0", PythonMultiPrecisionComplex(p.zero()))?;
        locals.set_item("n1", PythonMultiPrecisionComplex(p.i(1)))?;
        locals.set_item("error", PythonMultiPrecisionFloat(p.tolerance(35)))?;
        py.run(
            c"
import tempfile
assert integration.EvaluationOptions().boundary_error_strategy == 'automatic'
assert integration.EvaluationOptions(boundary_error_strategy='fundamental_matrix').boundary_error_strategy == 'fundamental_matrix'
assert integration.EvaluationOptions(boundary_error_strategy='scalar_norm').boundary_error_strategy == 'scalar_norm'
opts = integration.EvaluationOptions()
try:
    integration.EvaluationOptions(boundary_error_strategy='unchecked')
    raise AssertionError('invalid strategy accepted')
except integration.InvalidInputError:
    pass
flow = integration.KinematicTransport(eps, {x:[[rate,rate],[-rate,-rate]]}, [I1,I2], one,
    branch_domain='real nilpotent transport with finite input evidence', options=opts)
bank = integration.BoundaryCache()
flow.add_boundary(bank, {x:zero}, [[n1,n0]], 0, verified_digits=30,
    comparison_errors=[[error,error]], provenance='finite 30-digit analytic source')
answer = flow.evaluate(bank, {x:one}, 0, 0, admit_straight_path=True)
assert answer.fundamental_boundary_charts > 0
assert answer.fundamental_boundary_fallback is None
assert answer.fundamental_boundary_retry is not None
assert 20 <= answer.verified_digits <= answer.input_verified_digits == 30
assert not isinstance(answer.coefficients[0][0], (float,complex))
with tempfile.TemporaryDirectory() as directory:
    bank.save(directory)
    resumed = integration.BoundaryCache.load(directory)
    hit = flow.evaluate(resumed, {x:one}, 0, 0)
    assert hit.cache_hit and hit.coefficients == answer.coefficients
    assert hit.comparison_errors == answer.comparison_errors
",
            Some(&locals),
            Some(&locals),
        )?;
        Ok(())
    })
    .unwrap();
}

#[cfg(feature = "python_stubgen")]
#[test]
fn boundary_error_policy_and_evidence_diagnostics_have_generated_stubs() {
    Python::initialize();
    let info = pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica.community.hep.integration".into(),
        env!("CARGO_MANIFEST_DIR").into(),
    )
    .unwrap();
    let text = info
        .modules
        .get("symbolica.community.hep.integration")
        .unwrap()
        .to_string();
    assert!(
        text.contains("boundary_error_strategy: builtins.str = 'automatic'"),
        "{text}"
    );
    assert!(text.contains("def boundary_error_strategy(self) -> builtins.str:"));
    assert!(text.contains("def fundamental_boundary_charts(self) -> builtins.int:"));
    assert!(
        text.contains("def fundamental_boundary_fallback(self) -> typing.Optional[builtins.str]:")
    );
    assert!(
        text.contains("def fundamental_boundary_retry(self) -> typing.Optional[builtins.str]:")
    );
}
