#![cfg(feature = "python")]
#[path = "support/native_cut_graph.rs"]
mod native_cut_graph;
#[path = "support/native_mixed_cut_graph.rs"]
mod native_mixed_cut_graph;
use feynkit_py::PyFeynmanDiagram;
use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};
use symbolica::api::python::{PythonExpression, create_symbolica_module};
use symbolica::domains::float::PythonMultiPrecisionComplex;
use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn cut_bindings_share_native_objects_preserve_precision_and_cancel_without_the_gil() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "symbolica.community.hepkit")?;
        feynkit_py::initialize_feynkit(&hep)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let local = PyDict::new(py);
        local.set_item("integration", &integration)?;
        local.set_item("hep", &hep)?;
        local.set_item("diagram", PyFeynmanDiagram::from(native_cut_graph::diagram()))?;
        local.set_item("D", PythonExpression::from(native_cut_graph::kinematics().dimension().to_symbolic()))?;
        local.set_item("P", PythonExpression::from(feynkit_graph::symbols::external_momentum().call(0)))?;
        local.set_item("s", PythonExpression::from(Atom::num(25)))?;
        local.set_item("eps", PythonExpression::from(Atom::var(symbol!("native_cut_python::eps"))))?;
        local.set_item("one", PythonExpression::from(Atom::one()))?;
        local.set_item("zero", PythonExpression::from(Atom::zero()))?;
        local.set_item("sample", PythonExpression::from(Atom::num(Rational::from((1, 13)))))?;
        py.run(c"
import concurrent.futures
import time
kinematics = hep.Kinematics(D).with_scalar_product(P, P, s)
options = integration.EvaluationOptions(digits=30, guard_digits=40, workers=2)
evaluator = integration.IntegralEvaluator(options=options)
control = integration.ComputationControl()
selection = dict(cut_index=0, future_channel=[one], loop_prescriptions=['insensitive'])
values = evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {}, eps, [sample], control=control, **selection)
assert any('Sample' in event for event in control.poll())
result = evaluator.evaluate_cut_diagram(diagram, kinematics, {}, eps, **selection)
assert result.verified_digits == 30
assert result.evidence_kind == 'independent_fits'
assert result.validation_samples > 0
assert not isinstance(values[0], (float, complex))
try:
    evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {}, eps, [zero], **selection)
    raise AssertionError('zero epsilon was accepted')
except integration.InvalidInputError:
    pass
try:
    evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {}, eps, [sample],
        cut_index=0, future_channel=[one], loop_prescriptions=['wrong'])
    raise AssertionError('invalid prescription was accepted')
except integration.InvalidInputError:
    pass
try:
    evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {}, eps, [sample],
        cut_index=0, future_channel=[one], loop_prescriptions=['+i0'])
    raise AssertionError('invalid physical measure was accepted')
except integration.UnsupportedInputError:
    pass
zero_values = evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {}, eps, [sample],
    edge_powers={1: 0}, **selection)
assert zero_values[0] == 0
raised_values = evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {}, eps, [sample],
    edge_powers={1: 2}, **selection)
assert raised_values[0] != 0
control = integration.ComputationControl()
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
    future = pool.submit(evaluator.evaluate_cut_diagram_samples, diagram, kinematics, {}, eps,
        [sample] * 10000, control=control, **selection)
    deadline = time.monotonic() + 30
    while not control.poll():
        assert not future.done(), 'calculation completed before progress could be polled'
        assert time.monotonic() < deadline, 'GIL held during native calculation'
        time.sleep(0.001)
    control.cancel()
    try:
        future.result(timeout=30)
        raise AssertionError('cancelled calculation was accepted')
    except integration.CalculationCancelled:
        pass
assert control.cancelled
", Some(&local), Some(&local))?;
        let value = local.get_item("values")?.unwrap().get_item(0)?.extract::<PythonMultiPrecisionComplex>()?;
        let p = Precision::decimal(80).unwrap();
        let eps = Rational::from((1, 13));
        let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
        let base = p.div(&p.i(1), &p.scale(&pi, 8, 1));
        let expected = p.mul(&base, &p.mul(&p.pow(&p.scale(&pi, 4, 25), &p.rational(&eps)),
            &p.div(&p.gamma_real(&p.rational(&(Rational::one()-&eps)).re).unwrap(),
                &p.gamma_real(&p.rational(&(Rational::from(2)-Rational::from(2)*eps)).re).unwrap())));
        assert!(p.close(&value.0, &expected, 50));
        Ok(())
    }).unwrap();
}

#[test]
fn explicit_slot_options_roundtrip_and_evaluate_native_partial_cut_flow() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "symbolica.community.hepkit")?;
        feynkit_py::initialize_feynkit(&hep)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let local = PyDict::new(py);
        local.set_item("integration", &integration)?;
        local.set_item("hep", &hep)?;
        local.set_item("diagram", PyFeynmanDiagram::from(native_mixed_cut_graph::connected_diagram()))?;
        local.set_item("terminal", PyFeynmanDiagram::from(native_cut_graph::diagram()))?;
        for (name, value) in [
            ("D", native_cut_graph::kinematics().dimension().to_symbolic()),
            ("P", feynkit_graph::symbols::external_momentum().call(0)),
            ("s", Atom::num(25)),
            ("M", Atom::var(symbol!("UFO::M"))),
            ("eps", Atom::var(symbol!("native_cut_python_slots::eps"))),
            ("one", Atom::one()),
            ("sample", Atom::num(Rational::from((1, 13)))),
        ] {
            local.set_item(name, PythonExpression::from(value))?;
        }
        py.run(c"
kinematics = hep.Kinematics(D).with_scalar_product(P, P, s)
options = integration.EvaluationOptions(guard_digits=50, deformed_propagator_slots=[2])
assert options.mass_mode == 'explicit'
assert options.deformed_propagator_slots == [2]
assert 'deformed_propagator_slots=[2]' in repr(options)
evaluator = integration.IntegralEvaluator(options=options)
retained = evaluator.options
roundtrip = integration.EvaluationOptions(mass_mode=retained.mass_mode,
    deformed_propagator_slots=retained.deformed_propagator_slots, guard_digits=50)
assert roundtrip.mass_mode == options.mass_mode
assert roundtrip.deformed_propagator_slots == options.deformed_propagator_slots
for mode in ['automatic', 'all', 'mass', 'propagator', 'branch', 'loop']:
    named = integration.EvaluationOptions(mass_mode=mode)
    assert named.mass_mode == mode and named.deformed_propagator_slots is None
    try:
        integration.EvaluationOptions(mass_mode=mode, deformed_propagator_slots=[2])
        raise AssertionError('conflicting placement was silently overridden')
    except integration.InvalidInputError:
        pass
for invalid in [dict(mass_mode='explicit'), dict(deformed_propagator_slots=[]), dict(mass_mode='unknown')]:
    try:
        integration.EvaluationOptions(**invalid)
        raise AssertionError('invalid placement settings were accepted')
    except integration.InvalidInputError:
        pass
selection = dict(cut_index=0, future_channel=[one], loop_prescriptions=['insensitive', '+i0'])
partial = evaluator.evaluate_cut_diagram_samples(diagram, kinematics, {M: one}, eps, [sample], **selection)
all_lines = integration.IntegralEvaluator(options=integration.EvaluationOptions(guard_digits=50, mass_mode='all'))
full = all_lines.evaluate_cut_diagram_samples(diagram, kinematics, {M: one}, eps, [sample], **selection)
assert not isinstance(partial[0], (float, complex))
# Slots 0/1 are cuts; slots >=4 are ISPs or out of bounds in this four-pole family.
for slot in [0, 1, 4, 999]:
    invalid = integration.IntegralEvaluator(options=integration.EvaluationOptions(deformed_propagator_slots=[slot]))
    try:
        invalid.evaluate_cut_diagram_samples(diagram, kinematics, {M: one}, eps, [sample], **selection)
        raise AssertionError('cut or nonphysical slot was accepted')
    except integration.InvalidInputError:
        pass
# Terminal dispatch must not silently ignore an explicit cut selection either.
try:
    evaluator.evaluate_cut_diagram_samples(terminal, kinematics, {}, eps, [sample],
        cut_index=0, future_channel=[one], loop_prescriptions=['insensitive'])
    raise AssertionError('terminal ignored explicit nonphysical slot')
except integration.InvalidInputError:
    pass
", Some(&local), Some(&local))?;
        let partial = local.get_item("partial")?.unwrap().get_item(0)?.extract::<PythonMultiPrecisionComplex>()?;
        let full = local.get_item("full")?.unwrap().get_item(0)?.extract::<PythonMultiPrecisionComplex>()?;
        assert!(Precision::decimal(70).unwrap().close(&partial.0, &full.0, 20));
        Ok(())
    }).unwrap();
}

#[cfg(feature = "python_stubgen")]
#[test]
fn cut_bindings_have_descriptive_generated_stubs_with_native_types() {
    Python::initialize();
    // inventory includes both native owners in the same linked extension.
    let info = pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica.community.hep.integration".into(),
        env!("CARGO_MANIFEST_DIR").into(),
    )
    .unwrap();
    let source = info
        .modules
        .get("symbolica.community.hep.integration")
        .unwrap()
        .to_string();
    for method in ["evaluate_cut_diagram", "evaluate_cut_diagram_samples"] {
        assert!(source.contains(&format!("def {method}(")), "{source}");
    }
    assert!(source.contains("future_channel:"));
    assert!(source.contains("loop_prescriptions:"));
    assert!(source.contains("edge_powers:"));
    assert!(source.contains("deformed_propagator_slots:"));
    assert!(source.contains("def mass_mode("));
    assert!(source.contains("FeynmanDiagram"));
    assert!(source.contains("Kinematics"));
    assert!(!source.contains("class RustFlow"));
}
