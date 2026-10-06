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
fn prescribed_python_routes_share_native_branches_evidence_and_terminal_cache() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let local = PyDict::new(py);
        local.set_item("integration", integration)?;
        let p = Precision::decimal(90).unwrap();
        for key in ["s", "eps", "r", "z", "Y"] {
            local.set_item(key, PythonExpression::from(Atom::parse(key, "python_prescribed", Default::default()).unwrap()))?;
        }
        for (key, value) in [("one", 1), ("zero", 0), ("minus_one", -1), ("minus_two", -2)] {
            local.set_item(key, PythonExpression::from(Atom::num(value)))?;
        }
        local.set_item("half", PythonExpression::from(Atom::num(Rational::from((1, 2)))))?;
        local.set_item("n1", PythonMultiPrecisionComplex(p.i(1)))?;
        local.set_item("nm1", PythonMultiPrecisionComplex(p.i(-1)))?;
        local.set_item("n0", PythonMultiPrecisionComplex(p.zero()))?;
        local.set_item("error0", PythonMultiPrecisionFloat(p.real(0)))?;
        py.run(c"
import tempfile
opts = integration.EvaluationOptions(digits=20, guard_digits=30, series_order=48)
def declaration(side):
    return integration.ContinuationPrescription([(s, side)], domain=side + ' physical threshold', unprescribed_side=side)
def make_flow(side, canonical=False):
    continuation = declaration(side)
    assert continuation.prescriptions == [(s, side)]
    assert continuation.unprescribed_side == side
    if canonical:
        return integration.KinematicTransport.canonical(eps, [s], [r], [[[one]]], [Y], one,
            roots={r:s}, branch_domain='selected logarithm', options=opts, continuation=continuation)
    return integration.KinematicTransport(eps, {s:[[eps/s]]}, [Y], one,
        branch_domain='selected logarithm', options=opts, continuation=continuation)
def seed(flow, bank, at=one, values=None, sheets=None):
    flow.add_boundary(bank, {s:at}, values or [[n1], [n0], [n0]], 0,
        verified_digits=60, comparison_errors=[[error0] for _ in (values or [0,0,0])],
        provenance='exact analytic normalization', root_sheets=sheets)

upper, lower = make_flow('+i0'), make_flow('-i0')
assert upper.identity != lower.identity
upper_bank, lower_bank = integration.BoundaryCache(), integration.BoundaryCache()
seed(upper, upper_bank)
seed(lower, lower_bank)
try:
    upper.evaluate(upper_bank, {s:minus_one}, 0, 2, admit_straight_path=True)
    raise AssertionError('straight-path admission authorized a prescribed detour')
except integration.EvaluationError:
    pass
assert len(upper_bank) == 1
up = upper.evaluate(upper_bank, {s:minus_one}, 0, 2, admit_prescribed_path=True)
down = lower.evaluate(lower_bank, {s:minus_one}, 0, 2, admit_prescribed_path=True)
assert 20 <= up.verified_digits <= up.input_verified_digits == 60
assert 20 <= down.verified_digits <= down.input_verified_digits == 60
assert up.input_verified_digits == down.input_verified_digits == 60
assert up.steps > 0 and up.inserted_points > 0
repeat = upper.evaluate(upper_bank, {s:minus_one}, 0, 2)
assert repeat.cache_hit and repeat.steps == 0 and repeat.coefficients == up.coefficients
assert not isinstance(up.coefficients[1][0], (float, complex))
count = len(upper_bank)
try:
    lower.evaluate(upper_bank, {s:minus_one}, 0, 2, admit_prescribed_path=True)
    raise AssertionError('opposite prescription reused incompatible evidence')
except integration.EvaluationError:
    pass
assert len(upper_bank) == count
control = integration.ComputationControl()
control.cancel()
try:
    upper.evaluate(upper_bank, {s:minus_two}, 0, 2, admit_prescribed_path=True, control=control)
    raise AssertionError('cancellation was ignored')
except integration.CalculationCancelled:
    pass
assert len(upper_bank) == count

rooted = make_flow('+i0', True)
root_bank = integration.BoundaryCache()
seed(rooted, root_bank, sheets={r:1})
root_result = rooted.evaluate(root_bank, {s:minus_one}, 0, 2,
    root_sheets={r:1}, admit_prescribed_path=True)
assert root_result.root_sheets == {r:1}
root_count = len(root_bank)
try:
    rooted.evaluate(root_bank, {s:minus_one}, 0, 2,
        root_sheets={r:-1}, admit_prescribed_path=True)
    raise AssertionError('inconsistent destination root sheet was accepted')
except integration.EvaluationError:
    pass
assert len(root_bank) == root_count
rooted_lower = make_flow('-i0', True)
root_lower_bank = integration.BoundaryCache()
seed(rooted_lower, root_lower_bank, sheets={r:1})
root_lower = rooted_lower.evaluate(root_lower_bank, {s:minus_one}, 0, 2,
    root_sheets={r:-1}, admit_prescribed_path=True)
assert root_lower.root_sheets == {r:-1}

endpoint_flow = integration.KinematicTransport(eps, {s:[[one/s]]}, [Y], one,
    branch_domain='y=s on an upper approach', options=opts, continuation=declaration('+i0'))
endpoint_bank = integration.BoundaryCache()
seed(endpoint_flow, endpoint_bank, at=minus_one, values=[[nm1]])
route = integration.EndpointRoute(z, {s:z}, half, homotopy='positive radial approach after upper detour')
try:
    endpoint_flow.evaluate_endpoint(endpoint_bank, route, 0, 0, admit_endpoint=True)
    raise AssertionError('matching homotopy was not admitted')
except integration.EvaluationError:
    pass
assert len(endpoint_bank) == 1 and endpoint_bank.endpoint_entries() == []
endpoint = endpoint_flow.evaluate_endpoint(endpoint_bank, route, 0, 0,
    admit_matching_path=True, admit_endpoint=True)
assert 20 <= endpoint.verified_digits <= endpoint.input_verified_digits <= 60
assert endpoint.matching_boundary.steps > 0
assert len(endpoint_bank.endpoint_entries()) == 1
with tempfile.TemporaryDirectory() as directory:
    endpoint_bank.save(directory)
    restored = integration.BoundaryCache.load(directory)
    terminal = endpoint_flow.evaluate_endpoint(restored, route, 0, 0, admit_endpoint=True)
    assert terminal.cache_hit and terminal.coefficients == endpoint.coefficients
    assert terminal.comparison_errors == endpoint.comparison_errors
    assert terminal.matching_boundary.steps == 0
    merged = integration.BoundaryCache()
    merged.extend(restored)
    assert len(merged) == len(restored) and len(merged.endpoint_entries()) == 1
    try:
        endpoint_flow.evaluate_endpoint(restored, route, 0, 0)
        raise AssertionError('terminal cache hit bypassed endpoint admission')
    except integration.EvaluationError:
        pass

for bad in ['upper', 'insensitive']:
    try:
        integration.ContinuationPrescription([(s,bad)], domain='invalid')
        raise AssertionError('bad prescription accepted')
    except integration.InvalidInputError:
        pass
for polynomial in [eps, r, one/s, zero]:
    try:
        bad = integration.ContinuationPrescription([(polynomial,'+i0')], domain='invalid')
        integration.KinematicTransport(eps, {s:[[eps/s]]}, [Y], one,
            branch_domain='invalid', continuation=bad)
        raise AssertionError('native declaration admission was bypassed')
    except (integration.InvalidInputError, integration.UnsupportedInputError):
        pass
", Some(&local), Some(&local))?;
        let value = |name: &str, order: usize| -> PyResult<_> {
            local.get_item(name)?.unwrap().getattr("coefficients")?.get_item(order)?.get_item(0)?.extract::<PythonMultiPrecisionComplex>()
        };
        let ipi = symbolica_amflow::ComplexFloat::new(p.real(0), p.real(1).pi());
        assert!(p.close(&value("up",1)?.0, &ipi, 20));
        assert!(p.close(&value("down",1)?.0, &p.neg(&ipi), 20));
        assert!(p.close(&value("up",2)?.0, &p.scale(&p.mul(&ipi,&ipi),1,2),20));
        assert!(p.close(&value("root_result",1)?.0,&p.scale(&ipi,1,2),20));
        assert!(p.close(&value("root_lower",1)?.0,&p.scale(&ipi,-1,2),20));
        assert!(p.close(&value("endpoint",0)?.0,&p.zero(),20));
        Ok(())
    }).unwrap();
}

#[cfg(feature = "python_stubgen")]
#[test]
fn prescribed_transport_stubs_keep_native_expressions_and_descriptive_names() {
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
    for needle in [
        "class ContinuationPrescription",
        "continuation:",
        "admit_prescribed_path:",
        "class EndpointRoute",
        "evaluate_endpoint(",
    ] {
        assert!(text.contains(needle), "missing {needle}");
    }
    assert!(text.contains("from symbolica.core import Expression"));
    assert!(text.contains("prescriptions: typing.Sequence[tuple[Expression, builtins.str]]"));
}
