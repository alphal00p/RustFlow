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
fn exact_python_endpoint_relations_preserve_native_evidence_and_cache_identity() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let local = PyDict::new(py);
        local.set_item("integration", integration)?;
        for key in ["s", "eps", "z", "Y1", "Y2", "r"] {
            local.set_item(key, PythonExpression::from(Atom::parse(key, "python_constraints", Default::default()).unwrap()))?;
        }
        for (key, text) in [("one", "1"), ("zero", "0"), ("half", "1/2"), ("tiny", "1/10^300")] {
            local.set_item(key, PythonExpression::from(Atom::parse(text, "python_constraints", Default::default()).unwrap()))?;
        }
        local.set_item("imag", PythonExpression::from(Atom::num(Complex::new(Rational::zero(), Rational::one()))))?;
        let p = Precision::decimal(90).unwrap();
        let amplitude = symbolica_amflow::ComplexFloat::new(p.real(1), p.rational(&Rational::from((1, 3))).re);
        local.set_item("amplitude", PythonMultiPrecisionComplex(amplitude.clone()))?;
        local.set_item("n1", PythonMultiPrecisionComplex(p.i(1)))?;
        local.set_item("n0", PythonMultiPrecisionComplex(p.zero()))?;
        local.set_item("error0", PythonMultiPrecisionFloat(p.real(0)))?;
        py.run(c"
import tempfile
C, R, Constraints = integration.AsymptoticCoefficient, integration.AsymptoticRelation, integration.EndpointConstraints
select = C(0, zero, 1)
relation = R([(select, one)], zero)
assert (select.component, select.power, select.log_power) == (0, zero, 1)
assert relation.value == zero and relation.terms[0][1] == one
constraints = Constraints([relation], 'explicit absence of the logarithmic mode')
assert constraints.provenance == 'explicit absence of the logarithmic mode'
assert constraints.max_dimension == 64 and constraints.max_order == 256
assert constraints.max_coefficient_bits == 65536 and constraints.max_scalar_cells == 4000000
assert constraints.relations[0].terms[0][0].log_power == 1
opts = integration.EvaluationOptions(digits=20, guard_digits=30, series_order=48)
route = integration.EndpointRoute(z, {s:z}, half, homotopy='positive radial endpoint')
flow = integration.KinematicTransport(eps, {s:[[zero,one/s],[zero,zero]]}, [Y1,Y2], one,
    branch_domain='correlated finite endpoint', options=opts)
bank = integration.BoundaryCache()
def seed(f, b, values=None):
    values = values or [[amplitude,n0]]
    f.add_boundary(b, {s:one}, values, 0, verified_digits=60,
        comparison_errors=[[error0,error0] for row in values], provenance='independent constant solution')
seed(flow, bank)
def endpoint(f, b, asserted=constraints, last=0, **kw):
    return f.evaluate_endpoint(b, route, 0, last, constraints=asserted,
        admit_matching_path=True, admit_endpoint=True, **kw)
try:
    endpoint(flow, bank, None)
    raise AssertionError('rounded zero was treated as an exact relation')
except integration.AccuracyError:
    pass
assert len(bank) == 1 and bank.endpoint_entries() == []
answer = endpoint(flow, bank)
assert answer.constraints.provenance == constraints.provenance
assert answer.constraints.relations[0].terms[0][0].power == zero
assert 20 <= answer.verified_digits <= answer.input_verified_digits == 60
assert not answer.cache_hit and len(bank.endpoint_entries()) == 1
assert not isinstance(answer.coefficients[0][0], (float,complex))
with tempfile.TemporaryDirectory() as directory:
    path = directory + '/boundary.bin'
    bank.save(path)
    restored = integration.BoundaryCache.load(path)
    again = endpoint(flow, restored)
    assert again.cache_hit and again.coefficients == answer.coefficients
    assert again.comparison_errors == answer.comparison_errors
    assert again.constraints.provenance == constraints.provenance
    assert restored.endpoint_entries()[0].constraints.provenance == constraints.provenance
    # Limits govern work rather than mathematical identity of a completed record.
    limited = Constraints([relation], constraints.provenance, max_order=2)
    assert endpoint(flow, restored, limited).cache_hit
    distinct = Constraints([relation], 'separately asserted identical physical relation')
    assert not endpoint(flow, restored, distinct).cache_hit
    assert len(restored.endpoint_entries()) == 2
    count = len(restored)
    try:
        flow.evaluate_endpoint(restored, route, 0, 0, constraints=constraints)
        raise AssertionError('terminal cache hit bypassed endpoint admission')
    except integration.InvalidInputError:
        pass
    control = integration.ComputationControl()
    control.cancel()
    try:
        endpoint(flow, restored, control=control)
        raise AssertionError('cancelled request reused a terminal')
    except integration.CalculationCancelled:
        pass
    assert len(restored) == count and len(restored.endpoint_entries()) == 2

for bad, expected_error in [(Constraints([R([(select,one)],tiny)], 'tiny exact divergent coefficient'), integration.AccuracyError),
                            (Constraints([relation, R([(select,one)],one)], 'contradictory assertions'), integration.InvalidInputError)]:
    count = len(bank)
    try:
        endpoint(flow, bank, bad)
        raise AssertionError('inconsistent or divergent exact constraints accepted')
    except expected_error:
        pass
    assert len(bank) == count and len(bank.endpoint_entries()) == 1
mismatch = integration.BoundaryCache()
seed(flow, mismatch, [[amplitude,n1]])
try:
    endpoint(flow, mismatch)
    raise AssertionError('physical relations hid an inconsistent regular component')
except integration.AccuracyError:
    pass
assert len(mismatch) == 1 and mismatch.endpoint_entries() == []

for rows, provenance in [([], 'empty'), ([R([],zero)], 'empty row'), ([relation], '')]:
    try:
        Constraints(rows, provenance)
        raise AssertionError('malformed declaration accepted')
    except integration.InvalidInputError:
        pass
for bad in [R([(C(0,imag),one)],zero), R([(select,s)],zero), R([(select,one)],s)]:
    try:
        Constraints([bad], 'non-Gaussian declaration')
        raise AssertionError('unsupported exact field accepted')
    except integration.UnsupportedInputError:
        pass
try:
    Constraints([relation], 'resource policy', max_dimension=0)
    raise AssertionError('invalid resource limit accepted')
except integration.InvalidInputError:
    pass
outside = Constraints([R([(C(2,zero,1),one)],zero)], 'out-of-range endpoint component')
try:
    endpoint(flow, bank, outside)
    raise AssertionError('selector lower bound replaced full endpoint preflight')
except integration.InvalidInputError:
    pass

# Canonical connections use the same flattened epsilon hierarchy and exact core.
canonical = integration.KinematicTransport.canonical(eps, [s], [s], [[[zero,one],[zero,zero]]],
    [Y1,Y2], one, branch_domain='canonical hierarchy', options=opts)
canonical_bank = integration.BoundaryCache()
seed(canonical, canonical_bank, [[amplitude,n0],[n1,n0]])
hierarchy = Constraints([R([(C(2,zero,1),one)],zero)], 'zero leading second component')
canonical_answer = endpoint(canonical, canonical_bank, hierarchy, last=1)
assert canonical_answer.leading_power == 0 and canonical_answer.last_power == 1

prescription = integration.ContinuationPrescription([(s,'+i0')], domain='upper endpoint')
prescribed = integration.KinematicTransport(eps, {s:[[zero,one/s],[zero,zero]]}, [Y1,Y2], one,
    branch_domain='prescribed constrained endpoint', options=opts, continuation=prescription)
prescribed_bank = integration.BoundaryCache()
seed(prescribed, prescribed_bank)
prescribed_answer = endpoint(prescribed, prescribed_bank)
assert endpoint(prescribed, prescribed_bank).cache_hit

# Unsupported algebraic constrained spaces must remain an explicit native failure.
rooted = integration.KinematicTransport(eps, {s:[[zero,r/s],[zero,zero]]}, [Y1,Y2], one,
    roots={r:s}, branch_domain='root constraint admission', options=opts)
root_bank = integration.BoundaryCache()
rooted.add_boundary(root_bank, {s:one}, [[amplitude,n0]], 0, verified_digits=60,
    comparison_errors=[[error0,error0]], provenance='constant root-system solution', root_sheets={r:1})
root_route = integration.EndpointRoute(z, {s:z}, half, homotopy='positive root ray', root_sheets={r:1})
try:
    rooted.evaluate_endpoint(root_bank, root_route, 0, 0, constraints=constraints,
        admit_matching_path=True, admit_endpoint=True)
    raise AssertionError('unproved root constraint space was admitted')
except integration.UnsupportedInputError:
    pass
assert len(root_bank) == 1 and root_bank.endpoint_entries() == []
", Some(&local), Some(&local))?;
        for name in ["answer", "canonical_answer", "prescribed_answer"] {
            let row = local.get_item(name)?.unwrap().getattr("coefficients")?.get_item(0)?;
            let actual = row.get_item(0)?.extract::<PythonMultiPrecisionComplex>()?;
            let zero = row.get_item(1)?.extract::<PythonMultiPrecisionComplex>()?;
            assert!(p.close(&actual.0, &amplitude, 20), "{name}");
            assert!(p.close(&zero.0, &p.zero(), 20), "{name}");
        }
        Ok(())
    }).unwrap();
}

#[cfg(feature = "python_stubgen")]
#[test]
fn constrained_endpoint_stubs_use_native_exact_expression_types() {
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
        "class AsymptoticCoefficient",
        "class AsymptoticRelation",
        "class EndpointConstraints",
        "constraints: typing.Optional[EndpointConstraints] = None",
        "power: Expression",
        "value: Expression",
        "max_coefficient_bits: builtins.int = 65536",
        "terms: typing.Sequence[tuple[AsymptoticCoefficient, Expression]]",
    ] {
        assert!(text.contains(needle), "missing {needle}");
    }
}
