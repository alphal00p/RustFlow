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
fn general_rooted_python_connections_preserve_native_transport_and_evidence() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let local = PyDict::new(py);
        local.set_item("integration", integration)?;
        let p = Precision::decimal(100).unwrap();
        for key in ["s", "t", "eps", "r", "q", "z", "Y", "Z"] {
            local.set_item(key, PythonExpression::from(Atom::parse(key, "python_dense_roots", Default::default()).unwrap()))?;
        }
        for (key, value) in [("one", 1), ("zero", 0), ("two", 2), ("three", 3), ("four", 4), ("eight", 8), ("nine", 9), ("minus_one", -1)] {
            local.set_item(key, PythonExpression::from(Atom::num(value)))?;
        }
        local.set_item("quarter", PythonExpression::from(Atom::num(Rational::from((1, 4)))))?;
        local.set_item("sixteenth", PythonExpression::from(Atom::num(Rational::from((1, 16)))))?;
        for (key, value) in [("n0", p.zero()), ("n1", p.i(1)), ("n3", p.i(3)), ("exp1", p.exp(&p.i(1)))] {
            local.set_item(key, PythonMultiPrecisionComplex(value))?;
        }
        local.set_item("error0", PythonMultiPrecisionFloat(p.real(0)))?;
        py.run(c"
import inspect
import tempfile
opts = integration.EvaluationOptions(digits=20, guard_digits=30, series_order=48)
assert 'roots' in str(inspect.signature(integration.KinematicTransport))
def seed(flow, bank, at, values, sheets):
    return flow.add_boundary(bank, at, values, 0, verified_digits=60,
        comparison_errors=[[error0 for _ in row] for row in values],
        provenance='independent analytic fundamental matrix', root_sheets=sheets)
def declaration(side):
    return integration.ContinuationPrescription([(s, side)], domain=side+' analytic continuation', unprescribed_side=side)
def coupled(variables, radicand, continuation=None):
    matrix = [[eps/r, one/r], [zero, eps/r]]
    return integration.KinematicTransport(eps, {v:matrix for v in variables}, [Y,Z], one,
        roots={r:radicand}, branch_domain='exp(2 epsilon (r-r0)) times a triangular matrix',
        options=opts, continuation=continuation)

flow = coupled([s,t], s+t)
assert flow.dimension == 2 and flow.roots == {r:s+t}
bank = integration.BoundaryCache()
initial = [[n0,n1], [n0,n0], [n0,n0]]
seed(flow, bank, {s:one,t:zero}, initial, {r:1})
first = flow.evaluate(bank, {s:three,t:one}, 0, 2, root_sheets={r:1}, admit_straight_path=True)
assert 20 <= first.verified_digits <= first.input_verified_digits == 60
assert first.steps > 0 and first.inserted_points > 0
assert first.root_sheets == {r:1}
with tempfile.TemporaryDirectory() as directory:
    bank.save(directory)
    restored = integration.BoundaryCache.load(directory)
    nearby = flow.evaluate(restored, {s:eight,t:one}, 0, 2, root_sheets={r:1}, admit_straight_path=True)
    assert nearby.starting_coordinates == {s:three,t:one}
    assert nearby.input_verified_digits <= first.verified_digits
    hit = flow.evaluate(restored, {s:eight,t:one}, 0, 2, root_sheets={r:1})
    assert hit.cache_hit and hit.steps == 0 and hit.coefficients == nearby.coefficients
    assert hit.comparison_errors == nearby.comparison_errors
    count = len(restored)
    for sheets in [None, {}, {r:0}, {q:1}, {r:-1}]:
        try:
            flow.evaluate(restored, {s:eight,t:one}, 0, 2, root_sheets=sheets, admit_straight_path=True)
            raise AssertionError('missing or inconsistent root germ accepted')
        except integration.EvaluationError:
            pass
        assert len(restored) == count
    control = integration.ComputationControl()
    control.cancel()
    try:
        flow.evaluate(restored, {s:nine,t:one}, 0, 2, root_sheets={r:1}, admit_straight_path=True, control=control)
        raise AssertionError('cancelled dense rooted transport ran')
    except integration.CalculationCancelled:
        pass
    assert len(restored) == count

opposite_bank = integration.BoundaryCache()
seed(flow, opposite_bank, {s:one,t:zero}, initial, {r:-1})
opposite = flow.evaluate(opposite_bank, {s:three,t:one}, 0, 2, root_sheets={r:-1}, admit_straight_path=True)

# Nontrivial native quotient denominator, two roots and two physical derivatives.
quotient = integration.KinematicTransport(eps,
    {s:[[one/(two*r*(r+q))]], t:[[one/(two*q*(r+q))]]}, [Y], one,
    roots={r:s, q:t}, branch_domain='positive roots of both invariants', options=opts)
quotient_bank = integration.BoundaryCache()
seed(quotient, quotient_bank, {s:one,t:four}, [[n3]], {r:1,q:1})
quotient_result = quotient.evaluate(quotient_bank, {s:four,t:nine}, 0, 0,
    root_sheets={r:1,q:1}, admit_straight_path=True)
count = len(quotient_bank)
try:
    quotient.evaluate(quotient_bank, {s:four,t:four}, 0, 0,
        root_sheets={r:1,q:1}, admit_straight_path=True)
    raise AssertionError('original quotient norm domain was erased')
except integration.EvaluationError:
    pass
assert len(quotient_bank) == count

for side, sheets, result_name in [('+i0',{r:1},'upper'),('-i0',{r:-1},'lower')]:
    prescribed = coupled([s], s, declaration(side))
    prescribed_bank = integration.BoundaryCache()
    seed(prescribed, prescribed_bank, {s:one}, initial, {r:1})
    value = prescribed.evaluate(prescribed_bank, {s:minus_one}, 0, 2,
        root_sheets=sheets, admit_prescribed_path=True)
    assert value.root_sheets == sheets
    assert 20 <= value.verified_digits <= value.input_verified_digits == 60
    globals()[result_name] = value
    hit = prescribed.evaluate(prescribed_bank, {s:minus_one}, 0, 2, root_sheets=sheets)
    assert hit.cache_hit and hit.coefficients == value.coefficients

for continuation in [None, declaration('+i0')]:
    endpoint_flow = integration.KinematicTransport(eps, {s:[[one/r]]}, [Y], one,
        roots={r:s}, branch_domain='exp(2 sqrt(s))', options=opts, continuation=continuation)
    endpoint_bank = integration.BoundaryCache()
    seed(endpoint_flow, endpoint_bank, {s:quarter}, [[exp1]], {r:1})
    route = integration.EndpointRoute(z, {s:z}, sixteenth,
        root_sheets={r:1}, homotopy='positive radial approach')
    endpoint = endpoint_flow.evaluate_endpoint(endpoint_bank, route, 0, 0,
        admit_matching_path=True, admit_endpoint=True)
    assert 20 <= endpoint.verified_digits <= endpoint.input_verified_digits <= 60
    assert len(endpoint_bank.endpoint_entries()) == 1
    terminal = endpoint_flow.evaluate_endpoint(endpoint_bank, route, 0, 0, admit_endpoint=True)
    assert terminal.cache_hit and terminal.coefficients == endpoint.coefficients

rational_args = (eps, {s:[[eps/s]]}, [Y], one)
legacy = integration.KinematicTransport(*rational_args, branch_domain='rational')
empty = integration.KinematicTransport(*rational_args, branch_domain='rational', roots={})
assert legacy.identity == empty.identity and legacy.roots == empty.roots == {}
for bad_roots in [{r:eps}, {r:r}, {r:zero}, {s:s}, {one:s}]:
    try:
        integration.KinematicTransport(eps, {s:[[one/r]]}, [Y], one,
            roots=bad_roots, branch_domain='invalid roots')
        raise AssertionError('native root validation bypassed')
    except integration.EvaluationError:
        pass
", Some(&local), Some(&local))?;
        let coefficients = |name: &str| -> PyResult<Vec<Vec<PythonMultiPrecisionComplex>>> {
            local.get_item(name)?.unwrap().getattr("coefficients")?.extract()
        };
        for (name, expected) in [
            ("first", [[2,1],[4,2],[4,2]]),
            ("nearby", [[4,1],[16,4],[32,8]]),
            ("opposite", [[-2,1],[4,-2],[-4,2]]),
        ] {
            let actual = coefficients(name)?;
            assert_eq!(actual.len(), expected.len());
            for (row, want) in actual.iter().zip(expected) {
                assert_eq!(row.len(), want.len());
                for (a,b) in row.iter().zip(want) { assert!(p.close(&a.0,&p.i(b),20), "{name}"); }
            }
        }
        for (name, sign) in [("upper",1),("lower",-1)] {
            let shift = p.complex(-2,2*sign);
            let mut exponential = p.i(1);
            let actual = coefficients(name)?;
            assert_eq!(actual.len(), 3);
            for (order, row) in actual.iter().enumerate() {
                assert_eq!(row.len(), 2);
                assert!(p.close(&row[0].0, &p.mul(&shift,&exponential),20),"{name} order {order}");
                assert!(p.close(&row[1].0, &exponential,20),"{name} order {order}");
                exponential=p.scale(&p.mul(&exponential,&shift),1,(order+1) as i64);
            }
        }
        assert!(p.close(&coefficients("quotient_result")?[0][0].0,&p.i(5),20));
        assert!(p.close(&coefficients("endpoint")?[0][0].0,&p.i(1),20));
        Ok(())
    }).unwrap();
}

#[cfg(feature = "python_stubgen")]
#[test]
fn dense_root_declarations_appear_in_native_expression_stubs() {
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
    let constructor = text
        .lines()
        .find(|line| {
            line.trim_start()
                .starts_with("def __new__(cls, epsilon: Expression, derivatives:")
        })
        .unwrap();
    assert!(
        constructor
            .contains("roots: typing.Optional[typing.Mapping[Expression, Expression]] = None"),
        "dense constructor lacks native roots: {constructor}"
    );
    assert!(text.contains("def roots(self) -> builtins.dict[Expression, Expression]:"));
}
