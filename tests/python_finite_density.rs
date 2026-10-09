#![cfg(feature = "python")]
//! Execute the registered bindings with the same linked native dependency graph.
//! This does not alter or require an independently installed embedding host.
use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};
use symbolica::api::python::create_symbolica_module;

#[test]
fn native_density_binding_preserves_full_values_errors_and_cancellation() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "symbolica.community.hepkit")?;
        feynkit_py::initialize_feynkit(&hep)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let local = PyDict::new(py);
        local.set_item("integration", integration)?;
        local.set_item("one_loop", include_str!("../examples/finite_density/massive_one_loop_tadpole.json"))?;
        local.set_item("two_loop", include_str!("../examples/finite_density/massive_two_loop_sunset.json"))?;
        local.set_item("four_loop", include_str!("../examples/finite_density/chain_of_three_parallel_pairs.json"))?;
        py.run(c"
import concurrent.futures
import json
import time
from decimal import Decimal, localcontext

prepared = json.loads(integration.prepare_finite_density(one_loop))
assert prepared['input']['loops'] == 1
options = integration.EvaluationOptions(digits=18, guard_digits=40, series_order=60, mass_mode='all')
control = integration.ComputationControl()
report = json.loads(integration.evaluate_finite_density(one_loop, '1/2', options=options, control=control))
assert report['full_amplitude'] is True and report['normalization'] == 'unscaled Euclidean amplitude'
assert report['verified_digits'] is None
assert [part['cut_slots'] for part in report['contributions']] == [[], [0]]
assert report['occupied_reports'][0]['physical_arity'] == 4
assert report['occupied_reports'][0]['native_storage_capacity'] is None
assert control.poll(), 'native progress was not exposed through the binding'
with localcontext() as ctx:
    ctx.prec = 50
    pi = Decimal('3.1415926535897932384626433832795028841971693993751')
    for value in report['values']:
        assert isinstance(value['real'], str) and isinstance(value['imaginary'], str)
        assert abs(Decimal(value['imaginary'])) < Decimal('1e-20')
    assert abs(Decimal(report['values'][0]['real']) + 1/(4*pi)) < Decimal('1e-16')
    assert abs(Decimal(report['values'][1]['real'])) < Decimal('1e-20')
    assert abs(Decimal(report['contributions'][0]['values'][1]['real']) - 1/(4*pi)) < Decimal('1e-16')
    assert abs(Decimal(report['contributions'][1]['values'][1]['real']) + 1/(4*pi)) < Decimal('1e-16')

# No coefficient answer enters fitting: the native owner checks a second grid.
laurent = json.loads(integration.evaluate_finite_density(one_loop))
assert laurent['full_amplitude'] is True and laurent['operation'] == 'finite-density'
assert all(expansion['verified_digits'] >= 12 and expansion['validation_samples'] > 0
           for expansion in laurent['expansions'])
assert all(set(expansion['coefficients']) == {'-1', '0'} for expansion in laurent['expansions'])

for epsilon in ['0', 'unknown_epsilon']:
    try:
        integration.evaluate_finite_density(one_loop, epsilon)
        raise AssertionError('invalid regulator was accepted')
    except integration.InvalidInputError:
        pass
try:
    integration.evaluate_finite_density(one_loop, '1/2', start_scale=3)
    raise AssertionError('invalid start scale was accepted')
except integration.InvalidInputError:
    pass
massless = json.loads(four_loop)
for edge in massless['edges']:
    edge['mass_squared'] = '0'
try:
    integration.evaluate_finite_density(json.dumps(massless), '1/7')
    raise AssertionError('unsupported massless flowing sector was accepted')
except integration.UnsupportedInputError:
    pass

# A second Python thread must be able to observe and cancel native preparation.
control = integration.ComputationControl()
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
    future = pool.submit(integration.evaluate_finite_density, two_loop, '4/5', options=options, control=control)
    deadline = time.monotonic() + 30
    while not control.poll():
        assert not future.done(), 'calculation finished before progress was observable'
        assert time.monotonic() < deadline, 'GIL held during native evaluation'
        time.sleep(0.001)
    control.cancel()
    try:
        future.result(timeout=30)
        raise AssertionError('cancelled native calculation was accepted')
    except integration.CalculationCancelled:
        pass
assert control.cancelled
", Some(&local), Some(&local))?;
        Ok(())
    }).unwrap();
}
