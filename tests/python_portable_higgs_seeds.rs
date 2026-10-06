#![cfg(feature = "python_api")]
use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};
use symbolica::api::python::{PythonExpression, create_symbolica_module};
use symbolica::domains::float::{PythonMultiPrecisionComplex, PythonMultiPrecisionFloat};
use symbolica_amflow::Precision;

#[test]
fn portable_higgs_connections_preserve_supplied_evidence_and_restart() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let locals = PyDict::new(py);
        locals.set_item("integration", integration)?;
        locals.set_item("automatic", cfg!(feature = "automatic"))?;
        let p = Precision::decimal(100).unwrap();
        locals.set_item("zero", PythonMultiPrecisionComplex(p.zero()))?;
        locals.set_item("error", PythonMultiPrecisionFloat(p.tolerance(45)))?;
        py.run(
            c"
import tempfile
assert integration.automatic_boundary_generation_available == automatic
for name in ['IntegralEvaluator', 'PreparedIntegralFamily', 'ReductionTables']:
    assert hasattr(integration, name) == automatic
    assert (name in integration.__all__) == automatic
certificates = {
    'planar': 'b77f8d97fadc07001ebf3cd87e68ed5338b40a4b0d58fb06a8636956fbd87dc8',
    'nonplanar': 'ab5364501a91e2e7c0960f2f8b84fe14f039321153145e62e55c75ecc4521958',
}
for topology, certificate in certificates.items():
    system = integration.HiggsJetIntegralSystem(topology)
    assert system.automatic_boundary_generation_available == automatic
    assert hasattr(system, 'generate_boundary') == automatic
    flow = system.kinematic_transport(integration.EvaluationOptions(workers=1))
    assert isinstance(flow, integration.KinematicTransport)
    assert system.mathematical_fingerprint == certificate
    renamed = integration.HiggsJetIntegralSystem(topology, namespace='portable_renamed')
    assert renamed.mathematical_fingerprint == certificate
    assert renamed.kinematic_transport().identity != flow.identity
    configuration = system.configurations()[0]
    coefficients = [[zero] * system.dimension for _ in range(5)]
    errors = [[error] * system.dimension for _ in range(5)]
    cache = integration.BoundaryCache()
    supplied = flow.add_boundary(cache, configuration.start, coefficients, 0,
        verified_digits=40, comparison_errors=errors,
        provenance='synthetic supplied values for ownership/evidence plumbing only',
        root_sheets=configuration.root_sheets)
    assert supplied.verified_digits == supplied.input_verified_digits == 40
    assert supplied.coefficients == coefficients and supplied.comparison_errors == errors
    hit = flow.evaluate(cache, configuration.start, 0, 4, root_sheets=configuration.root_sheets)
    assert hit.cache_hit and hit.steps == 0 and hit.coefficients == coefficients
    count = len(cache)
    for bad in [None, {}]:
        try:
            flow.add_boundary(cache, configuration.start, coefficients, 0,
                verified_digits=40, comparison_errors=errors,
                provenance='invalid missing root germ', root_sheets=bad)
            raise AssertionError('root requirements bypassed')
        except integration.InvalidInputError:
            pass
        assert len(cache) == count
    with tempfile.TemporaryDirectory() as directory:
        cache.save(directory)
        restored = integration.BoundaryCache.load(directory)
        same = system.kinematic_transport().evaluate(restored, configuration.start, 0, 4,
            root_sheets=configuration.root_sheets)
        assert same.cache_hit and same.comparison_errors == errors
    for invalid in [integration.EvaluationOptions(dimension=6),
                    integration.EvaluationOptions(prescription='-i0')]:
        try:
            system.kinematic_transport(invalid)
            raise AssertionError('fixed mathematical conventions changed')
        except integration.UnsupportedInputError:
            pass
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
fn portable_higgs_adapter_has_native_generated_stubs() {
    Python::initialize();
    let info = pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica.community.hep.integration".into(),
        env!("CARGO_MANIFEST_DIR").into(),
    )
    .unwrap();
    let text = info.modules["symbolica.community.hep.integration"].to_string();
    assert!(text.contains("automatic_boundary_generation_available: builtins.bool"));
    assert!(text.contains("def mathematical_fingerprint(self) -> builtins.str:"));
    assert!(text.contains("def kinematic_transport(self, options: typing.Optional[EvaluationOptions] = None) -> KinematicTransport:"), "{text}");
}

#[test]
#[ignore = "requires the independently generated community notebook loader and scientific bundle"]
fn all_packaged_higgs_values_enter_the_current_native_evidence_path() {
    let loader = std::env::var("RUSTFLOW_PORTABLE_BUNDLE_LOADER").expect("loader path");
    let bundle = std::env::var("RUSTFLOW_PORTABLE_BUNDLE_DIRECTORY").expect("bundle directory");
    let origin = std::env::var("RUSTFLOW_PORTABLE_ORIGIN_CACHE").expect("origin cache directory");
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let locals = PyDict::new(py);
        locals.set_item("Float", core.getattr("Float")?)?;
        locals.set_item("ComplexFloat", core.getattr("ComplexFloat")?)?;
        locals.set_item("core", core)?;
        locals.set_item("integration", integration)?;
        locals.set_item("automatic", cfg!(feature = "automatic"))?;
        locals.set_item("Expression", py.get_type::<PythonExpression>())?;
        locals.set_item("loader_path", loader)?;
        locals.set_item("bundle_path", bundle)?;
        locals.set_item("origin_cache", origin)?;
        py.run(c"
import gzip, importlib.util, json, pathlib, sys, tempfile, types
# Expose the single statically linked runtime to the independently maintained
# notebook loader; never import an installed wheel with a second Symbolica owner.
package = types.ModuleType('symbolica')
package.__path__ = []
package.Float, package.ComplexFloat, package.E = Float, ComplexFloat, Expression.parse
sys.modules['symbolica'] = package
sys.modules['symbolica.core'] = core
for name in ['symbolica.community', 'symbolica.community.hep']:
    parent = types.ModuleType(name)
    parent.__path__ = []
    sys.modules[name] = parent
sys.modules['symbolica.community.hep.integration'] = integration
spec = importlib.util.spec_from_file_location('portable_bundle_loader', loader_path)
loader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(loader)
systems = {name: integration.HiggsJetIntegralSystem(name) for name in ['planar', 'nonplanar']}
try:
    integration.BoundaryCache.load(origin_cache)
    raise AssertionError('prior-build binary cache admitted by current runtime')
except integration.BoundaryCacheError:
    pass
cache, results = loader.load_boundary_bundle(bundle_path, systems)
assert len(cache) == len(results) == 16
records = json.loads(gzip.decompress((pathlib.Path(bundle_path)/'boundaries.json.gz').read_bytes()))['boundaries']
count = 0
for record in records:
    result = results[record['label']]
    assert result.cache_hit and result.steps == 0
    assert result.verified_digits == result.input_verified_digits == 40
    assert result.working_bits == 415
    assert systems[record['topology']].kinematic_transport().identity != record['origin_identity']
    for row, encoded in zip(result.coefficients, record['coefficients'], strict=True):
        for value, pair in zip(row, encoded, strict=True):
            assert value.real.as_integer_ratio() == tuple(map(int, pair[0][:2]))
            assert value.imag.as_integer_ratio() == tuple(map(int, pair[1][:2]))
            assert value.real.precision == pair[0][2]
            assert value.imag.precision == pair[1][2]
            count += 1
    for row, encoded in zip(result.comparison_errors, record['comparison_errors'], strict=True):
        for value, number in zip(row, encoded, strict=True):
            assert value.as_integer_ratio() == tuple(map(int, number[:2]))
            assert value.precision == number[2]
assert count == 4360
with tempfile.TemporaryDirectory() as directory:
    stored_provenance = [(entry.identity, entry.coordinates, entry.provenance) for entry in cache.entries()]
    cache.save(directory)
    restored = integration.BoundaryCache.load(directory)
    assert [(entry.identity, entry.coordinates, entry.provenance) for entry in restored.entries()] == stored_provenance
    for system in systems.values():
        for configuration in system.configurations():
            hit = system.kinematic_transport().evaluate(restored, configuration.start, 0, 4,
                root_sheets=configuration.root_sheets)
            assert hit.cache_hit and hit.steps == 0
            expected_result = results[configuration.label]
            for attribute in ['verified_digits', 'input_verified_digits', 'working_bits']:
                assert getattr(hit, attribute) == getattr(expected_result, attribute), attribute
            assert expected_result.provenance in hit.provenance
            for row, (actual, expected) in enumerate(zip(hit.coefficients, results[configuration.label].coefficients, strict=True)):
                for column, (a, b) in enumerate(zip(actual, expected, strict=True)):
                    for part in ['real', 'imag']:
                        x, y = getattr(a, part), getattr(b, part)
                        assert x.as_integer_ratio() == y.as_integer_ratio(), (configuration.label, row, column, part, x.as_integer_ratio(), y.as_integer_ratio(), x.precision, y.precision)
                        assert x.precision == y.precision, (configuration.label, row, column, part, x.precision, y.precision)
            for actual, expected in zip(hit.comparison_errors, results[configuration.label].comparison_errors, strict=True):
                for x, y in zip(actual, expected, strict=True):
                    assert x.as_integer_ratio() == y.as_integer_ratio()
                    assert x.precision == y.precision
            assert hit.coefficients == results[configuration.label].coefficients
", Some(&locals), Some(&locals))?;
        Ok(())
    }).unwrap();
}
