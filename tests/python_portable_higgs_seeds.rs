#![cfg(feature = "python_api")]
use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};
use symbolica::api::python::{PythonExpression, create_symbolica_module};
use symbolica::domains::float::{PythonMultiPrecisionComplex, PythonMultiPrecisionFloat};
use symbolica_amflow::Precision;

#[test]
fn numerical_method_citations_are_process_isolated_and_cumulative() {
    const CHILD: &str = "RUSTFLOW_CITATION_TEST_OPERATION";
    let Ok(operation) = std::env::var(CHILD) else {
        for operation in [
            "import",
            "ordinary",
            "kinematic",
            "higgs",
            "model",
            "projector",
            "automatic",
        ] {
            if operation == "automatic" && !cfg!(feature = "automatic") {
                continue;
            }
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "numerical_method_citations_are_process_isolated_and_cumulative",
                    "--nocapture",
                ])
                .env(CHILD, operation)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{operation}: {}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
        }
        return;
    };
    assert!(symbolica_amflow::python::get_citations().is_empty());
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "symbolica.community.hepkit")?;
        feynkit_py::initialize_feynkit(&hep)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("hep", hep)?;
        locals.set_item("integration", integration)?;
        locals.set_item("operation", &operation)?;
        assert!(symbolica_amflow::python::get_citations().is_empty());
        py.run(
            c"
integration.EvaluationOptions()
integration.BoundaryCache()
integration.ComputationControl()
for invalid in [lambda: integration.HiggsJetIntegralSystem('invalid'),
                lambda: integration.DifferentialSystem(core.S('citation_x'), []),
                lambda: integration.HiggsJetAmplitude.with_form_factor_vertices(hep.Model.phi4())]:
    try:
        invalid()
    except integration.InvalidInputError:
        pass
    else:
        raise AssertionError('invalid API preparation succeeded')
",
            Some(&locals),
            Some(&locals),
        )?;
        assert!(symbolica_amflow::python::get_citations().is_empty());
        py.run(
            c"
if operation == 'ordinary':
    integration.DifferentialSystem(core.S('citation_x'), [[core.E('0')]])
elif operation == 'kinematic':
    integration.KinematicTransport(core.S('citation_eps'), {core.S('citation_x'):[[core.E('0')]]},
        [core.S('citation_master')], core.E('1'), branch_domain='citation preparation')
elif operation == 'higgs':
    integration.HiggsJetIntegralSystem('planar')
elif operation == 'model':
    integration.HiggsJetAmplitude.with_form_factor_vertices(hep.Model.standard_model())
elif operation == 'projector':
    integration.HiggsJetFormFactorProjector()
elif operation == 'automatic':
    integration.IntegralEvaluator()
",
            Some(&locals),
            Some(&locals),
        )?;
        Ok(())
    })
    .unwrap();
    let snapshot = || {
        symbolica_amflow::python::get_citations()
            .into_iter()
            .map(|c| (c.id, c.reference, c.bibtex, c.reasons, c.description, c.url))
            .collect::<Vec<_>>()
    };
    let first = snapshot();
    let expected = match operation.as_str() {
        "import" => vec![],
        "ordinary" | "kinematic" => vec!["arXiv:2607.08477", "arXiv:2006.05510"],
        "automatic" => vec!["arXiv:2607.08477", "arXiv:2006.05510", "arXiv:2201.11669"],
        "higgs" => vec!["arXiv:2607.08477", "arXiv:2006.05510", "arXiv:2112.07578"],
        "model" | "projector" => vec!["arXiv:2112.07578"],
        _ => unreachable!(),
    };
    assert_eq!(
        first.iter().map(|c| c.0.as_str()).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        snapshot(),
        first,
        "reading citations reset cumulative usage"
    );
    if let (Some(amflow2), Some(diffexp)) = (
        first.iter().find(|c| c.0 == "arXiv:2607.08477"),
        first.iter().find(|c| c.0 == "arXiv:2006.05510"),
    ) {
        assert_ne!(amflow2.3, diffexp.3, "method-specific citation reasons");
        assert_ne!(amflow2.4, diffexp.4, "paper-specific descriptions");
    }
    for (id, reference, bibtex, reasons, description, url) in first {
        assert!(!reference.is_empty() && bibtex.starts_with("@article{"));
        assert!(!reasons.is_empty());
        assert!(!description.is_empty());
        assert_eq!(
            url,
            format!("https://arxiv.org/abs/{}", id.strip_prefix("arXiv:").unwrap())
        );
    }
}

#[test]
fn higgs_form_factor_extension_returns_the_shared_native_model() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "symbolica.community.hepkit")?;
        feynkit_py::initialize_feynkit(&hep)?;
        let integration = PyModule::new(py, "symbolica.community.hep.integration")?;
        symbolica_amflow::python::register(&integration)?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("hep", hep)?;
        locals.set_item("integration", integration)?;
        py.run(
            c"
import json
model = hep.Model.standard_model()
original = model.to_json()
extended = integration.HiggsJetAmplitude.with_form_factor_vertices(model)
assert type(extended) is hep.Model
assert model.to_json() == original
assert extended.particle_by_pdg(21).name == 'g'
assert extended.particle_by_pdg(25).name == 'H'
assert isinstance(extended.parameter('GGGHEWWW_ForFac1_RE').symbol, core.Expression)
assert isinstance(extended.coupling('GGGH_HEFT_C1').expression, core.Expression)
for name in ['GGGHEWWW', 'GGGHEWZZ', 'GGGHHEFT']:
    assert extended.vertex_rule(name).particles == ['g', 'g', 'g', 'H']
# The returned object is accepted directly by the native process owner.
process = extended.process(['g', 'g'], ['g', 'H'],
    vertex_allow=['GGGHEWWW', 'GGGHEWZZ', 'GGGHHEFT'])
assert isinstance(process, hep.Process)
extended_json = extended.to_json()
try:
    integration.HiggsJetAmplitude.with_form_factor_vertices(extended)
    raise AssertionError('conflicting effective declarations were overwritten')
except integration.InvalidInputError as error:
    assert 'already exists' in str(error)
assert extended.to_json() == extended_json and model.to_json() == original
try:
    integration.HiggsJetAmplitude.with_form_factor_vertices(hep.Model.phi4())
    raise AssertionError('a model without the required native species was accepted')
except integration.InvalidInputError:
    pass
restored = hep.Model.from_json(extended_json)
assert json.loads(restored.to_json()) == json.loads(extended_json)
",
            Some(&locals),
            Some(&locals),
        )?;
        let extended = locals.get_item("extended")?.unwrap();
        let native = extended.extract::<PyRef<'_, feynkit_py::PyModel>>()?;
        let owner = std::sync::Arc::clone(native.as_model());
        let second_view = feynkit_py::PyModel::from(std::sync::Arc::clone(&owner));
        assert!(std::sync::Arc::ptr_eq(
            native.as_model(),
            second_view.as_model()
        ));
        assert_eq!(owner.particle_by_pdg(21).unwrap().name, "g");
        let original = locals.get_item("model")?.unwrap();
        let original = original.extract::<PyRef<'_, feynkit_py::PyModel>>()?;
        assert_eq!(
            owner.vertex_rules().len(),
            original.as_model().vertex_rules().len() + 3
        );
        Ok(())
    })
    .unwrap();
}

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
    assert!(
        text.contains("@staticmethod\n    def with_form_factor_vertices(model: Model) -> Model:"),
        "{text}"
    );
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
