"""Copy the frozen runtime harness and join its two prediction phases.

The phase bodies are copied verbatim. Only their duplicate Run/flow preparation
is removed; both use the same freshly prepared public native flow. No reference
or persisted rule is loaded. Generation alone does not start a calculation.
"""
from pathlib import Path
import hashlib
import json

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
ORIGINAL = ROOT / 'tests/finite_density_runtime_flow.rs'

def digest(data):
    return hashlib.sha256(data).hexdigest()

def phase(source, name):
    marker = 'fn ' + name + '() {'
    assert source.count(marker) == 1
    tail = source.split(marker, 1)[1]
    body = tail.split('\n#[test]', 1)[0].rstrip()
    assert body.endswith('}')
    body = body[:-1]
    prefix = '\n    let run = Run::new();\n    let flow = run.full();'
    assert body.startswith(prefix)
    return body[len(prefix):]

def main():
    source = ORIGINAL.read_text()
    fixed = phase(source, 'runtime_graph_full_amplitude')
    laurent = phase(source, 'runtime_graph_full_laurent')
    suffix = '''

#[test]
#[ignore = "fresh full three-loop preparation shared by fixed-D and Laurent predictions"]
fn runtime_graph_joint_predictions() {
    let run = Run::new();
    run.save("joint-harness.json", json!({
        "prepared_flow_count": 1,
        "phase_order": ["fixed_dimension", "laurent"],
        "original_runtime_phase_bodies": true,
        "persisted_rules_loaded": false,
        "reference_values_loaded": false,
        "phase_configuration": "configuration.json"
    }));
    let flow = run.full();
    run.context.emit(Progress::Stage { name: "joint harness fixed-dimension phase".into() }).unwrap();
    {
''' + fixed + '''
    }
    run.save("fixed-phase-completed.json", json!({"predictions_saved": true, "reference_comparison_performed": false}));
    run.context.emit(Progress::Stage { name: "joint harness Laurent phase".into() }).unwrap();
    {
''' + laurent + '''
    }
    run.save("joint-phases-completed.json", json!({"both_phases_completed": true, "reference_comparison_performed": false}));
}
'''
    generated = source + suffix
    (BASE / 'joint-runtime.rs').write_text(generated)
    binding = {
        'original_path': str(ORIGINAL.relative_to(ROOT)),
        'original_sha256': digest(source.encode()),
        'generated_sha256': digest(generated.encode()),
        'generator_sha256': digest(Path(__file__).read_bytes()),
        'fixed_phase_body_sha256': digest(fixed.encode()),
        'laurent_phase_body_sha256': digest(laurent.encode()),
        'scope': 'verbatim runtime source plus one joined test; one Run::new and one Run::full before unchanged phase bodies',
        'references_loaded': False,
        'run_started': False,
    }
    (BASE / 'harness-binding.json').write_text(json.dumps(binding, indent=2) + '\n')
    print(json.dumps(binding))

if __name__ == '__main__':
    if (BASE / 'harness-binding.json').exists() or (BASE / 'joint-runtime.rs').exists():
        raise ValueError('preserve existing harness; use a new report for changed source')
    main()
