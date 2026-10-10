"""Bind completed coherent lower-loop checks as launch prerequisites only."""
import hashlib
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / '2026-10-10-singleton-physical-zero-integration-v2'
PREFIX = 'singleton-physical-zero-v2'

def h(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def main():
    assert not (HERE / 'numerical-prerequisites.json').exists()
    cap = json.loads((HERE / 'capsule-binding.json').read_text())
    rows = []
    for name, refs, refinements in [('two-loop-fixed', 40, 30), ('two-loop-laurent', 30, 24), ('massive-regression', 10, 0)]:
        stem = PREFIX + '-' + name
        cp = PARENT / (stem + '-comparison.json')
        bp = PARENT / (stem + '-comparison-binding.json')
        rp = PARENT / (stem + '-binding.json')
        comparison = json.loads(cp.read_text())
        binding = json.loads(bp.read_text())
        run = json.loads(rp.read_text())
        assert comparison['status'] == binding['status'] == 'passed'
        assert binding['exit_code'] == 0 and binding['comparison_sha256'] == h(cp)
        assert comparison['reference_comparison_count'] == refs
        assert comparison['independent_refinement_comparison_count'] == refinements
        assert run['exit_code'] == 0 and run['post_run_source_and_executable_unchanged'] is True
        assert run['source_snapshot_sha256'] == cap['source_snapshot_sha256']
        assert run['build_provenance_sha256'] == cap['build_provenance_sha256']
        assert binding['native_provenance']['source_snapshot_sha256'] == cap['source_snapshot_sha256']
        if name == 'massive-regression':
            assert comparison['historical_same_profile_comparison_count'] == 10
            assert len(comparison['assembly_checks']) == 2
        for path in [cp, bp, rp]:
            copy = HERE / 'evidence' / path.name
            assert not copy.exists()
            shutil.copy2(path, copy)
            assert h(copy) == h(path)
            rows.append({'source_path': str(path), 'copy_path': str(copy.relative_to(HERE)),
                         'sha256': h(copy), 'bytes': copy.stat().st_size})
    record = {
        'scope': 'Completed production-library lower-loop comparisons are launch prerequisites. Their values are never passed to the joint native prediction process.',
        'captured_files': rows, 'source_snapshot_sha256': cap['source_snapshot_sha256'],
        'build_provenance_sha256': cap['build_provenance_sha256'],
        'all_three_comparisons_passed': True, 'capture_script_sha256': h(Path(__file__)),
    }
    (HERE / 'numerical-prerequisites.json').write_text(json.dumps(record, indent=2) + '\n')
    print(h(HERE / 'numerical-prerequisites.json'))

if __name__ == '__main__':
    main()
