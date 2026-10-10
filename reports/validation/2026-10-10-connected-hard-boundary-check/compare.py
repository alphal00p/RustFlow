#!/usr/bin/env python3
"""Saved-output comparison after all native profiles have completed."""
from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
sys.path.insert(0, str(ROOT / 'tools/finite_density'))
from compare_massive_reference import complex_decimal


def digest(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    output = BASE / 'comparison.json'
    assert not output.exists(), 'preserve existing comparisons'
    build_path, run_path = BASE / 'native-build-binding.json', BASE / 'native-run-binding.json'
    build, run = read(build_path), read(run_path)
    assert build['exit_code'] == run['exit_code'] == 0
    assert run['post_run_bindings_unchanged'] is True
    for binding in [build, run]:
        for name, sha in binding['inputs_sha256'].items():
            assert digest(Path(name)) == sha
    resources_path = BASE / 'native-resources.json'
    assert read(resources_path)['exit_code'] == 0
    directory = BASE / 'native-predictions'
    assert not (directory / 'failure.json').exists()
    configurations = [(18, 60), (28, 60), (28, 80)]
    paths = [directory / f'prediction-{digits}-{order}.json' for digits, order in configurations]
    assert set(directory.glob('prediction-*.json')) == set(paths)
    records = [read(path) for path in paths]
    for row, (digits, order) in zip(records, configurations):
        assert row['dimension'] == '12/5' and row['epsilon'] == '4/5'
        assert row['powers'] == [1, 1, 1]
        assert row['digits'] == digits and row['series_order'] == order and row['guard_digits'] == 40
        assert row['terminal_policy'] == 'TadpolesOnly' and row['bubble_subloops'] is False
        assert row['native_preparations'] > 0 and row['reference_files_read'] == 0
    # Reference access starts only after completed native output checks.
    reference_path = BASE / 'reference.json'
    reference = read(reference_path)
    assert reference['source_sha256'] == digest(BASE / 'reference.py')
    assert reference['oracle_records_read'] == reference['native_prediction_files_read'] == 0
    comparisons, refinements = [], []
    with localcontext() as context:
        context.prec = 80
        expected = D(reference['reference_native'])
        previous = None
        for row in records:
            real, imaginary = complex_decimal(row['value'])
            relative = abs(real - expected) / abs(expected)
            assert relative < D('1e-10') and abs(imaginary) < D('1e-12')
            comparisons.append({'digits': row['digits'], 'series_order': row['series_order'],
                                'relative_real_error': str(relative), 'absolute_imaginary_error': str(abs(imaginary))})
            if previous is not None:
                change = max(abs(real - previous[0]), abs(imaginary - previous[1])) / abs(expected)
                assert change < D('1e-12')
                refinements.append(str(change))
            previous = (real, imaginary)
    checked = [build_path, run_path, resources_path, *paths, reference_path,
               BASE / 'reference.py', BASE / 'native.rs', BASE / 'derivation.md', Path(__file__),
               ROOT / 'tools/finite_density/compare_massive_reference.py']
    output.write_text(json.dumps({'status': 'passed',
        'scope': 'Independent 10-digit check of one convergent ordinary connected hard coefficient, evaluated recursively without bubble or sunset terminals. No full three-loop finite-density acceptance claim.',
        'comparisons': comparisons, 'relative_native_refinements': refinements,
        'inputs_sha256': {str(p.relative_to(ROOT)): digest(p) for p in checked},
        'definitions_or_tolerances_fitted': False,
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
