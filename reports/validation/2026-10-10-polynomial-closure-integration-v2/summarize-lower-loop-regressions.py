#!/usr/bin/env python3
"""Summarize completed, independently compared lower-loop regressions only."""
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PREFIX = 'polynomial-closure-v2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def maximum(rows, field):
    values = [Decimal(row[field]) for row in rows if row.get(field) is not None]
    return str(max(values)) if values else None


def main():
    output = BASE / (PREFIX + '-lower-loop-regressions.json')
    assert not output.exists(), 'preserve completed reports'
    snapshot = BASE / (PREFIX + '-source-hashes.json')
    build = BASE / (PREFIX + '-build-provenance.json')
    gates = BASE / (PREFIX + '-root-gates.json')
    rows = []
    bindings = {}
    for tag in ('two-loop-fixed', 'two-loop-laurent', 'massive-regression'):
        stem = PREFIX + '-' + tag
        paths = {kind: BASE / (stem + suffix) for kind, suffix in {
            'binding': '-binding.json', 'resources': '-resources.json',
            'comparison': '-comparison.json',
            'comparison_binding': '-comparison-binding.json',
            'resource_provenance': '-resources.provenance.json',
            'log': '-resources.log',
        }.items()}
        data = {kind: json.loads(path.read_text()) for kind, path in paths.items()
                if kind != 'log'}
        assert data['binding']['exit_code'] == data['resources']['exit_code'] == 0
        assert data['binding']['source_snapshot_sha256'] == sha(snapshot)
        assert data['binding']['build_provenance_sha256'] == sha(build)
        assert data['binding']['post_run_source_and_executable_unchanged']
        assert data['comparison_binding']['status'] == data['comparison']['status'] == 'passed'
        comparison = data['comparison']
        references = comparison.get('comparisons', comparison.get('reference_comparisons'))
        refinements = comparison.get('refinements', [])
        historical = comparison.get('historical_regressions', [])
        assert all(row['passed'] for row in references + refinements + historical
                   + comparison['assembly_checks'])
        counts = (len(references), len(refinements), len(historical))
        assert counts == {'two-loop-fixed': (40, 30, 0),
                          'two-loop-laurent': (30, 24, 0),
                          'massive-regression': (10, 0, 10)}[tag]
        matches = re.findall(r'test result: ok\. 1 passed; 0 failed;.*?finished in ([0-9.]+)s',
                             paths['log'].read_text())
        assert len(matches) == 1
        rows.append({
            'gate': tag, 'status': 'passed',
            'harness_seconds': float(matches[0]),
            'wall_seconds': data['resources']['wall_seconds'],
            'peak_child_rss_kib': data['resources']['peak_child_rss_kib'],
            'reference_comparisons': len(references),
            'independent_profile_refinements': len(refinements),
            'historical_same_profile_comparisons': len(historical),
            'assembly_checks': len(comparison['assembly_checks']),
            'largest_reference_relative_difference': maximum(references, 'relative_difference'),
            'largest_reference_absolute_difference': maximum(references, 'absolute_difference'),
            'largest_refinement_relative_difference': maximum(refinements, 'relative_difference'),
            'largest_historical_relative_difference': maximum(historical, 'relative_difference'),
            'bindings': {kind: {'path': str(path.relative_to(ROOT)), 'sha256': sha(path)}
                         for kind, path in paths.items()},
        })
    for path in (snapshot, build, gates, Path(__file__)):
        bindings[str(path.relative_to(ROOT))] = sha(path)
    report = {
        'schema': 1, 'status': 'passed',
        'scope': 'Lower-loop numerical regression of optional polynomial-closure-v1 using unchanged native RustRed. Complete massless two-loop fixed-D and Laurent profiles, plus one massive two-loop profile. No three-loop or four-loop numerical acceptance follows from this report.',
        'source_policy': 'polynomial-closure-v1',
        'predictions_saved_before_reference_comparison': True,
        'supplied_oracle_numerical_records_compared': 0,
        'reference_uncertainty': 'Empirical independent precision/quadrature refinement and exact formula checks; no rigorous numerical error interval.',
        'performance_scope': 'Prebuilt process resources only, excluding compilation and Nix setup. Concurrent independent runs prevent a controlled speed comparison.',
        'massive_scope': 'Single existing 28-digit/order-80/start-scale-12 profile; no new massive refinement claim.',
        'gates': rows, 'bindings_sha256': bindings,
    }
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'gates': len(rows), 'output': str(output)}))


if __name__ == '__main__':
    main()
