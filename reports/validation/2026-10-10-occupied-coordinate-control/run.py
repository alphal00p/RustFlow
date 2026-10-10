#!/usr/bin/env python3
"""Exact unimodular input rerouting; unchanged native selected-cut diagnostic."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
K = ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration'
sys.path.insert(0, str(K))
from frozen_sources import digest, verify_snapshot


def main():
    source = ROOT / 'examples/finite_density/massless_three_loop_chain.json'
    original = json.loads(source.read_text())
    # Old momenta P = R Q. R is its own inverse with determinant -1.
    routing = [[1, 0, -1], [0, 1, -1], [0, 0, -1]]
    assert [[sum(routing[i][k] * routing[k][j] for k in range(3))
             for j in range(3)] for i in range(3)] == [[1,0,0],[0,1,0],[0,0,1]]
    permutation = [3, 4, 2, 0, 1]
    transformed = json.loads(source.read_text())
    transformed['name'] += '_exact_occupied_coordinate_control'
    transformed['edges'] = []
    for old_slot in permutation:
        edge = dict(original['edges'][old_slot])
        old_row = list(map(int, edge['routing']))
        edge['routing'] = [str(sum(old_row[k] * routing[k][j] for k in range(3)))
                           for j in range(3)]
        transformed['edges'].append(edge)
    transformed['targets'][0]['powers'] = [original['targets'][0]['powers'][i] for i in permutation]
    transformed['targets'][1]['powers'] = [original['targets'][1]['powers'][i] for i in permutation]
    assert original['targets'][1]['numerator'] == 'g1_2+u1*u2'
    transformed['targets'][1]['numerator'] = 'g1_2-g1_3-g2_3+g3_3+(u1-u3)*(u2-u3)'
    # The only charged old loop is P1; P1 = Q1-Q3, and Q3 has zero charge.
    assert original['loop_charges'] == [[1],[0],[0]]
    assert [[sum(routing[i][j]*transformed['loop_charges'][j][0] for j in range(3))]
            for i in range(3)] == original['loop_charges']
    input_path = HERE / 'rerouted-input.json'
    binding_path = HERE / 'run-binding.json'
    assert not input_path.exists() and not binding_path.exists()
    input_path.write_text(json.dumps(transformed, indent=2)+'\n')
    snapshot_path = K / 'requested-ray-point-source-hashes.json'
    build_path = K / 'requested-ray-point-build-provenance.json'
    verify_snapshot(snapshot_path)
    build = json.loads(build_path.read_text())
    executable, = [name for name in build['files'] if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert digest(ROOT/executable) == build['files'][executable]['sha256']
    settings = dict(json.loads((K/'requested-ray-point-three-loop-fixed-binding.json').read_text())['settings'])
    settings.update({
        'RUSTFLOW_WEIGHTED_INPUT': str(input_path.relative_to(ROOT)),
        'RUSTFLOW_WEIGHTED_CUTS': '0',
        'RUSTFLOW_DENSITY_FLOW_PROFILES': '18:60:8',
        'RUSTFLOW_DENSITY_FLOW_REPORT': str((HERE/'selected-cut').relative_to(ROOT)),
    })
    inherited = {key for key in os.environ if key.startswith(('RUSTFLOW_WEIGHTED_', 'RUSTFLOW_DENSITY_FLOW_'))}
    environment = {key:value for key,value in os.environ.items() if key not in inherited}
    environment.update(settings)
    command = ['timeout', '600', executable, '--ignored', '--exact',
               'runtime_graph_occupied_flow', '--nocapture', '--test-threads=1']
    binding = {
        'scope': 'Exact input rerouting and physical-edge permutation of original cut [3], represented as new cut [0]. Both original targets retained. No full-amplitude claim or speed comparison; no source/reducer changes or imported rules.',
        'original_input': str(source.relative_to(ROOT)), 'original_sha256': digest(source),
        'rerouted_input_sha256': digest(input_path),
        'old_momenta_in_new': routing, 'determinant': -1, 'absolute_jacobian': 1,
        'old_physical_slot_for_new': permutation,
        'numerator_substitution': {'P1':'Q1-Q3','P2':'Q2-Q3','P3':'-Q3'},
        'graph_incidence_preserved_per_edge': True, 'chemical_assignments_preserved': True,
        'build_provenance_sha256': digest(build_path), 'source_snapshot_sha256': digest(snapshot_path),
        'launcher_sha256': digest(Path(__file__)), 'executable': executable,
        'executable_sha256': digest(ROOT/executable), 'settings': settings, 'command': command,
        'cleared_inherited_runtime_variable_names': sorted(inherited),
    }
    binding_path.write_text(json.dumps(binding, indent=2)+'\n')
    result = subprocess.run([sys.executable, str(K/'run-resource-command.py'),
                             str(HERE/'resources.json'), *command], cwd=ROOT, env=environment)
    verify_snapshot(snapshot_path)
    assert digest(ROOT/executable) == binding['executable_sha256']
    binding.update(exit_code=result.returncode, post_run_source_and_executable_unchanged=True)
    binding_path.write_text(json.dumps(binding, indent=2)+'\n')
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
