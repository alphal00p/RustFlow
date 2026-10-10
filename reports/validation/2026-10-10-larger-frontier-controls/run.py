#!/usr/bin/env python3
"""Bounded frontier-cap control; same original input and native source policy."""
import argparse
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cut', choices=['single3', 'double'])
    args = parser.parse_args()
    snapshot_path = K / 'requested-ray-point-source-hashes.json'
    build_path = K / 'requested-ray-point-build-provenance.json'
    verify_snapshot(snapshot_path)
    build = json.loads(build_path.read_text())
    executable, = [name for name in build['files'] if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert digest(ROOT/executable) == build['files'][executable]['sha256']
    predecessor = K / 'requested-ray-point-three-loop-fixed-binding.json'
    settings = dict(json.loads(predecessor.read_text())['settings'])
    changes = {
        'RUSTFLOW_WEIGHTED_FRONTIER': ('1024', '4096'),
        'RUSTFLOW_WEIGHTED_REQUESTED': ('4096', '16384'),
        'RUSTFLOW_WEIGHTED_DOMAINS': ('8192', '32768'),
    }
    for key, (old, new) in changes.items():
        assert settings[key] == old
        settings[key] = new
    settings.update({
        'RUSTFLOW_WEIGHTED_CUTS': '3' if args.cut == 'single3' else '0,3',
        'RUSTFLOW_DENSITY_FLOW_PROFILES': '18:60:8',
        'RUSTFLOW_DENSITY_FLOW_REPORT': str((HERE/args.cut).relative_to(ROOT)),
    })
    binding_path = HERE / (args.cut+'-binding.json')
    assert not binding_path.exists() and not (HERE/args.cut).exists()
    command = ['timeout', '1200', executable, '--ignored', '--exact',
               'runtime_graph_occupied_flow', '--nocapture', '--test-threads=1']
    inherited = {key for key in os.environ if key.startswith(('RUSTFLOW_WEIGHTED_', 'RUSTFLOW_DENSITY_FLOW_'))}
    environment = {key:value for key,value in os.environ.items() if key not in inherited}
    environment.update(settings)
    binding = {
        'scope': 'Original selected cut, unchanged sources/native search. Increase frontier and associated history/domain caps by four to test whether the earlier frontier truncates later contraction. No rule imports, full-amplitude acceptance or speed comparison.',
        'predecessor_binding_sha256': digest(predecessor),
        'build_provenance_sha256': digest(build_path), 'source_snapshot_sha256': digest(snapshot_path),
        'launcher_sha256': digest(Path(__file__)), 'executable': executable,
        'executable_sha256': digest(ROOT/executable),
        'input_sha256': digest(ROOT/settings['RUSTFLOW_WEIGHTED_INPUT']),
        'changes': {key:{'before':old,'after':new} for key,(old,new) in changes.items()},
        'settings': settings, 'command': command,
        'cleared_inherited_runtime_variable_names': sorted(inherited),
    }
    binding_path.write_text(json.dumps(binding, indent=2)+'\n')
    result = subprocess.run([sys.executable, str(K/'run-resource-command.py'),
                             str(HERE/(args.cut+'-resources.json')), *command], cwd=ROOT, env=environment)
    verify_snapshot(snapshot_path)
    assert digest(ROOT/executable) == binding['executable_sha256']
    binding.update(exit_code=result.returncode, post_run_source_and_executable_unchanged=True)
    binding_path.write_text(json.dumps(binding, indent=2)+'\n')
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
