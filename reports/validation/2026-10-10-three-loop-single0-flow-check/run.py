#!/usr/bin/env python3
"""Fresh original cut-0 native closure, boundary and transport diagnostic."""
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
K = ROOT/'reports/validation/2026-10-10-requested-ray-point-integration'
sys.path.insert(0,str(K))
from frozen_sources import digest, verify_snapshot


def main():
    snapshot_path = K/'requested-ray-point-source-hashes.json'
    build_path = K/'requested-ray-point-build-provenance.json'
    predecessor = K/'requested-ray-point-three-loop-fixed-binding.json'
    verify_snapshot(snapshot_path)
    build=json.loads(build_path.read_text())
    executable, = [name for name in build['files'] if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert digest(ROOT/executable)==build['files'][executable]['sha256']
    baseline=json.loads(predecessor.read_text())
    settings=dict(baseline['settings'])
    settings.update({'RUSTFLOW_WEIGHTED_CUTS':'0',
        'RUSTFLOW_DENSITY_FLOW_PROFILES':'18:60:8',
        'RUSTFLOW_DENSITY_FLOW_REPORT':str((HERE/'selected-cut').relative_to(ROOT))})
    assert settings['RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT']=='true'
    assert settings['RUSTFLOW_WEIGHTED_RAY_DOMAINS']==settings['RUSTFLOW_WEIGHTED_POINT_DOMAINS']=='1'
    assert settings['RUSTFLOW_WEIGHTED_FRONTIER']=='1024'
    assert settings['RUSTFLOW_WEIGHTED_INPUT']=='examples/finite_density/massless_three_loop_chain.json'
    binding_path=HERE/'run-binding.json'
    assert not binding_path.exists() and not (HERE/'selected-cut').exists()
    command=['timeout','600',executable,'--ignored','--exact','runtime_graph_occupied_flow','--nocapture','--test-threads=1']
    inherited={key for key in os.environ if key.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_FLOW_'))}
    environment={key:value for key,value in os.environ.items() if key not in inherited}
    environment.update(settings)
    changes={key:{'before':baseline['settings'].get(key),'after':settings.get(key)}
        for key in sorted(set(settings)|set(baseline['settings'])) if settings.get(key)!=baseline['settings'].get(key)}
    binding={'scope':'Selected original cut [0], both original targets, fresh public factory and native closure/boundary/transport. No prior proof import. One 18:60:8 profile and 600 s cap. Zero-sector agreement alone cannot validate a nonzero hard boundary, full amplitude, precision refinement or performance.',
        'predecessor_binding_sha256':digest(predecessor),
        'build_provenance_sha256':digest(build_path),'source_snapshot_sha256':digest(snapshot_path),
        'launcher_sha256':digest(Path(__file__)),'executable':executable,'executable_sha256':digest(ROOT/executable),
        'input_sha256':digest(ROOT/settings['RUSTFLOW_WEIGHTED_INPUT']),
        'settings':settings,'changes_from_K_full':changes,'command':command,
        'selected_test':'runtime_graph_occupied_flow','full_amplitude':False,
        'cleared_inherited_runtime_variable_names':sorted(inherited)}
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    result=subprocess.run([sys.executable,str(K/'run-resource-command.py'),str(HERE/'resources.json'),*command],cwd=ROOT,env=environment)
    verify_snapshot(snapshot_path)
    assert digest(ROOT/executable)==binding['executable_sha256']
    binding.update(exit_code=result.returncode,post_run_source_and_executable_unchanged=True)
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
