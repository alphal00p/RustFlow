#!/usr/bin/env python3
"""Fresh longer bounded singleton run after a timed-out contracting frontier."""
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
K=ROOT/'reports/validation/2026-10-10-requested-ray-point-integration'
PREVIOUS=ROOT/'reports/validation/2026-10-10-larger-frontier-controls'
sys.path.insert(0,str(K))
from frozen_sources import digest,verify_snapshot


def main():
    prior_path=PREVIOUS/'single3-binding.json'
    prior=json.loads(prior_path.read_text())
    assert prior['exit_code']==124 and prior['post_run_source_and_executable_unchanged']
    verify_snapshot(K/'requested-ray-point-source-hashes.json')
    executable=prior['executable']
    assert digest(ROOT/executable)==prior['executable_sha256']
    settings=dict(prior['settings'])
    assert settings['RUSTFLOW_WEIGHTED_ROUNDS']=='16'
    settings['RUSTFLOW_WEIGHTED_ROUNDS']='24'
    settings['RUSTFLOW_DENSITY_FLOW_REPORT']=str((HERE/'single3').relative_to(ROOT))
    binding_path=HERE/'run-binding.json'
    assert not binding_path.exists() and not (HERE/'single3').exists()
    command=['timeout','3600',executable,'--ignored','--exact',
             'runtime_graph_occupied_flow','--nocapture','--test-threads=1']
    inherited={key for key in os.environ if key.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_FLOW_'))}
    environment={key:value for key,value in os.environ.items() if key not in inherited}
    environment.update(settings)
    binding={'scope':'Fresh bounded production run, not a checkpoint resume. The prior1200-second run contracted from2510 to872 frontier labels before timeout. Raise time cap to3600 seconds and round cap16→24; keep all physical input, sources, native discovery and other budgets unchanged. All prior cost/failure remains recorded. No source/rule imports, speed or full-amplitude claim.',
             'predecessor_binding':{'path':str(prior_path.relative_to(ROOT)),'sha256':digest(prior_path)},
             'predecessor_resources_sha256':digest(PREVIOUS/'single3-resources.json'),
             'source_snapshot_sha256':prior['source_snapshot_sha256'],
             'build_provenance_sha256':prior['build_provenance_sha256'],
             'input_sha256':prior['input_sha256'],
             'executable':executable,'executable_sha256':prior['executable_sha256'],
             'launcher_sha256':digest(Path(__file__)),'command':command,'settings':settings,
             'cleared_inherited_runtime_variable_names':sorted(inherited),
             'checkpoint_resume':False,'original_rules_imported':False}
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    result=subprocess.run([sys.executable,str(K/'run-resource-command.py'),str(HERE/'resources.json'),*command],cwd=ROOT,env=environment)
    verify_snapshot(K/'requested-ray-point-source-hashes.json')
    assert digest(ROOT/executable)==binding['executable_sha256']
    binding.update(exit_code=result.returncode,post_run_source_and_executable_unchanged=True)
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    raise SystemExit(result.returncode)


if __name__=='__main__':main()
