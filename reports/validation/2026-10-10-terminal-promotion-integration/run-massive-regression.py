#!/usr/bin/env python3
"""Save the existing single massive profile with its unchanged source policy.

The frozen prior command supplies the input, epsilon, precision, order, start
scale and ordinary closure defaults. The same explicit polynomial source policy is retained, and the requested
ray/point schedule stays disabled. Only the report destination and rebuilt
executable change. This launcher reads no reference values.
"""
import json
import os
import subprocess
import sys
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot

PREFIX='terminal-promotion'

def main():
    snapshot_path=BASE/(PREFIX+'-source-hashes.json')
    snapshot=verify_snapshot(snapshot_path)
    build_path=BASE/(PREFIX+'-build-provenance.json')
    build=json.loads(build_path.read_text())
    assert build['source_snapshot_sha256']==digest(snapshot_path)
    gates_path=BASE/(PREFIX+'-root-gates.json')
    gates=json.loads(gates_path.read_text())
    assert gates['status']=='passed'and gates['build_provenance']['sha256']==digest(build_path)
    assert any(row['gate']=='polynomial-sources'and row['status']=='passed'for row in gates['gates'])
    prior=ROOT/'reports/validation/2026-10-09-finite-density-native-assembly/requested-rays-massive-regression-resources.provenance.json'
    command=json.loads(prior.read_text())['command']
    name=PREFIX+'-massive-regression'
    report=BASE/name
    resources=BASE/(name+'-resources.json')
    binding_path=BASE/(name+'-binding.json')
    if report.exists()or binding_path.exists()or any(BASE.glob(name+'-resources.*')):
        raise ValueError('preserve prior outputs')
    assignment=next(x for x in command if x.startswith('RUSTFLOW_DENSITY_FLOW_REPORT='))
    command[command.index(assignment)]='RUSTFLOW_DENSITY_FLOW_REPORT='+str(report.relative_to(ROOT))
    assert not any(x.startswith('RUSTFLOW_WEIGHTED_SOURCE_POLICY=')for x in command)
    native_executable=next(x for x in command if x.startswith('target/release/deps/finite_density_runtime_flow-'))
    matches=[name for name in build['files']if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert len(matches)==1
    executable=matches[0]
    command[command.index(native_executable)]=executable
    command.insert(command.index('env')+1,'RUSTFLOW_WEIGHTED_SOURCE_POLICY=polynomial-closure-v1')
    assert digest(ROOT/executable)==build['files'][executable]['sha256']
    settings={x.split('=',1)[0]:x.split('=',1)[1]for x in command if x.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_'))}
    assert settings['RUSTFLOW_WEIGHTED_INPUT']=='examples/finite_density/massive_two_loop_sunset.json'
    assert settings['RUSTFLOW_DENSITY_FLOW_EPSILON']=='4/5'
    assert settings['RUSTFLOW_DENSITY_FLOW_PROFILES']=='28:80:12'
    input_path=ROOT/settings['RUSTFLOW_WEIGHTED_INPUT']
    assert digest(input_path)==snapshot['files'][settings['RUSTFLOW_WEIGHTED_INPUT']]
    environment={key:value for key,value in os.environ.items()if not key.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_'))}
    binding={'scope':'One existing massive two-loop profile with optional polynomial-closure-v1; predictions only. Independent/historical comparison follows saved outputs. This single profile is not a new precision-refinement gate.',
        'source_snapshot_sha256':digest(snapshot_path),'build_provenance_sha256':digest(build_path),
        'root_gates_sha256':digest(gates_path),'compile_identities':build['compile_identities'],
        'input':{'path':settings['RUSTFLOW_WEIGHTED_INPUT'],'sha256':digest(input_path)},
        'executable':{'path':executable,**build['files'][executable]},
        'baseline_command':{'path':str(prior.relative_to(ROOT)),'sha256':digest(prior)},
        'settings':settings,'command':command,'launcher_sha256':digest(Path(__file__)),
        'inherited_weighted_and_density_environment_removed':True,
        'requested_ray_point_schedule_enabled':False,
        'predictions_before_reference_access':True,'experimental_native_portfolio_enabled':False}
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    result=subprocess.run([sys.executable,str(BASE/'run-resource-command.py'),str(resources),*command],cwd=ROOT,env=environment)
    verify_snapshot(snapshot_path)
    assert digest(ROOT/executable)==build['files'][executable]['sha256']
    binding['exit_code']=result.returncode
    binding['post_run_source_and_executable_unchanged']=True
    binding['resources_sha256']=digest(resources)
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
