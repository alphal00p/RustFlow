#!/usr/bin/env python3
"""Save production three-loop predictions with the optional requested-ray-point schedule.

Uses only the physical input and actual native source factory. No experimental
program or rule is imported. The 3600-second full fixed-D cap allows two
singleton cuts and the double cut after a source-only singleton needed about
664 seconds. This is a resource allowance, not a speed comparison.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot

PREFIX='requested-ray-point'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gate',choices=['fixed','laurent','double'])
    args=parser.parse_args()
    snapshot_path=BASE/(PREFIX+'-source-hashes.json')
    snapshot=verify_snapshot(snapshot_path)
    build_path=BASE/(PREFIX+'-build-provenance.json')
    build=json.loads(build_path.read_text())
    assert build['source_snapshot_sha256']==digest(snapshot_path)
    gates_path=BASE/(PREFIX+'-root-gates.json')
    gates=json.loads(gates_path.read_text())
    assert gates['status']=='passed'and gates['build_provenance']['sha256']==digest(build_path)
    assert any(row['gate']=='polynomial-sources'and row['status']=='passed'for row in gates['gates'])
    matches=[name for name in build['files']if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert len(matches)==1
    executable=matches[0]
    assert digest(ROOT/executable)==build['files'][executable]['sha256']
    name=PREFIX+'-three-loop-'+args.gate
    report=BASE/name
    resources=BASE/(name+'-resources.json')
    binding_path=BASE/(name+'-binding.json')
    if report.exists()or binding_path.exists()or any(BASE.glob(name+'-resources.*')):
        raise ValueError('preserve previous attempts; use a distinct reviewed launcher prefix')
    settings={
        'RUSTFLOW_WEIGHTED_INPUT':'examples/finite_density/massless_three_loop_chain.json',
        'RUSTFLOW_WEIGHTED_SOURCE_POLICY':'polynomial-closure-v1',
        'RUSTFLOW_WEIGHTED_REQUESTED_RAYS':'true','RUSTFLOW_WEIGHTED_FREE_VIRTUAL_ZEROS':'true',
        'RUSTFLOW_WEIGHTED_ZERO_FACES':'true','RUSTFLOW_WEIGHTED_PRIORITIZE_REQUESTED':'false',
        'RUSTFLOW_WEIGHTED_REUSED_RULES':'65536','RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'0',
        'RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT':'true','RUSTFLOW_WEIGHTED_RAY_DOMAINS':'1',
        'RUSTFLOW_WEIGHTED_POINT_DOMAINS':'1',
        'RUSTFLOW_WEIGHTED_ACTIVE_TARGET_CLOSURE':'true','RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'8192',
        'RUSTFLOW_WEIGHTED_GUARD_PASSES':'0','RUSTFLOW_WEIGHTED_DEPTH':'3','RUSTFLOW_WEIGHTED_DOMAINS':'8192',
        'RUSTFLOW_WEIGHTED_ROUNDS':'16','RUSTFLOW_WEIGHTED_FRONTIER':'1024','RUSTFLOW_WEIGHTED_REQUESTED':'4096',
        'RUSTFLOW_DENSITY_FLOW_EPSILON':'-5/4',
        'RUSTFLOW_DENSITY_FLOW_PROFILES':'18:60:8'if args.gate=='double'else'18:60:8,28:60:8,28:80:8,28:80:12',
        'RUSTFLOW_DENSITY_FLOW_REPORT':str(report.relative_to(ROOT)),
    }
    if args.gate=='double':settings['RUSTFLOW_WEIGHTED_CUTS']='0,3'
    input_path=ROOT/settings['RUSTFLOW_WEIGHTED_INPUT']
    assert digest(input_path)==snapshot['files'][settings['RUSTFLOW_WEIGHTED_INPUT']]
    test={'fixed':'runtime_graph_full_amplitude','laurent':'runtime_graph_full_laurent',
          'double':'runtime_graph_occupied_flow'}[args.gate]
    timeout_seconds={'fixed':3600,'laurent':1800,'double':600}[args.gate]
    command=['timeout',str(timeout_seconds),executable,
             '--ignored','--exact',test,'--nocapture','--test-threads=1']
    inherited={key:value for key,value in os.environ.items()
               if key.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_FLOW_'))}
    environment={key:value for key,value in os.environ.items()if key not in inherited}
    environment.update(settings)
    binding={'source_snapshot_sha256':digest(snapshot_path),'build_provenance_sha256':digest(build_path),
        'root_gates_sha256':digest(gates_path),'launcher_sha256':digest(Path(__file__)),
        'input_sha256':digest(input_path),'executable':executable,
        'executable_sha256':digest(ROOT/executable),'settings':settings,'command':command,
        'cleared_inherited_runtime_variable_names':sorted(inherited),
        'scope':'Bounded actual RustFlow flow with polynomial-closure-v1 sources, optional requested-ray-point scheduling and unchanged native reducer. No experimental rule/program imports. Three-loop numerical profiles remain unchanged. The 3600-second full fixed-D allowance accounts for two singleton cuts and a double cut after a source-only singleton required approximately 664 seconds; no speed comparison is claimed.',
        'timeout_seconds':timeout_seconds,
        'requested_ray_point_schedule_enabled':True,
        'experimental_programs_or_rules_imported':False,
        'predictions_before_reference_access':True,'experimental_native_portfolio_enabled':False}
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    result=subprocess.run([sys.executable,str(BASE/'run-resource-command.py'),str(resources),*command],cwd=ROOT,env=environment)
    verify_snapshot(snapshot_path)
    assert digest(ROOT/executable)==build['files'][executable]['sha256']
    binding['exit_code']=result.returncode
    binding['post_run_source_and_executable_unchanged']=True
    binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
