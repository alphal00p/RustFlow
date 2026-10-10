#!/usr/bin/env python3
"""Save bounded production predictions using the optional polynomial source policy.

No reference values, experimental rules or isolated native portfolio are loaded.
The algebraic and numerical limits/profiles are retained. The full fixed-D
timeout is 1800 seconds because a successful source-only singleton preparation
already needed about 664 seconds; this run is not a matched speed benchmark.
Comparison is a separate step after saved predictions.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot

PREFIX='polynomial-closure-v2'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gate',choices=['fixed'])
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
    name=PREFIX+'-three-loop-'+args.gate+'-guards3-per-residual1'
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
        'RUSTFLOW_WEIGHTED_REUSED_RULES':'65536','RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'1',
        'RUSTFLOW_WEIGHTED_ACTIVE_TARGET_CLOSURE':'true','RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'262144',
        'RUSTFLOW_WEIGHTED_GUARD_PASSES':'3','RUSTFLOW_WEIGHTED_DEPTH':'3','RUSTFLOW_WEIGHTED_DOMAINS':'8192',
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
    timeout_seconds=600 if args.gate=='double'else 1800
    command=['timeout',str(timeout_seconds),executable,
             '--ignored','--exact',test,'--nocapture','--test-threads=1']
    inherited={key:value for key,value in os.environ.items()
               if key.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_FLOW_'))}
    environment={key:value for key,value in os.environ.items()if key not in inherited}
    environment.update(settings)
    predecessor=BASE/(PREFIX+'-three-loop-fixed-guards3-binding.json')
    predecessor_resources=BASE/(PREFIX+'-three-loop-fixed-guards3-resources.json')
    failed=json.loads(predecessor.read_text())
    assert failed['exit_code']!=0 and failed['post_run_source_and_executable_unchanged']
    assert failed['source_snapshot_sha256']==digest(snapshot_path)
    assert failed['build_provenance_sha256']==digest(build_path)
    expected=dict(failed['settings']);expected['RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL']='1';expected['RUSTFLOW_DENSITY_FLOW_REPORT']=settings['RUSTFLOW_DENSITY_FLOW_REPORT']
    assert expected==settings
    assert failed['command']==command
    binding={'predecessor_failed_binding':{'path':str(predecessor.relative_to(ROOT)),'sha256':digest(predecessor)},
        'predecessor_failed_resources':{'path':str(predecessor_resources.relative_to(ROOT)),'sha256':digest(predecessor_resources)},
        'retry_reason':'The previous three-pass attempt exhausted the shared domain allocation before later exact-point obligations were reached. This existing-option scheduling control changes only the per-residual domain budget from 32 to 1, retaining the 8192 shared budget and three conditional-point passes. It does not add new search or fallback logic.',
        'changed_search_setting':{'RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':{'before':'32','after':'1'}},
        'source_snapshot_sha256':digest(snapshot_path),'build_provenance_sha256':digest(build_path),
        'root_gates_sha256':digest(gates_path),'launcher_sha256':digest(Path(__file__)),
        'input_sha256':digest(input_path),'executable':executable,
        'executable_sha256':digest(ROOT/executable),'settings':settings,'command':command,
        'cleared_inherited_runtime_variable_names':sorted(inherited),
        'scope':'Bounded actual RustFlow production flow with optional polynomial-closure-v1 sources and unchanged native reducer. The previous guards3 retry settings are retained except that domains per residual is reduced from 32 to 1; numerical profiles are unchanged. The full fixed-D timeout is enlarged to 1800 seconds based on the approximately 664-second successful source-only singleton cost; no matched speed claim is made.',
        'timeout_seconds':timeout_seconds,
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
