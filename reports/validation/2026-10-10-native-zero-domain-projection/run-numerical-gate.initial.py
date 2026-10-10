#!/usr/bin/env python3
"""Run bounded predictions from a frozen build. No reference data are loaded."""
import argparse,json,subprocess,sys
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('gate',choices=['two-loop-fixed','two-loop-laurent','three-loop-fixed','three-loop-priority-fixed','three-loop-double','three-loop-laurent'])
args=parser.parse_args()
snapshot_path=BASE/'zero-projection-source-hashes.json'
snapshot=verify_snapshot(snapshot_path)
build_path=BASE/'zero-projection-build-provenance.json';build=json.loads(build_path.read_text())
assert build['source_snapshot_sha256']==digest(snapshot_path)
exe=[name for name in build['files'] if Path(name).name.startswith('finite_density_runtime_flow-')]
assert len(exe)==1;exe=exe[0]
assert digest(ROOT/exe)==build['files'][exe]['sha256']
three=args.gate.startswith('three-');laurent=args.gate.endswith('laurent');double=args.gate.endswith('double');priority='priority' in args.gate
name='zero-projection-'+args.gate;report=BASE/name;resources=BASE/(name+'-resources.json')
if report.exists() or any(BASE.glob(name+'-resources.*')):raise ValueError('preserve existing artifacts; select a new gate name')
settings={
 'RUSTFLOW_WEIGHTED_INPUT':'examples/finite_density/massless_'+('three_loop_chain' if three else 'two_loop_sunset')+'.json',
 'RUSTFLOW_WEIGHTED_REQUESTED_RAYS':'true','RUSTFLOW_WEIGHTED_FREE_VIRTUAL_ZEROS':'true',
 'RUSTFLOW_WEIGHTED_ZERO_FACES':'true','RUSTFLOW_WEIGHTED_PRIORITIZE_REQUESTED':str(priority).lower(),
 'RUSTFLOW_WEIGHTED_REUSED_RULES':'65536','RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'32',
 'RUSTFLOW_WEIGHTED_ACTIVE_TARGET_CLOSURE':'true','RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'262144',
 'RUSTFLOW_WEIGHTED_GUARD_PASSES':'0','RUSTFLOW_WEIGHTED_DEPTH':'3','RUSTFLOW_WEIGHTED_DOMAINS':'8192',
 'RUSTFLOW_WEIGHTED_ROUNDS':'16','RUSTFLOW_WEIGHTED_FRONTIER':'1024','RUSTFLOW_WEIGHTED_REQUESTED':'4096',
 'RUSTFLOW_DENSITY_FLOW_EPSILON':'-5/4',
 'RUSTFLOW_DENSITY_FLOW_PROFILES':'18:60:8' if double else '18:60:8,28:60:8,28:80:8,28:80:12',
 'RUSTFLOW_DENSITY_FLOW_REPORT':str(report.relative_to(ROOT)),
}
if double:settings['RUSTFLOW_WEIGHTED_CUTS']='0,3'
input_path=ROOT/settings['RUSTFLOW_WEIGHTED_INPUT']
assert digest(input_path)==snapshot['files'][settings['RUSTFLOW_WEIGHTED_INPUT']]
test='runtime_graph_occupied_flow' if double else 'runtime_graph_full_laurent' if laurent else 'runtime_graph_full_amplitude'
command=['timeout','1800' if laurent and three else '600','env',*(k+'='+v for k,v in settings.items()),exe,'--ignored','--exact',test,'--nocapture','--test-threads=1']
binding={'source_snapshot_sha256':digest(snapshot_path),'build_provenance_sha256':digest(build_path),'launcher_sha256':digest(Path(__file__)),'input_sha256':digest(input_path),'executable':exe,'executable_sha256':digest(ROOT/exe),'settings':settings,'command':command,'predictions_before_reference_access':True}
binding_path=BASE/(name+'-binding.json');binding_path.write_text(json.dumps(binding,indent=2)+'\n')
result=subprocess.run([sys.executable,str(BASE/'run-resource-command.py'),str(resources),*command],cwd=ROOT)
verify_snapshot(snapshot_path);assert digest(ROOT/exe)==build['files'][exe]['sha256']
binding['exit_code']=result.returncode;binding['post_run_source_and_executable_unchanged']=True
binding_path.write_text(json.dumps(binding,indent=2)+'\n')
sys.exit(result.returncode)
