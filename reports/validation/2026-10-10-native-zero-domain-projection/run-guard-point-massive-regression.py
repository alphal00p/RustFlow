"""Generate one frozen-build massive profile without reading reference values."""
import json,subprocess,sys
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
PREFIX='guard-point'
snapshot=BASE/(PREFIX+'-source-hashes.json');verify_snapshot(snapshot)
build_path=BASE/(PREFIX+'-build-provenance.json');build=json.loads(build_path.read_text())
assert build['source_snapshot_sha256']==digest(snapshot)
summary=BASE/(PREFIX+'-root-gates.json');gates=json.loads(summary.read_text())
assert gates['status']=='passed' and gates['build_provenance']['sha256']==digest(build_path)
assert {x['gate'] for x in gates['gates']}=={'native','lib','flow-boundary','massless-sources','source-fingerprint','python','capacity'}
old=ROOT/'reports/validation/2026-10-09-finite-density-native-assembly/requested-rays-massive-regression-resources.provenance.json'
command=json.loads(old.read_text())['command']
report=BASE/'guard-point-massive-regression';resources=BASE/'guard-point-massive-regression-resources.json'
if report.exists() or resources.exists():raise ValueError('preserve prior run outputs')
old_assignment=next(x for x in command if x.startswith('RUSTFLOW_DENSITY_FLOW_REPORT='))
command[command.index(old_assignment)]='RUSTFLOW_DENSITY_FLOW_REPORT='+str(report.relative_to(ROOT))
executable=next(x for x in command if x.startswith('target/release/deps/finite_density_runtime_flow-'))
assert digest(ROOT/executable)==build['files'][executable]['sha256']
input_name=next(x.split('=',1)[1] for x in command if x.startswith('RUSTFLOW_WEIGHTED_INPUT='))
assert digest(ROOT/input_name)==json.loads(snapshot.read_text())['files'][input_name]
import os
environment={k:v for k,v in os.environ.items() if not k.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_'))}
binding={'scope':'One current-source massive two-loop profile; prediction generation only, no reference values read. Numerical comparison follows saved output.','source_snapshot_sha256':digest(snapshot),'build_provenance_sha256':digest(build_path),'compile_identities':build['compile_identities'],'input':{'path':input_name,'sha256':digest(ROOT/input_name)},'executable':{'path':executable,**build['files'][executable]},'baseline_command':{'path':str(old.relative_to(ROOT)),'sha256':digest(old)},'command':command,'launcher_sha256':digest(Path(__file__)),'inherited_weighted_and_density_environment_removed':True}
(BASE/'guard-point-massive-regression-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
result=subprocess.run([sys.executable,str(BASE/'run-resource-command.py'),str(resources),*command],cwd=ROOT,env=environment)
verify_snapshot(snapshot);assert digest(ROOT/executable)==build['files'][executable]['sha256']
(BASE/'guard-point-massive-regression-post-run-binding.json').write_text(json.dumps({'exit_code':result.returncode,'source_snapshot_sha256':digest(snapshot),'build_provenance_sha256':digest(build_path),'executable_sha256':digest(ROOT/executable),'unchanged_source_and_executable_verified':True,'resources_sha256':digest(resources)},indent=2)+'\n')
sys.exit(result.returncode)
