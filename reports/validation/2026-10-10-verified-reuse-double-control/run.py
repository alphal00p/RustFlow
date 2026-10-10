"""One fresh selected double-cut run using a captured immutable VI executable."""
import datetime,hashlib,json,os,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
h=lambda p: hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<8*1024*1024 else file_hash(p)
def file_hash(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
cap_record=HERE/'capsule-binding.json';capsule=json.loads(cap_record.read_text());CAP=Path(capsule['capsule_path'])
def verify_static():
 for row in capsule['static_files']:assert file_hash(Path(row['path']))==row['sha256'],row['path']
 for row in capsule['evidence_copies']:assert file_hash(HERE/row['copy_path'])==row['sha256'],row['copy_path']
 assert file_hash(HERE/capsule['source_archive']['path'])==capsule['source_archive']['sha256']
verify_static()
assert not (HERE/'run-binding.json').exists() and not (HERE/'double').exists()
report=HERE/'double'
settings={
 'RUSTFLOW_WEIGHTED_INPUT':str(CAP/'input.json'),
 'RUSTFLOW_WEIGHTED_SOURCE_POLICY':'polynomial-closure-v1',
 'RUSTFLOW_WEIGHTED_REQUESTED_RAYS':'true','RUSTFLOW_WEIGHTED_FREE_VIRTUAL_ZEROS':'true',
 'RUSTFLOW_WEIGHTED_ZERO_FACES':'true','RUSTFLOW_WEIGHTED_PRIORITIZE_REQUESTED':'false',
 'RUSTFLOW_WEIGHTED_REUSED_RULES':'65536','RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'0',
 'RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT':'true','RUSTFLOW_WEIGHTED_RAY_DOMAINS':'1',
 'RUSTFLOW_WEIGHTED_POINT_DOMAINS':'1','RUSTFLOW_WEIGHTED_ACTIVE_TARGET_CLOSURE':'true',
 'RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'8192','RUSTFLOW_WEIGHTED_GUARD_PASSES':'0',
 'RUSTFLOW_WEIGHTED_DEPTH':'3','RUSTFLOW_WEIGHTED_DOMAINS':'32768',
 'RUSTFLOW_WEIGHTED_ROUNDS':'32','RUSTFLOW_WEIGHTED_FRONTIER':'8192','RUSTFLOW_WEIGHTED_REQUESTED':'16384',
 'RUSTFLOW_WEIGHTED_CUTS':'0,3','RUSTFLOW_DENSITY_FLOW_EPSILON':'-5/4',
 'RUSTFLOW_DENSITY_FLOW_PROFILES':'18:60:8','RUSTFLOW_DENSITY_FLOW_REPORT':str(report),
}
inherited={k for k in os.environ if k.startswith(('RUSTFLOW_WEIGHTED_','RUSTFLOW_DENSITY_FLOW_'))}
environment={k:v for k,v in os.environ.items()if k not in inherited};environment.update(settings)
command=['timeout','1800',str(CAP/'runtime'),'--ignored','--exact','runtime_graph_occupied_flow','--nocapture','--test-threads=1']
record={'scope':'Fresh original three-loop selected cut [0,3], all uncut physical factors deformed, original two targets. No checkpoint/rules/resume/reference values imported. One fixed-D profile only; no complete three-loop/four-loop acceptance or refinement/speed claim.',
 'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'capsule_binding_sha256':h(cap_record),'launcher_sha256':h(Path(__file__)),
 'source_snapshot_sha256':capsule['source_snapshot_sha256'],'build_provenance_sha256':capsule['build_provenance_sha256'],'root_gates_sha256':capsule['root_gates_sha256'],
 'input_sha256':h(CAP/'input.json'),'settings':settings,'command':command,'cwd':str(CAP),'cleared_inherited_runtime_variable_names':sorted(inherited),
 'global_caps':{'timeout_seconds':1800,'rounds':32,'frontier':8192,'requested':16384,'domains':32768,'reused_rules':65536,'direct_zero_attempts':8192},
 'requested_schedule':{'ray_domains':1,'point_domains':1,'depth':3},'predictions_before_reference_access':True,
 'post_run_checks_scope':'Captured static files/source archive/evidence only, not live production sources or shared target artifacts.'}
(HERE/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n')
result=subprocess.run([sys.executable,str(CAP/'run-resource-command.py'),str(HERE/'resources.json'),*command],cwd=CAP,env=environment)
verify_static()
record.update(exit_code=result.returncode,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),post_run_static_capsule_unchanged=True)
(HERE/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n')
raise SystemExit(result.returncode)
