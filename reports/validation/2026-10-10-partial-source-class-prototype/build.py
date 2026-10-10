#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,time,sys,resource,os
root=pathlib.Path(__file__).resolve().parent
repo=root.parents[2]
old=json.loads((repo/'reports/validation/2026-10-10-partial-origin-capability-prototype/build-05.json').read_text())
common=repo/'target/isolated-partial-boundary-validation-libs-20261010/deps'
attempt=sys.argv[1]
out=pathlib.Path('/tmp/rustflow-partial-source-class-20261010');out.mkdir(exist_ok=True)
def sha(p):return hashlib.file_digest(open(p,'rb'),'sha256').hexdigest()
cmd=old['command'][:]
for i,s in enumerate(cmd):
 if s=='partial_origin_prototype':cmd[i]='partial_source_class_prototype'
 if s.startswith('dependency=') and 'target/release/deps' in s:cmd[i]='dependency='+str(common)
 if s.startswith('/tmp/rustflow-partial-origin-20261010/probe.rs'):cmd[i]=str(root/'source/probe.rs')
 if s=='/tmp/rustflow-partial-origin-20261010/prototype-tests':cmd[i]=str(out/('tests-'+attempt))
 if '=/' in s and 'target/release/deps/' in s:
  name,path=s.split('=',1);cmd[i]=name+'='+str(common/pathlib.Path(path).name)
libs={}
for s in cmd:
 if '=/' in s and s.endswith('.rlib'):
  name,path=s.split('=',1);p=pathlib.Path(path);expected=old['libraries'][name]['sha256'];actual=sha(p);assert actual==expected,(name,actual,expected)
  libs[name]={'path':str(p),'sha256':actual,'bytes':p.stat().st_size}
sources={str(p.relative_to(root)):sha(p) for p in sorted((root/'source').rglob('*.rs'))}
report={'command':cmd,'libraries':libs,'source_sha256':sources,'compile_environment':{'CARGO_CRATE_NAME':'partial_source_class_prototype','CARGO_MANIFEST_DIR':str(repo)},'cargo_invoked':False,'shared_target_written':False,'input':{'path':str(repo/'examples/finite_density/massless_three_loop_chain.json'),'sha256':sha(repo/'examples/finite_density/massless_three_loop_chain.json')}}
(root/('build-'+attempt+'.json')).write_text(json.dumps(report,indent=2)+'\n')
t=time.monotonic()
with open(root/('build-'+attempt+'.log'),'w') as log:p=subprocess.run(cmd,cwd=repo,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'CARGO_CRATE_NAME':'partial_source_class_prototype','CARGO_MANIFEST_DIR':str(repo)})
report.update(exit_code=p.returncode,wall_seconds=time.monotonic()-t,max_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
if p.returncode==0:report['executable']={'path':cmd[-1],'sha256':sha(cmd[-1])}
(root/('build-'+attempt+'.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['exit_code','wall_seconds','max_rss_kib']}));sys.exit(p.returncode)
