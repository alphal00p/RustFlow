"""Freeze a completed unsuccessful baseline without promoting its frontier."""
import gzip,hashlib,json,subprocess,re
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
assert not (HERE/'artifact-manifest.json').exists()
binding=read(HERE/'run-binding.json');resources=read(HERE/'resources.json')
assert binding['exit_code']==resources['exit_code'] and binding['exit_code']!=0
assert binding['post_run_source_and_executable_unchanged']
folder=HERE/'single3';corpus=folder/'native-closure'
assert not list(folder.glob('prediction-*.json')) and not list(corpus.glob('*-closed.json'))
provisionals=sorted(corpus.glob('round-*-provisional.json'));assert provisionals
last=provisionals[-1];meta=read(last);binary=last.with_suffix('.bin')
K=HERE.parent/'2026-10-10-requested-ray-point-integration'
toolrecord=K/'native-digest-build.json';tool=read(toolrecord);toolpath=Path(tool['executable']['path']);assert h(toolpath)==tool['executable']['sha256']
b3,n=subprocess.check_output([str(toolpath),str(binary)],text=True).split();assert b3==meta['program_blake3'] and int(n)==meta['program_bytes']==binary.stat().st_size
link=meta['active_state']['requested_discovery'];assert re.fullmatch('requested-discovery-[0-9]+[.]json',link['path'])
txpath=corpus/link['path'];tx=read(txpath)
assert subprocess.check_output([str(toolpath),str(txpath)],text=True).split()[0]==link['blake3']
assert tx['native_program']=={'blake3':b3,'bytes':int(n)} and tx['status']=='complete'
summary={'status':'timed_out_before_closed_checkpoint'if binding['exit_code']==124 else'failed_before_closed_checkpoint',
 'exit_code':binding['exit_code'],'wall_cap_seconds':3600,'completed_provisional_rounds':len(provisionals),
 'frontier_history':[len(read(p)['frontier'])for p in provisionals],
 'last_provisional_round':meta['round'],'last_provisional_frontier_size':len(meta['frontier']),
 'native_rules':meta['native_rule_count'],'historical_requests':len(meta['active_state']['historical_requests']),
 'last_submitted_requests':len(meta['active_state']['submitted_requests']),
 'unsubmitted_requests':len(meta['active_state']['deferred_requests']),
 'closed_system_claim':False,'basis_claim':False,'predictions':0,'full_three_loop_acceptance':False,
 'resources':resources,'scope':'Fresh unchanged production baseline at enlarged bounds. The last frontier is provisional; lack of new discovery requests does not certify completed final replay, target/derivative audit, historical candidate collection, boundaries or numerical evaluation. Post-frontier timing attribution was not instrumented in this executable. Separate profiling is independent. All earlier failed costs are retained; no resume or speed claim.',
 'bindings_sha256':{str(p.relative_to(ROOT)):h(p)for p in [HERE/'run.py',HERE/'run-binding.json',HERE/'resources.json',HERE/'resources.provenance.json',HERE/'resources.log',last,binary,txpath,toolrecord]}}
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(HERE/'README.md').write_text(f"The fresh unchanged baseline ended with exit {binding['exit_code']} at its 3,600-second cap. It saved {len(provisionals)} provisional rounds; the frontier contracted from {max(summary['frontier_history'])} to {len(meta['frontier'])} labels, with {meta['native_rule_count']} rules and {summary['historical_requests']} historical requests. No final closed checkpoint or prediction was produced. The last frontier is not an accepted basis.\n\nThe executable did not separately time its final source replay, original target/basis-derivative audit, and exhaustive historical candidate-map collection. A separate diagnostic investigates those phases without turning this timed-out run into a closure claim. The one-hour cost remains explicit. The post-run source and executable checks passed.\n\nLatest provisional proof/discovery metadata and the initial proof remain raw; completed older payloads and the log are compressed losslessly with byte/hash restoration checks. There is no checkpoint-resume API and no continuation claim.\n")
preserve={last,binary,txpath,provisionals[0],provisionals[0].with_suffix('.bin')}
archives=[]
for p in sorted([p for p in corpus.rglob('*')if p.is_file()and p not in preserve]+[HERE/'resources.log']):
 raw=p.read_bytes();dst=p.with_name(p.name+'.gz');assert not dst.exists();dst.write_bytes(gzip.compress(raw,compresslevel=9,mtime=0));assert gzip.decompress(dst.read_bytes())==raw
 archives.append({'original_path':str(p.relative_to(HERE)),'archive_path':str(dst.relative_to(HERE)),'original_sha256':hashlib.sha256(raw).hexdigest(),'archive_sha256':h(dst),'original_bytes':len(raw),'archive_bytes':dst.stat().st_size});p.unlink()
(HERE/'archive-map.json').write_text(json.dumps({'gzip_mtime':0,'restoration_checked':True,'entries':archives,'preserved_raw':[str(p.relative_to(HERE))for p in sorted(preserve)]},indent=2)+'\n')
mapping={row['original_path']:row for row in archives}
for name,expected in summary['bindings_sha256'].items():
 p=ROOT/name;raw=p.read_bytes()if p.exists()else gzip.decompress((HERE/mapping[str(p.relative_to(HERE))]['archive_path']).read_bytes());assert hashlib.sha256(raw).hexdigest()==expected
files={str(p.relative_to(HERE)):{'sha256':h(p),'bytes':p.stat().st_size}for p in sorted(HERE.rglob('*'))if p.is_file()and'__pycache__'not in p.parts}
(HERE/'artifact-manifest.json').write_text(json.dumps({'scope':summary['scope'],'files':files,'file_count':len(files),'archive_and_binding_checks':'passed'},indent=2)+'\n')
print(json.dumps({'status':summary['status'],'files':len(files),'manifest_sha256':h(HERE/'artifact-manifest.json')}))
