#!/usr/bin/env python3
"""Freeze completed fresh runs and move the large executable outside the repo."""
import hashlib
import json
from pathlib import Path
import shutil

BASE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for mode in ['double-slots12','double-slots24']:
    binding=json.loads((BASE/f'{mode}-binding.json').read_text())
    resource=json.loads((BASE/f'{mode}-resources.json').read_text())
    assert binding['capsule_unchanged'] and resource['exit_code']==binding['exit_code']
    progress=json.loads((BASE/mode/'progress.json').read_text())
    result_path=BASE/mode/'result.json'
    result=json.loads(result_path.read_text()) if result_path.exists() else None
    status=result['status'] if result is not None else 'wall_timeout_unclosed'
    if result is None: assert resource['exit_code']==124
    rows.append({'mode':mode,'status':status,'resources':resource,
                 'completed_rounds':len(progress),'frontier_history':[r['frontier'] for r in progress],
                 'last_completed_round':progress[-1],
                 'closed':result.get('closed') if result else None,
                 'binding_sha256':sha(BASE/f'{mode}-binding.json')})
capsule=json.loads((BASE/'capsule-manifest.json').read_text())
source=BASE/'capsule/pilot'
item=next(i for i in capsule['files'] if i['capsule']=='capsule/pilot')
assert source.stat().st_size==item['bytes'] and sha(source)==item['sha256']
destination=Path('/tmp/rustflow-partial-placement-longer-double-20261010/pilot')
destination.parent.mkdir(parents=True,exist_ok=True)
assert not destination.exists()
shutil.move(source,destination)
assert destination.stat().st_size==item['bytes'] and sha(destination)==item['sha256']
relocation={'scope':'Byte-identical static executable relocated only after both process and post-run capsule checks completed. Original runtime manifests are unchanged.',
            'original_path':str(source),'relative_original_path':'capsule/pilot',
            'current_path':str(destination),'bytes':item['bytes'],'sha256':item['sha256']}
(BASE/'executable-relocation.json').write_text(json.dumps(relocation,indent=2)+'\n')
summary={'scope':'Fresh larger-bound algebraic closure controls only; no old rule import, boundary values, transport, references or physical admission.',
         'bounds':{'wall_seconds':1200,'rounds':24,'requests':16384,'frontier':4096,'rules':65536,
                   'ray_domains':1,'point_domains':1,'depth':3,'zero_attempts_per_round':8192},
         'unchanged_source_and_algorithm':True,'arms':rows,
         'first_launcher_error':'One shell invocation used a mistyped Python path and exited127 before any native process or output directory. The successful launch used the recorded installed interpreter.',
         'runtime_library_files_required':False}
(BASE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
files={str(p.relative_to(BASE)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(BASE.rglob('*')) if p.is_file() and p.name!='artifact-manifest.json'}
manifest={'files':files,'file_count':len(files),'total_bytes':sum(i['bytes'] for i in files.values()),'external_byte_identical_executable':relocation}
(BASE/'artifact-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'arms':[{'mode':r['mode'],'status':r['status'],'rounds':r['completed_rounds']} for r in rows],
                  'files':len(files),'manifest_sha256':sha(BASE/'artifact-manifest.json')}))
