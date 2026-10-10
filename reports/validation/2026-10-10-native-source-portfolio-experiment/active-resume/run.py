#!/usr/bin/env python3
"""Continue only a completed timeout, retaining cumulative rounds and search limits."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from checkpoint import validate, SETTINGS

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('engine',choices=['baseline','proposal'])
p.add_argument('parent',type=Path)
p.add_argument('tag')
a=p.parse_args()
assert a.tag and all(c.isalnum()or c=='-'for c in a.tag)
base=Path(__file__).resolve().parent
root=base.parents[3]
parent=a.parent.resolve()
out=base/a.tag
assert not out.exists()
read=lambda path:json.loads(path.read_text())
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
config=read(parent/'configuration.json')
prior=read(parent/'run-binding.json')
prior_resources=read(parent/'resources.json')
progress=read(parent/'progress.json')
assert progress
last=progress[-1]['round']
round_json=parent/f'round-{last:03}.json'
round_bin=parent/f'round-{last:03}.bin'
state=read(round_json)
next_round=validate(a.engine,config,prior,prior_resources,progress,state,
                    result_exists=(parent/'result.json').exists(),closed_exists=(parent/'closed.bin').exists())
assert round_bin.is_file()
build_path=base/(a.engine+'-build.json')
build=read(build_path)
assert build['exit_code']==0 and build['engine']==a.engine
exe=Path(build['executable']['path'])
assert sha(exe)==build['executable']['sha256']
assert sha(Path(build['source']['path']))==build['source']['sha256']
assert all(sha(Path(value['path']))==value['sha256']for value in build['dependencies'].values())
corpus=Path(config['source_program']).resolve()
assert sha(corpus)==prior['source_corpus_sha256']
assert prior_resources['command']==prior['command']
assert Path(prior['command'][3]).resolve()==corpus
if 'resume_build_binding_path'in prior:
    parent_build_path=Path(prior['resume_build_binding_path'])
    parent_build=read(parent_build_path)
else:
    parent_build_path=(base.parent/'active-pilot/build-binding.json'if a.engine=='proposal'else
        root/'reports/validation/2026-10-10-global-boost-source-experiment/native-active-pilot/pilot-build-binding.json')
    parent_build=read(parent_build_path)
assert sha(parent_build_path)==prior['native_binding_sha256']
assert parent_build['dependencies']['rustred']['sha256']==build['dependencies']['rustred']['sha256']
prior_executable=Path(prior['command'][2])
assert sha(prior_executable)==prior['executable_sha256']
prior_provenance=read(parent/'resources.provenance.json')
assert prior_provenance['command']==prior['command']
assert prior_provenance['argument_files'][str(prior_executable)]['sha256']==prior['executable_sha256']
assert prior_provenance['argument_files'][str(corpus)]['sha256']==prior['source_corpus_sha256']
bound=[parent/name for name in ['configuration.json','run-binding.json','resources.json',
    'resources.provenance.json','progress.json']]+[round_json,round_bin,corpus,build_path,
    parent_build_path,exe,Path(build['source']['path']),Path(__file__),base/'checkpoint.py']
hashes={str(path):sha(path)for path in bound}
prior_wall=prior.get('cumulative_wall_seconds_prior',0)+prior_resources['wall_seconds']
command=['timeout','600',str(exe),str(corpus),str(out),str(parent)]
record={'scope':'Bounded continuation from the last committed completed round of a timed-out run. Original global limits retained; unfinished work is repeated and its prior cost remains included. No physical numerical acceptance.',
    'engine':a.engine,'command':command,'settings':SETTINGS,'source_corpus_sha256':sha(corpus),
    'executable_sha256':sha(exe),'native_binding_sha256':sha(build_path),
    'resume_build_binding_path':str(build_path),'build_binding_sha256':sha(build_path),
    'parent':str(parent),'first_continuation_round':next_round,
    'cumulative_wall_seconds_prior':prior_wall,'continuation_timeout_seconds':600,
    'frozen_parent_and_build_artifact_sha256':hashes,
    'discarded_uncommitted_round_files':[str(path)for path in sorted(parent.glob('round-*'))
        if path.suffix in ['.json','.bin']and int(path.stem.split('-')[-1])>last],
    'original_corpus_rules_imported':False,'parent_discovered_rules_restored_with_native_replay':True}
out.mkdir()
# Preserve the commit marker and its last complete pair in the continuation
# directory too. A timeout before its first new round can then still resume
# this same authenticated state, without inventing progress or resetting caps.
for path in [parent/'configuration.json',parent/'progress.json',round_json,round_bin]:
    shutil.copyfile(path,out/path.name)
binding_path=out/'run-binding.json'
binding_path.write_text(json.dumps(record,indent=2)+'\n')
result=subprocess.run([sys.executable,str(root/'reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py'),str(out/'resources.json'),*command],
                      cwd=root,env={**os.environ,**SETTINGS})
assert all(sha(Path(path))==value for path,value in hashes.items())
record['exit_code']=result.returncode
record['post_run_inputs_unchanged']=True
record['cumulative_wall_seconds']=prior_wall+read(out/'resources.json')['wall_seconds']
binding_path.write_text(json.dumps(record,indent=2)+'\n')
raise SystemExit(result.returncode)
