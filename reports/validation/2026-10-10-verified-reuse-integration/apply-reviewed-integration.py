#!/usr/bin/env python3
"""Apply the reviewed reuse changes only after the frozen baseline has ended."""
import hashlib,json,subprocess
from pathlib import Path
ROOT=Path.cwd(); BASE=Path(__file__).resolve().parent
PROTOTYPE=ROOT/'reports/validation/2026-10-10-verified-program-reuse'
NATIVE=ROOT/'vendor/rustred/crates/rustred-core'
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())

def main():
    record=BASE/'applied-source-binding.json'
    assert not record.exists(),'preserve the prior source application'
    baseline=ROOT/'reports/validation/2026-10-10-longer-singleton-control/run-binding.json'
    done=read(baseline)
    assert 'exit_code' in done and done['post_run_source_and_executable_unchanged'], 'baseline is still active or its freeze failed'
    build=read(PROTOTYPE/'library-build.json'); tests=read(PROTOTYPE/'native-tests-resources.json')
    assert build['exit_code']==tests['exit_code']==0
    assert build['production_sources_unchanged'] and build['isolated_sources_unchanged']
    assert all(h(NATIVE/name)==sha for name,sha in build['base_sources'].items())
    native_patch=PROTOTYPE/'isolated-native.patch'
    assert h(native_patch)==build['patch_sha256']
    draft=BASE/'integration-draft'; mapping=read(draft/'source-map.json')
    assert all(h(ROOT/row['path'])==row['base_sha256'] for row in mapping['files'])
    audit=mapping['unchanged_final_audit']; assert h(ROOT/audit['path'])==audit['sha256']
    adapter_patch=draft/'verified-reuse-integration.patch';assert h(adapter_patch)==mapping['patch_sha256']
    actions=[['git','apply','--directory=vendor/rustred/crates/rustred-core',str(native_patch)],['git','apply',str(adapter_patch)]]
    for action in actions:subprocess.run(action[:2]+['--check']+action[2:],check=True)
    for action in actions:subprocess.run(action,check=True)
    paths=[NATIVE/'src/solver/guarded'/n for n in ('lifecycle.rs','persistence.rs','verified_tests.rs')]+[ROOT/row['path'] for row in mapping['files']]
    before_format={str(p.relative_to(ROOT)):h(p) for p in paths}
    subprocess.run(['rustfmt','--edition=2024','--config','skip_children=true',*[str(p) for p in paths]],check=True)
    assert h(ROOT/audit['path'])==audit['sha256']
    subprocess.run(['git','diff','--check','--',*[str(p) for p in paths]],check=True)
    record.write_text(json.dumps({'scope':'Reviewed immutable verified-program reuse only. Existing sources, physical admission, discovery and final closure audit unchanged.',
        'completed_baseline_sha256':h(baseline),'baseline_exit_code':done['exit_code'],
        'prototype_manifest_sha256':h(PROTOTYPE/'artifact-manifest.json'),
        'native_patch_sha256':h(native_patch),'adapter_patch_sha256':h(adapter_patch),
        'before_format_sha256':before_format,'production_source_sha256':{str(p.relative_to(ROOT)):h(p) for p in paths},
        'unchanged_final_audit':audit,'coherent_build_pending':True,'launcher_sha256':h(Path(__file__))},indent=2)+'\n')
    print('Applied reviewed sources; final source snapshot and coherent release gates remain required.')
if __name__=='__main__':main()
