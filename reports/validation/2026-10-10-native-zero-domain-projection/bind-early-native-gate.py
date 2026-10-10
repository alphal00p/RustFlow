"""Bind the already executed fresh native suite to the completed shared build.

This does not execute or repeat any test. Every executable/source binding must
match the pre-launch records of the early native run.
"""
import argparse,json
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix',required=True)
    args=parser.parse_args()
    prefix=args.prefix
    if not prefix or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789-' for c in prefix):raise ValueError('unsafe prefix')
    snapshot=BASE/(prefix+'-source-hashes.json');verify_snapshot(snapshot)
    build_path=BASE/(prefix+'-build-provenance.json');build=json.loads(build_path.read_text())
    early_path=BASE/(prefix+'-native-early-binding.json');early=json.loads(early_path.read_text())
    resource_path=BASE/(prefix+'-native-resources.json');resources=json.loads(resource_path.read_text())
    launch_path=resource_path.with_suffix('.provenance.json');launch=json.loads(launch_path.read_text())
    name=early['executable']['path'];expected=early['executable']['sha256']
    assert early['source_snapshot_sha256']==build['source_snapshot_sha256']==digest(snapshot)
    assert expected==build['files'][name]['sha256']==digest(ROOT/name)==launch['argument_files'][name]['sha256']
    assert resources['command']==launch['command']
    assert resources['exit_code']==0
    assert resources['command']==['timeout','180',name,'solver::guarded','--nocapture','--test-threads=1']
    output=BASE/(prefix+'-native-binding.json')
    if output.exists():raise ValueError('preserve existing gate binding')
    report={'source_snapshot_sha256':digest(snapshot),'build_provenance_sha256':digest(build_path),
        'launcher_sha256':digest(Path(__file__)),'gate':'native','command':resources['command'],
        'exit_code':resources['exit_code'],'native_numerical_acceptance_run':False,
        'scope':'The fresh native suite ran once after its executable linked; this record binds that exact execution to the subsequently completed shared build.',
        'early_binding':{'path':str(early_path.relative_to(ROOT)),'sha256':digest(early_path)},
        'resources_sha256':digest(resource_path),'launch_sha256':digest(launch_path)}
    output.write_text(json.dumps(report,indent=2)+'\n')
    print('Already executed native gate is bound to the completed shared build; no rerun.')

if __name__=='__main__':main()
