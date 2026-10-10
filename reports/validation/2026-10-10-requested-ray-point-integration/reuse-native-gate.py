#!/usr/bin/env python3
"""Reuse prior native assertions only when source, features and binary hashes match."""
import argparse
import json
import re
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix',default='requested-ray-point')
    args=parser.parse_args()
    assert re.fullmatch(r'[a-z0-9-]+',args.prefix)
    output=BASE/(args.prefix+'-native-gate-reuse.json')
    assert not output.exists()
    prior=ROOT/'reports/validation/2026-10-10-native-zero-domain-projection'
    old_snapshot_path=prior/'guard-point-source-hashes.json'
    new_snapshot_path=BASE/(args.prefix+'-source-hashes.json')
    old=json.loads(old_snapshot_path.read_text())
    new=verify_snapshot(new_snapshot_path)
    native=lambda snap:{key:value for key,value in snap['files'].items()
        if key.startswith('vendor/rustred/')or key in ['Cargo.toml','Cargo.lock']}
    assert native(old)==native(new),'native source or dependency declarations changed; rerun native gate'
    old_build_path=prior/'guard-point-build-provenance.json'
    new_build_path=BASE/(args.prefix+'-build-provenance.json')
    old_build=json.loads(old_build_path.read_text());new_build=json.loads(new_build_path.read_text())
    assert new_build['source_snapshot_sha256']==digest(new_snapshot_path)
    for key in ['rustred_unit','rustred_library']:
        assert new_build['features'][key]==old_build['features'][key]
    binaries=[key for key in old_build['files']if Path(key).name.startswith(('rustred-','librustred-'))]
    assert len(binaries)==2
    for path in binaries:
        assert new_build['files'][path]['sha256']==old_build['files'][path]['sha256']==digest(ROOT/path)
    resources=prior/'guard-point-native-resources.json';log=prior/'guard-point-native-resources.log'
    binding=prior/'guard-point-native-binding.json';record=json.loads(binding.read_text())
    assert record['exit_code']==json.loads(resources.read_text())['exit_code']==0
    assert record['source_snapshot_sha256']==digest(old_snapshot_path)
    assert record['build_provenance_sha256']==digest(old_build_path)
    results=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',log.read_text())
    assert len(results)==1 and results[0][1]=='0'
    verify_snapshot(new_snapshot_path)
    report={'status':'prior_native_gate_reused_by_exact_identity','scope':'No native test rerun. All native sources, dependency declarations, native features, test executable and library hashes match the prior passing gate.',
        'reused_passed_tests':int(results[0][0]),'reused_ignored_tests':int(results[0][2]),
        'new_native_test_invocations':0,'new_numerical_runs':0,
        'native_source_files':len(native(new)),
        'native_binaries':{path:new_build['files'][path]['sha256']for path in binaries},
        'artifacts':{str(path.relative_to(ROOT)):digest(path)for path in
            [old_snapshot_path,new_snapshot_path,old_build_path,new_build_path,resources,log,binding,Path(__file__)]}}
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'reused_passed_tests':report['reused_passed_tests']}))

if __name__=='__main__':main()
