#!/usr/bin/env python3
"""Use externally bound baseline selection checks; never import foreign proofs."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

PARENT = Path(__file__).resolve().parent
HERE = PARENT/'source-only'
ROOT = PARENT.parents[2]
DEPTH = ROOT/'reports/validation/2026-10-10-targeted-native-depth-controls'
K = ROOT/'reports/validation/2026-10-10-requested-ray-point-integration'
sys.path.insert(0,str(K))
from frozen_sources import digest, verify_snapshot


def meta(path):
    return {'path':str(path),'sha256':digest(path)}


def main():
    HERE.mkdir()
    prior = json.loads((PARENT/'build-binding.json').read_text())
    assert prior['exit_code'] == 0
    assert digest(PARENT/'probe.rs') == prior['source']['sha256']
    source = (PARENT/'probe.rs').read_text()
    start = source.index('    let historical=GuardedProgram::decode_generated')
    end = source.index('    let out=std::path::PathBuf', start)
    source = source[:start] + '''    // Historical status is independently proved by the bound baseline run.
    // Never import its proof into the isolated portfolio schema.
    let baseline:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[6]).unwrap()).unwrap();
    assert_eq!(baseline["depth"],3);assert_eq!(baseline["ray_domains"],1);assert_eq!(baseline["point_domains"],1);
    assert_eq!(baseline["source_identity"],record.measure);
    let observations=baseline["selection_observations"].as_array().unwrap().clone();
    assert_eq!(observations.len(),selected.len());
    for (point,check) in selected.iter().zip(&observations) {
        assert_eq!(point["point"],check["point"]);assert_eq!(point["kind"],check["kind"]);
        if point["kind"]=="actual-frontier" { assert_eq!(check["status"],"Unresolved(NoApplicableRule)"); }
        else { assert!(check["status"].as_str().unwrap().starts_with("Applied")); }
    }
''' + source[end:]
    source = source.replace('Historical program is replayed ONLY to verify selection status, then discarded before fresh discovery;',
                            'Historical selection was verified by the independently bound baseline run; no foreign proof is imported;')
    (HERE/'probe.rs').write_text(source)
    for sector in ['single3','double']:
        case=HERE/sector
        case.mkdir()
        baseline_binding=json.loads((DEPTH/sector/'depth-3-binding.json').read_text())
        assert baseline_binding['exit_code']==0 and baseline_binding['inputs_unchanged']
        assert baseline_binding['build_binding_sha256']==digest(DEPTH/'build-binding.json')
        for path, expected in baseline_binding['inputs'].items():
            assert digest(Path(path))==expected
        for name in ['input.bin','input.json','points.json']:
            shutil.copyfile(PARENT/sector/name,case/name)
            assert digest(case/name)==digest(DEPTH/sector/name)
        shutil.copyfile(DEPTH/sector/'depth-3/result.json',case/'baseline-result.json')
    executable=Path('/tmp/rustflow-targeted-portfolio-source-only-20261010')
    assert not executable.exists()
    command=list(prior['command'])
    command[command.index(str(PARENT/'probe.rs'))]=str(HERE/'probe.rs')
    command[-1]=str(executable)
    record={'scope':'Source-only isolated native portfolio, with independently bound production-native historical selection verification. Foreign proof schemas remain rejected.',
            'prior_failed_attempts':{sector:meta(PARENT/sector/'run-binding.json') for sector in ['single3','double']},
            'source':meta(HERE/'probe.rs'),'parent_source':meta(PARENT/'probe.rs'),
            'launcher':meta(Path(__file__)),'dependencies':prior['dependencies'],'command':command,'cargo_invoked':False}
    for item in record['dependencies'].values():
        assert digest(Path(item['path']))==item['sha256']
    path=HERE/'build-binding.json'
    path.write_text(json.dumps(record,indent=2)+'\n')
    started=time.monotonic()
    with (HERE/'build.log').open('w') as log:
        result=subprocess.run(command,cwd=ROOT,env={**os.environ,'CARGO_CRATE_NAME':'targeted_portfolio_probe'},stdout=log,stderr=subprocess.STDOUT)
    record.update(exit_code=result.returncode,wall_seconds=time.monotonic()-started)
    if result.returncode==0:record['executable']=meta(executable)
    path.write_text(json.dumps(record,indent=2)+'\n')
    if result.returncode:raise SystemExit(result.returncode)
    for sector in ['single3','double']:
        case=HERE/sector
        command=['timeout','300',str(executable),str(case/'input.bin'),str(case/'points.json'),
                 str(case/'proposal-depth-3'),'original','3',str(case/'baseline-result.json')]
        inputs=[case/name for name in ['input.bin','input.json','points.json','baseline-result.json']]
        inputs += [DEPTH/sector/'depth-3-binding.json',DEPTH/'build-binding.json',DEPTH/'manifest.json']
        binding={'engine':'isolated-source-portfolio','scope':record['scope'],'command':command,
                 'build_binding':meta(path),'executable':meta(executable),
                 'inputs':{str(p):digest(p) for p in inputs}}
        run_path=case/'run-binding.json'
        run_path.write_text(json.dumps(binding,indent=2)+'\n')
        result=subprocess.run([sys.executable,str(K/'run-resource-command.py'),str(case/'resources.json'),*command],cwd=ROOT)
        assert all(digest(Path(p))==h for p,h in binding['inputs'].items())
        binding.update(exit_code=result.returncode,inputs_unchanged=True)
        run_path.write_text(json.dumps(binding,indent=2)+'\n')
        print(json.dumps({'sector':sector,'exit_code':result.returncode}),flush=True)
    verify_snapshot(K/'requested-ray-point-source-hashes.json')


if __name__=='__main__':main()
