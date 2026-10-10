#!/usr/bin/env python3
"""Matched small searches against the previously tested isolated native portfolio."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
DEPTH = ROOT / 'reports/validation/2026-10-10-targeted-native-depth-controls'
PORTFOLIO = ROOT / 'reports/validation/2026-10-10-native-source-portfolio-experiment'
K = ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration'
sys.path.insert(0, str(K))
from frozen_sources import digest, verify_snapshot


def meta(path):
    return {'path': str(path), 'sha256': digest(path)}


def main():
    assert not (HERE/'build-binding.json').exists()
    depth_build = json.loads((DEPTH/'build-binding.json').read_text())
    assert depth_build['exit_code'] == 0
    assert digest(DEPTH/'probe.rs') == depth_build['source']['sha256']
    predecessor = PORTFOLIO/'proposal-generic-probe-build.json'
    old = json.loads(predecessor.read_text())
    assert old['exit_code'] == 0
    dependencies = {name:Path(item['path']) for name,item in old['dependencies'].items()}
    for name,path in dependencies.items():
        assert digest(path) == old['dependencies'][name]['sha256']
    source = (DEPTH/'probe.rs').read_text()
    needle = '"proof_sources":program.rules()'
    assert source.count(needle) == 1
    source = source.replace(needle, '"native_portfolio":program.rules().iter().map(experiment_stats).collect::<Vec<_>>(),'+needle)
    (HERE/'probe.rs').write_text(source)
    for sector in ['single3', 'double']:
        (HERE/sector).mkdir()
        for name in ['input.bin', 'input.json', 'points.json']:
            shutil.copyfile(DEPTH/sector/name, HERE/sector/name)
    executable = Path('/tmp/rustflow-targeted-portfolio-probe-20261010')
    assert not executable.exists()
    command = ['rustc', '--edition=2024', '--crate-name', 'targeted_portfolio_probe',
               '-C', 'opt-level=0', '-C', 'debuginfo=0', '--cfg', 'source_portfolio']
    for name,path in dependencies.items():
        command += ['--extern', name+'='+str(path)]
    command += ['-L', 'dependency='+str(ROOT/'target/release/deps')]
    for name in ['gmp-mpfr-sys-2e6d668657be102d','gmp-mpfr-sys-c789de89ed67d3db']:
        command += ['-L', 'native='+str(ROOT/'target/release/build'/name/'out/lib')]
    command += [str(HERE/'probe.rs'), '-o', str(executable)]
    build = {
        'scope': 'Validation-only copied native search driver, using the frozen isolated source-portfolio RustRed library. Production sources, build and reducer unchanged; no speed comparison.',
        'parent_probe_build': meta(DEPTH/'build-binding.json'),
        'parent_probe_source': meta(DEPTH/'probe.rs'),
        'parent_portfolio_build': meta(predecessor),
        'source_change': 'Only expose already-existing native portfolio statistics in per-point diagnostic output.',
        'source': meta(HERE/'probe.rs'), 'launcher': meta(Path(__file__)),
        'dependencies': {name:meta(path) for name,path in dependencies.items()},
        'command': command, 'cargo_invoked': False,
    }
    (HERE/'build-binding.json').write_text(json.dumps(build,indent=2)+'\n')
    started = time.monotonic()
    with (HERE/'build.log').open('w') as log:
        result = subprocess.run(command,cwd=ROOT,env={**os.environ,'CARGO_CRATE_NAME':'targeted_portfolio_probe'},stdout=log,stderr=subprocess.STDOUT)
    build.update(exit_code=result.returncode, wall_seconds=time.monotonic()-started)
    if result.returncode == 0:
        build['executable'] = meta(executable)
    (HERE/'build-binding.json').write_text(json.dumps(build,indent=2)+'\n')
    if result.returncode:
        raise SystemExit(result.returncode)
    for sector in ['single3','double']:
        case = HERE/sector
        command = ['timeout','300',str(executable),str(case/'input.bin'),str(case/'points.json'),
                   str(case/'proposal-depth-3'),'original','3']
        binding = {'scope':'Matched depth3/ray1/point1 control with isolated source portfolio; original complete source corpus, domains and zero boxes preserved. Historical native program is replayed only to verify selected statuses and discarded before fresh discovery.',
                   'engine':'isolated-source-portfolio', 'command':command,
                   'build_binding':meta(HERE/'build-binding.json'),
                   'executable':meta(executable),
                   'inputs':{name:meta(case/name) for name in ['input.bin','input.json','points.json']}}
        path = case/'run-binding.json'
        path.write_text(json.dumps(binding,indent=2)+'\n')
        result = subprocess.run([sys.executable,str(K/'run-resource-command.py'),str(case/'resources.json'),*command],cwd=ROOT)
        for item in binding['inputs'].values():
            assert digest(Path(item['path'])) == item['sha256']
        assert meta(executable) == build['executable']
        binding.update(exit_code=result.returncode, inputs_unchanged=True)
        path.write_text(json.dumps(binding,indent=2)+'\n')
        print(json.dumps({'sector':sector,'exit_code':result.returncode}),flush=True)
    verify_snapshot(K/'requested-ray-point-source-hashes.json')


if __name__ == '__main__':
    main()
