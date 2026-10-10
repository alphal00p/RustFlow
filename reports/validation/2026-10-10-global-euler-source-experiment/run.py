#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
mode=sys.argv[1];assert mode in ['original','global-euler']
build=json.loads((BASE/'build-binding.json').read_text())
executable=Path(build['executable']['path'])
assert digest(executable)==build['executable']['sha256']
assert digest(BASE/'probe.rs')==build['source']['sha256']
source=ROOT/'reports/validation/2026-10-10-targeted-native-depth-controls/double/input.bin'
points=BASE/'points.json';physical=ROOT/'examples/finite_density/massless_three_loop_chain.json'
out=BASE/mode;assert not out.exists()
command=['timeout','300',str(executable),str(source),str(points),str(out),mode,'3',str(physical)]
binding=BASE/f'{mode}-binding.json';assert not binding.exists()
inputs={str(p):digest(p) for p in [source,points,physical,executable,BASE/'probe.rs']}
report={'scope':'Fresh source-only global Euler presentation control, historical exact replay solely checks actual NoRule selection; no old rules imported into discovery; no periods or closure claim',
        'command':command,'inputs_sha256':inputs,'build_binding_sha256':digest(BASE/'build-binding.json'),'launcher_sha256':digest(Path(__file__))}
binding.write_text(json.dumps(report,indent=2)+'\n')
r=subprocess.run([sys.executable,str(ROOT/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(BASE/f'{mode}-resources.json'),*command],cwd=ROOT)
report.update(exit_code=r.returncode,inputs_unchanged=all(digest(Path(p))==h for p,h in inputs.items()))
assert report['inputs_unchanged'];binding.write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(r.returncode)
