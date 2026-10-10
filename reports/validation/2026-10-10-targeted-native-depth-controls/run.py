#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
mode=sys.argv[1]; depth=int(sys.argv[2])
assert mode in ['single3','double'] and depth in [3,4,5]
source=BASE/mode/'input.bin'; points=BASE/mode/'points.json'
executable=Path('/tmp/rustflow-targeted-depth-probe'); out=BASE/mode/f'depth-{depth}'
assert not out.exists()
build=json.loads((BASE/'build-binding.json').read_text())
assert digest(executable)==build['executable']['sha256']
command=['timeout','300',str(executable),str(source),str(points),str(out),'original',str(depth)]
inputs={str(p):digest(p)for p in [executable,source,points]}
report={'scope':'Targeted depth control; historical replay ONLY verifies actual NoRule selection, then original source definitions alone drive new ray1/point1 searches. No periods or old-rule union.',
        'inputs':inputs,'command':command,'build_binding_sha256':digest(BASE/'build-binding.json'),'launcher_sha256':digest(Path(__file__))}
binding=BASE/mode/f'depth-{depth}-binding.json'; assert not binding.exists()
binding.write_text(json.dumps(report,indent=2)+'\n')
result=subprocess.run([sys.executable,str(ROOT/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(BASE/mode/f'depth-{depth}-resources.json'),*command],cwd=ROOT)
report['exit_code']=result.returncode;report['inputs_unchanged']=all(digest(Path(p))==h for p,h in inputs.items())
binding.write_text(json.dumps(report,indent=2)+'\n');assert report['inputs_unchanged'];raise SystemExit(result.returncode)
