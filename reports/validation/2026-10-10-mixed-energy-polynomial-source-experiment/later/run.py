#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]
BASE = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

mode = sys.argv[1]
assert mode in ['original', 'slot6', 'slot8', 'both', 'both-reverse', 'temporal-only', 'boost-only']
source = BASE.parent / mode / 'program.bin'
executable = Path('/tmp/rustflow-mixed-energy-later-probe')
points = BASE / 'points.json'
out = BASE / f'{mode}.json'
assert not out.exists()
build = json.loads((BASE / 'build-binding.json').read_text())
assert digest(executable) == build['executable']['sha256']
assert digest(BASE / 'probe.rs') == build['source']['sha256']
command = ['timeout', '300', str(executable), str(source), str(points), str(out), 'original']
inputs = {str(path): digest(path) for path in [executable, source, points]}
report = {'scope': '32 deterministic actual downstream cut3 labels, original source corpora plus polynomial shifts. Fresh native ray cap1/point cap1 discovery only; original program proofs discarded. No period values.',
          'command': command, 'mode': mode, 'inputs': inputs,
          'launcher_sha256': digest(Path(__file__)),
          'selection_sha256': digest(BASE / 'selection.json'),
          'build_binding_sha256': digest(BASE / 'build-binding.json')}
binding = BASE / f'{mode}-run-binding.json'
assert not binding.exists()
binding.write_text(json.dumps(report, indent=2) + '\n')
result = subprocess.run([sys.executable, str(ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),
                         str(BASE / f'{mode}-resources.json'), *command], cwd=ROOT)
report['exit_code'] = result.returncode
report['inputs_unchanged'] = all(digest(Path(path)) == value for path, value in inputs.items())
binding.write_text(json.dumps(report, indent=2) + '\n')
assert report['inputs_unchanged']
raise SystemExit(result.returncode)
