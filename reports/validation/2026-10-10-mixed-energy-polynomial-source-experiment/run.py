#!/usr/bin/env python3
"""Run one source-only control against the frozen native library."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

mode = sys.argv[1]
assert mode in ['slot6', 'slot8', 'both', 'both-reverse', 'temporal-only', 'boost-only']
source = ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration/requested-ray-point-three-loop-fixed/native-closure/cut-3/round-000-provisional.bin'
executable = Path('/tmp/rustflow-mixed-energy-polynomial-probe')
points = BASE / 'points.json'
physical = ROOT / 'examples/finite_density/massless_three_loop_chain.json'
out = BASE / mode
assert not out.exists()
build = json.loads((BASE / 'build-binding.json').read_text())
assert digest(executable) == build['executable']['sha256']
assert digest(BASE / 'probe.rs') == build['source']['sha256']
command = ['timeout', '300', str(executable), str(source), str(points), str(physical), str(out), mode]
inputs = {str(path): digest(path) for path in [executable, source, points, physical]}
report = {'scope': 'Fresh native exact point discovery from complete original cut3 source metadata plus exact polynomial whole-row shifts; no historical rule imported, no numerical values.',
          'command': command, 'inputs': inputs,
          'launcher_sha256': digest(Path(__file__)),
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
