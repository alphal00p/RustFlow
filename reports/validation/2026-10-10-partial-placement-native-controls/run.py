#!/usr/bin/env python3
import hashlib
import json
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


mode = sys.argv[1]
assert mode in ['single-all', 'single-slot0', 'double-all', 'double-slots12', 'double-slots24']
build_path = BASE / 'build-binding.json'
build = json.loads(build_path.read_text())
assert build['exit_code'] == 0
exe = Path(build['executable']['path'])
assert sha(exe) == build['executable']['sha256']
assert sha(BASE / 'probe.rs') == build['source']['sha256']
for dependency in build['dependencies'].values():
    assert sha(Path(dependency['path'])) == dependency['sha256']
proof_paths = [BASE / 'chart_check.py', BASE / 'chart-checks.json', BASE / 'origin-proof.md']
physical = ROOT / 'examples/finite_density/massless_three_loop_chain.json'
proof = {str(p.relative_to(ROOT)): sha(p) for p in proof_paths}
proof_digest = hashlib.sha256(json.dumps(proof, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
out = BASE / mode
assert not out.exists()
binding = BASE / f'{mode}-binding.json'
assert not binding.exists()
command = ['timeout', '300', str(exe), str(physical), str(out), mode, proof_digest]
inputs = {str(p): sha(p) for p in [physical, exe, BASE / 'probe.rs', *proof_paths]}
report = {'scope': 'Fresh source-only restricted partial-placement comparison; no sealed production permit, old proof reuse, numerical values, or closure claim',
          'mode': mode, 'command': command, 'inputs_sha256': inputs, 'proof_inputs_sha256': proof,
          'proof_binding_sha256': proof_digest, 'build_binding_sha256': sha(build_path),
          'launcher_sha256': sha(Path(__file__)),
          'independent_proof_audit': 'acceptance_manifest read-only review: no blocker for full-positive-support lower-origin boxes; excludes deleted supports and eta endpoint'}
binding.write_text(json.dumps(report, indent=2) + '\n')
wrapper = ROOT / 'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'
result = subprocess.run([sys.executable, str(wrapper), str(BASE / f'{mode}-resources.json'), *command], cwd=ROOT)
report.update(exit_code=result.returncode, inputs_unchanged=all(sha(Path(p)) == h for p, h in inputs.items()))
assert report['inputs_unchanged']
binding.write_text(json.dumps(report, indent=2) + '\n')
raise SystemExit(result.returncode)
