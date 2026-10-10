#!/usr/bin/env python3
"""Run a matched source-only point control against an explicitly hash-bound corpus."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('tag', choices=['baseline', 'proposal'])
p.add_argument('label')
p.add_argument('corpus', type=Path)
p.add_argument('corpus_sha256')
a = p.parse_args()
assert a.label and all(c.isalnum() or c in '-_' for c in a.label)
base = Path(__file__).resolve().parent
root = base.parents[2]
name = a.tag + '-' + a.label + '-19-points'
output = base / (name + '.json')
resources = base / (name + '-resources.json')
binding_path = base / (name + '-binding.json')
assert not any(x.exists() for x in [output, resources, binding_path])
build_path = base / (a.tag + '-generic-probe-build.json')
build = json.loads(build_path.read_text())
assert build['exit_code'] == 0
exe = Path(build['executable']['path'])
fixture = base / 'points.json'
corpus = a.corpus.resolve()
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
assert sha(exe) == build['executable']['sha256']
assert sha(corpus) == a.corpus_sha256
source = Path(build['source']['path'])
assert sha(source) == build['source']['sha256']
command = ['timeout', '600', str(exe), str(corpus), str(fixture), str(output), 'original']
binding = {'scope': 'Fresh exact native source discovery; explicit known v1/v2 source-only decoding. No imported historical rules, physical values or closure claim.',
           'build_binding_sha256': sha(build_path), 'executable_sha256': sha(exe),
           'source_corpus': str(corpus), 'source_corpus_sha256': sha(corpus),
           'fixture_sha256': sha(fixture), 'probe_source_sha256': sha(source),
           'launcher_sha256': sha(Path(__file__)), 'command': command}
binding_path.write_text(json.dumps(binding, indent=2) + '\n')
r = subprocess.run([sys.executable, str(root / 'reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py'), str(resources), *command], cwd=root)
assert sha(exe) == binding['executable_sha256'] and sha(corpus) == binding['source_corpus_sha256']
assert sha(fixture) == binding['fixture_sha256'] and sha(source) == binding['probe_source_sha256']
binding['exit_code'] = r.returncode
binding['post_run_inputs_unchanged'] = True
if output.exists(): binding['output_sha256'] = sha(output)
binding_path.write_text(json.dumps(binding, indent=2) + '\n')
raise SystemExit(r.returncode)
