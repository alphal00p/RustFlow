#!/usr/bin/env python3
"""Deterministic label-only sample and paired source-only native point discovery."""
import collections
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

base = Path(__file__).resolve().parent
root = base.parents[2]
portfolio = root / 'reports/validation/2026-10-10-native-source-portfolio-experiment'
frontier_path = portfolio / 'active-pilot/polynomial62-baseline/round-007.json'
build_path = portfolio / 'baseline-generic-probe-build.json'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
frontier = json.loads(frontier_path.read_text())['frontier']
assert len(frontier) == 290 and len({tuple(x) for x in frontier}) == 290
groups = collections.defaultdict(list)
for label in frontier:
    assert len(label) == 12 and label[11] == 0 and label[6] <= 0
    groups[(label[0], label[6])].append(label)
key = lambda x: (x[9], x[7], x[8], tuple(x))
selected, strata = [], []
for group, labels in sorted(groups.items()):
    ordered = sorted(labels, key=key)
    chosen = [ordered[0], ordered[-1]]
    assert chosen[0] != chosen[1]
    selected.extend(chosen)
    strata.append({'cut_power': group[0], 'occupied_energy_index': group[1],
                   'population': len(labels), 'selected': chosen})
largest = min(groups, key=lambda group: (-len(groups[group]), group))
ordered = sorted(groups[largest], key=key)
extra = ordered[(len(ordered)-1)//2]
assert extra not in selected
selected.append(extra)
assert len(selected) == 19 and len({tuple(x) for x in selected}) == 19
fixture = base / 'points.json'
fixture.write_text(json.dumps(selected, indent=2) + '\n')
(base / 'selection.json').write_text(json.dumps({
    'frontier_path': str(frontier_path), 'frontier_sha256': sha(frontier_path),
    'frontier_count': len(frontier), 'round': 7,
    'policy': 'Stratify by exact (required-cut index0, occupied-energy index6). Choose the first and last in each stratum sorted by (upper H9, virtual-energy indices7,8, full label), plus the lower median of the most populous stratum (ties by stratum). Selection uses labels only, no coefficients, rule outcome or reference values.',
    'strata': strata, 'extra_median_stratum': list(largest), 'extra_median_point': extra,
    'selected_count': len(selected), 'points_sha256': sha(fixture),
}, indent=2) + '\n')
build = json.loads(build_path.read_text())
exe, source = Path(build['executable']['path']), Path(build['source']['path'])
assert build['exit_code'] == 0
assert sha(exe) == build['executable']['sha256']
assert sha(source) == build['source']['sha256']
shutil.copy2(source, base / 'generic-source-probe.rs')
shutil.copy2(build_path, base / 'baseline-generic-probe-build.json')
corpora = {
    'polynomial62': root / 'reports/validation/2026-10-10-polynomial-raw-ward-attribution/native-source-probe/strongest62-boost-ward/program.bin',
    'reciprocal56': root / 'reports/validation/2026-10-10-normalized-boost-reciprocal-experiment/native-source-probe/normalized-first/program.bin',
}
assert sha(corpora['polynomial62']) == 'f482d60e6aee6f3187be2eac039b26fec8fb8cc7cef185d577b5d2f8ddad14a2'
for mode, corpus in corpora.items():
    output = base / f'{mode}.json'
    resources = base / f'{mode}-resources.json'
    binding_path = base / f'{mode}-binding.json'
    assert not any(p.exists() for p in (output, resources, binding_path))
    command = ['timeout', '300', str(exe), str(corpus), str(fixture), str(output), 'original']
    binding = {
        'scope': 'Fresh source-only native point discovery and encoded/decoded replay. No historical rules, numerical references, full closure or reciprocal necessity claim.',
        'build_binding_sha256': sha(build_path), 'executable_sha256': sha(exe),
        'source_corpus': str(corpus), 'source_corpus_sha256': sha(corpus),
        'fixture_sha256': sha(fixture), 'probe_source_sha256': sha(source),
        'launcher_sha256': sha(Path(__file__)), 'command': command,
    }
    binding_path.write_text(json.dumps(binding, indent=2) + '\n')
    r = subprocess.run([sys.executable, str(root / 'reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py'), str(resources), *command], cwd=root)
    assert sha(exe) == binding['executable_sha256'] and sha(corpus) == binding['source_corpus_sha256']
    assert sha(fixture) == binding['fixture_sha256'] and sha(source) == binding['probe_source_sha256']
    binding['exit_code'] = r.returncode
    binding['post_run_inputs_unchanged'] = True
    if output.exists(): binding['output_sha256'] = sha(output)
    binding_path.write_text(json.dumps(binding, indent=2) + '\n')
    print(mode, r.returncode, flush=True)
