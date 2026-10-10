#!/usr/bin/env python3
"""Bounded native prediction launch; never reads the reference artifact."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]


def digest(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def main():
    binding_path = BASE / 'native-build-binding.json'
    binding = json.loads(binding_path.read_text())
    assert binding['exit_code'] == 0
    for name, sha in binding['inputs_sha256'].items():
        assert digest(Path(name)) == sha
    executable = Path(binding['executable'])
    assert digest(executable) == binding['executable_sha256']
    predictions = BASE / 'native-predictions'
    assert not predictions.exists(), 'preserve previous runs'
    runner = ROOT / 'reports/validation/2026-10-10-polynomial-closure-integration/run-resource-command.py'
    command = [sys.executable, str(runner), str(BASE / 'native-resources.json'),
               'timeout', '300', str(executable), str(predictions)]
    before = {str(p): digest(p) for p in [binding_path, executable, Path(__file__), runner]}
    result = subprocess.run(command, cwd=ROOT)
    assert all(digest(Path(name)) == sha for name, sha in before.items())
    assert all(digest(Path(name)) == sha for name, sha in binding['inputs_sha256'].items())
    (BASE / 'native-run-binding.json').write_text(json.dumps({
        'scope': 'Native predictions saved before comparison; no reference file opened.',
        'inputs_sha256': before, 'command': command, 'exit_code': result.returncode,
        'post_run_bindings_unchanged': True, 'reference_files_read': 0,
    }, indent=2) + '\n')
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
