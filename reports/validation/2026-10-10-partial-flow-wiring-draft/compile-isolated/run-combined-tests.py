"""Run the new proof and wiring tests in the one complete isolated crate."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PREFIXES = (
    'finite_density::partial_origin::',
    'finite_density::source_class::',
    'finite_density::flow::partial_tests::',
    'finite_density::boundary::partial_tests::',
    'finite_density::boundary::virtual_soft::',
)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, data):
    assert not path.exists(), path
    path.write_text(json.dumps(data, indent=2) + '\n')


def main(revision, build_attempt, gate_attempt):
    for name in [revision, build_attempt, gate_attempt]:
        assert name and all(c.isalnum() or c in '-_' for c in name)
    state = HERE / revision
    prepared_path = state / 'prepared.json'
    prepared = json.loads(prepared_path.read_text())
    result_path = state / build_attempt / 'result.json'
    result = json.loads(result_path.read_text())
    assert result['exit_code'] == 0 and result['post_integrity_passed']
    assert result['mode'] == 'test' and result['prepared_sha256'] == digest(prepared_path)
    executable = Path(result['output']['path'])
    assert digest(executable) == result['output']['sha256']
    spec = importlib.util.spec_from_file_location('full_crate_builder', HERE / 'build-full-crate.py')
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    assert digest(HERE / 'build-full-crate.py') == prepared['builder_sha256']
    builder.verify(prepared)
    folder = state / gate_attempt
    folder.mkdir()
    listing = subprocess.run([str(executable), '--list'], text=True, capture_output=True, check=True)
    (folder / 'test-list.log').write_text(listing.stdout + listing.stderr)
    names = sorted(line.removesuffix(': test') for line in listing.stdout.splitlines()
                   if line.endswith(': test') and line.startswith(PREFIXES))
    assert len(names) == len(set(names)) and 40 <= len(names) <= 100
    assert all(any(name.startswith(prefix) for name in names) for prefix in PREFIXES)
    assert any(name.endswith('final_audit_includes_raw_derivatives_candidates_and_combined_pole_weights')
               for name in names)
    assert any(name.endswith('actual_partial_regions_use_bound_virtual_zeros_and_ordinary_hard_seeds')
               for name in names)
    command = [shutil.which('timeout'), '--signal=TERM', '--kill-after=30', '900',
               str(executable), '--exact', *names, '--test-threads=1', '--nocapture']
    assert command[0] is not None
    wrapper = builder.BASE / 'run-resource-command.py'
    resources = folder / 'resources.json'
    # These proof tests derive their own fixtures. No numerical predictions or
    # oracle values are supplied by the launcher environment.
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith('RUSTFLOW_')}
    environment['CARGO_MANIFEST_DIR'] = prepared['source_root']
    save(folder / 'launch.json', {
        'scope': 'Selected proof, source-image, placement and boundary tests in one full copied crate; no production/numerical admission.',
        'prepared_sha256': digest(prepared_path), 'build_result_sha256': digest(result_path),
        'executable': result['output'], 'selected_prefixes': PREFIXES, 'selected_tests': names,
        'test_list_sha256': digest(folder / 'test-list.log'), 'command': command,
        'resource_wrapper': str(wrapper), 'resource_wrapper_sha256': digest(wrapper),
        'launcher_sha256': digest(Path(__file__)), 'environment_policy': 'Remove inherited RUSTFLOW_*; use copied crate manifest directory.'})
    completed = subprocess.run([sys.executable, str(wrapper), str(resources), *command],
                               cwd=prepared['source_root'], env=environment)
    integrity_error = None
    try:
        builder.verify(prepared)
        assert digest(executable) == result['output']['sha256']
    except Exception as error:
        integrity_error = repr(error)
    log = resources.with_suffix('.log').read_text()
    rows = re.findall(r'^test (\S+) \.\.\. (ok|FAILED|ignored)$', log, re.MULTILINE)
    observed = {name: status for name, status in rows}
    summary = re.findall(r'test result: (ok|FAILED)\. (\d+) passed; (\d+) failed; (\d+) ignored;', log)
    parsed = summary[-1] if summary else None
    valid_count = (len(rows) == len(names) and set(observed) == set(names)
                   and parsed is not None and int(parsed[1]) == len(names)
                   and int(parsed[2]) == 0 and int(parsed[3]) == 0)
    passed = completed.returncode == 0 and integrity_error is None and valid_count
    save(folder / 'summary.json', {
        'status': 'passed' if passed else 'failed', 'selected_count': len(names),
        'observed_cases': observed, 'parsed_harness_summary': parsed,
        'all_selected_cases_executed_once_and_passed': valid_count,
        'exit_code': completed.returncode, 'post_integrity_passed': integrity_error is None,
        'post_integrity_error': integrity_error, 'resources': json.loads(resources.read_text()),
        'artifacts': {p.name: digest(p) for p in sorted(folder.iterdir()) if p.is_file()},
        'scope_limit': 'No partial-flow numerical amplitude or three-/four-loop acceptance is established.'})
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--build-attempt', required=True)
    parser.add_argument('--gate-attempt', required=True)
    args = parser.parse_args()
    main(args.revision, args.build_attempt, args.gate_attempt)
