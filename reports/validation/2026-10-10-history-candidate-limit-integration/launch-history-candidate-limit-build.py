"""Prepare the coherent shared release only after the source owner announces freeze.

Importing this module captures no snapshot and starts no build. This batch adds
an optional active-only historical candidate-map limit after unchanged mandatory
source replay and all original-target/basis-derivative audits. Native proof and
application semantics, source policy and scheduler remain unchanged.
"""
import argparse
import subprocess
import sys
from pathlib import Path
from frozen_sources import digest

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PREFIX='history-candidate-limit'
SOURCE_FIXTURE='reports/validation/2026-10-10-polynomial-raw-ward-attribution/native-source-probe/strongest62-boost-ward/program.bin'
SOURCE_FIXTURE_SHA256='f482d60e6aee6f3187be2eac039b26fec8fb8cc7cef185d577b5d2f8ddad14a2'
REQUESTED_FIXTURE='fixtures/finite_density/requested_point_source_rows.json'
REQUESTED_FIXTURE_SHA256='789aacca33f653e1e8df7a02e5434f1ea3690f18c057253fa8789aa2a6e4753c'
CARGO=['cargo','test','--locked','--offline','--release','-p','symbolica-amflow','-p','rustred',
    '--features','symbolica-amflow/python','--lib','--test','finite_density_runtime_flow',
    '--test','finite_density_flow_boundary','--test','finite_density_massless_sources',
    '--test','source_fingerprint','--test','python_finite_density','--test','finite_density_capacity',
    '--test','finite_density_polynomial_sources','--no-run']

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-ready-confirmed',action='store_true')
    args=parser.parse_args()
    if not args.source_ready_confirmed:parser.error('explicit source owner readiness and freeze is required')
    if any(BASE.glob(PREFIX+'-source-hashes.json'))or any(BASE.glob(PREFIX+'-build-resources.*')):
        raise ValueError('preserve existing attempts; use a distinct prefix for a changed source batch')
    required=['vendor/rustred/crates/rustred-core/src/solver/guarded/verified_tests.rs','src/finite_density/reduction/requested.rs','src/finite_density/reduction/requested_tests.rs',
        'src/finite_density/reduction/active.rs','src/finite_density/guarded.rs',
        'tests/finite_density_runtime_flow.rs','tests/finite_density_polynomial_sources.rs']
    for name in required:
        if not(ROOT/name).is_file():raise ValueError('required frozen source missing: '+name)
    assert 'pub max_history_candidate_maps: Option<usize>' in (ROOT/'src/finite_density/reduction.rs').read_text(), 'reviewed candidate-limit option must be applied before freeze'
    if not(ROOT/'tests/finite_density_polynomial_sources.rs').is_file():
        raise ValueError('the exact generated-source comparison gate must be integrated before freeze')
    if digest(ROOT/REQUESTED_FIXTURE)!=REQUESTED_FIXTURE_SHA256:
        raise ValueError('frozen requested-point source fixture changed')
    if digest(ROOT/SOURCE_FIXTURE)!=SOURCE_FIXTURE_SHA256:
        raise ValueError('frozen source-only comparison corpus changed')
    subprocess.run([sys.executable,str(BASE/'capture-source-snapshot.py'),'--prefix',PREFIX,
        '--confirm-frozen','--include-source',SOURCE_FIXTURE,'--include-source',REQUESTED_FIXTURE,
        '--scope','Frozen optional active-only historical candidate-map limit, coverage diagnostics and phase events. Default None retains full collection; mandatory source replay, targets and basis derivatives are unchanged. No native algorithm or source-policy change.'],cwd=ROOT,check=True)
    command=['nix','develop','--command',sys.executable,str(BASE/'run-resource-command.py'),
             str(BASE/(PREFIX+'-build-resources.json')),*CARGO]
    raise SystemExit(subprocess.run(command,cwd=ROOT).returncode)

if __name__=='__main__':main()
