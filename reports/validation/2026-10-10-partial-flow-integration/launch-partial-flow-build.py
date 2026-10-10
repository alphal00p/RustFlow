"""Capture and build the frozen private partial-origin/source/endpoint/boundary flow integration.
Importing does not snapshot or build; explicit source-owner freeze is required.
"""
import argparse
import subprocess
import sys
from pathlib import Path
from frozen_sources import digest

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PREFIX='partial-flow'
SOURCE_FIXTURE='reports/validation/2026-10-10-polynomial-raw-ward-attribution/native-source-probe/strongest62-boost-ward/program.bin'
SOURCE_FIXTURE_SHA256='f482d60e6aee6f3187be2eac039b26fec8fb8cc7cef185d577b5d2f8ddad14a2'
REQUESTED_FIXTURE='fixtures/finite_density/requested_point_source_rows.json'
REQUESTED_FIXTURE_SHA256='789aacca33f653e1e8df7a02e5434f1ea3690f18c057253fa8789aa2a6e4753c'
CARGO=['cargo','test','--locked','--offline','--release','-p','symbolica-amflow','-p','rustred',
    '--features','symbolica-amflow/python','--lib','--test','finite_density_runtime_flow',
    '--test','finite_density_flow_boundary','--test','finite_density_massless_sources',
    '--test','finite_density_massless_flow',
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
        'tests/finite_density_runtime_flow.rs','tests/finite_density_polynomial_sources.rs',
        'src/finite_density/singleton_zero.rs',
        'vendor/rustred/crates/rustred-core/src/solver/guarded/lifecycle/unit_memo.rs',
        'vendor/rustred/crates/rustred-core/src/solver/guarded/lifecycle/unit_memo/tests.rs',
        'src/finite_density/partial_origin.rs','src/finite_density/partial_resources.rs',
        'src/finite_density/partial_endpoint.rs','src/finite_density/partial_continuation.rs',
        'src/finite_density/source_class.rs','src/finite_density/source_class/tests.rs',
        'src/finite_density/boundary/virtual_soft.rs','src/finite_density/boundary/partial_tests.rs',
        'src/finite_density/flow_partial_tests.rs']
    for name in required:
        if not(ROOT/name).is_file():raise ValueError('required frozen source missing: '+name)
    assert 'unit_reduction_memo' in (ROOT/'src/finite_density/reduction.rs').read_text(), 'reviewed memo option must be applied before freeze'
    assert 'RUSTFLOW_WEIGHTED_UNIT_MEMO' in (ROOT/'tests/finite_density_runtime_flow.rs').read_text(), 'memo runtime control must be included'
    assert 'pub struct PreparedSingletonZero' in (ROOT/'src/finite_density/singleton_zero.rs').read_text(), 'reviewed singleton physical-zero owner must be applied before freeze'
    assert 'audit_and_bind_reduced_system' in (ROOT/'src/finite_density/flow.rs').read_text()
    assert 'RUSTFLOW_WEIGHTED_SHIFTED_SLOTS' in (ROOT/'tests/finite_density_runtime_flow.rs').read_text()
    if not(ROOT/'tests/finite_density_polynomial_sources.rs').is_file():
        raise ValueError('the exact generated-source comparison gate must be integrated before freeze')
    if digest(ROOT/REQUESTED_FIXTURE)!=REQUESTED_FIXTURE_SHA256:
        raise ValueError('frozen requested-point source fixture changed')
    if digest(ROOT/SOURCE_FIXTURE)!=SOURCE_FIXTURE_SHA256:
        raise ValueError('frozen source-only comparison corpus changed')
    subprocess.run([sys.executable,str(BASE/'capture-source-snapshot.py'),'--prefix',PREFIX,
        '--confirm-frozen','--include-source',SOURCE_FIXTURE,'--include-source',REQUESTED_FIXTURE,
        '--scope','Frozen partial occupied-flow integration with private origin, source-image, finite-label endpoint and virtual-soft boundary owners. Complete native replay and actual combined-target transport remain mandatory. Existing all-shift regressions retain source/discovery settings and numerical thresholds; no partial numerical result is presumed.'],cwd=ROOT,check=True)
    command=['nix','develop','--command',sys.executable,str(BASE/'run-resource-command.py'),
             str(BASE/(PREFIX+'-build-resources.json')),*CARGO]
    raise SystemExit(subprocess.run(command,cwd=ROOT).returncode)

if __name__=='__main__':main()
