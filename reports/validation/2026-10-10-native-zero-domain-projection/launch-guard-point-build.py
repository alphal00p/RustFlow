"""Launch the next shared release only after an explicit owner freeze notice.

Preparing/importing this module captures no snapshot and starts no build. Old
zero-projection artifacts and shared helper scripts remain byte-for-byte intact.
"""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PREFIX='guard-point'
CARGO=['cargo','test','--locked','--offline','--release','-p','symbolica-amflow','-p','rustred',
    '--features','symbolica-amflow/python','--lib','--test','finite_density_runtime_flow',
    '--test','finite_density_flow_boundary','--test','finite_density_massless_sources',
    '--test','source_fingerprint','--test','python_finite_density','--test','finite_density_capacity','--no-run']

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-ready-confirmed',action='store_true',help='Use only after the source owner explicitly confirms the batch is ready and frozen')
    args=parser.parse_args()
    if not args.source_ready_confirmed:parser.error('source owner readiness/freeze is required')
    if any(BASE.glob(PREFIX+'-source-hashes.json')) or any(BASE.glob(PREFIX+'-build-resources.*')):
        raise ValueError('preserve prior attempts; a changed source batch requires another prefix')
    subprocess.run([sys.executable,str(BASE/'capture-source-snapshot.py'),'--prefix',PREFIX,'--confirm-frozen',
        '--scope','Frozen bounded native point refinement for conditional-rule failures; includes existing conservative zero-domain projection and guarded schema v2. Runtime numerical guard_passes=1 is a separately recorded configuration change.'],cwd=ROOT,check=True)
    command=['nix','develop','--command',sys.executable,str(BASE/'run-resource-command.py'),str(BASE/(PREFIX+'-build-resources.json')),*CARGO]
    result=subprocess.run(command,cwd=ROOT)
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
