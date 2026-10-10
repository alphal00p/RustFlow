#!/usr/bin/env python3
"""Run fresh coherent native and RustFlow gates after capture; never invokes Cargo."""
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
PREFIX = 'verified-reuse'
GATES = ['native', 'lib', 'flow-boundary', 'massless-sources',
         'source-fingerprint', 'python', 'capacity', 'polynomial-sources']


def main():
    # Native ownership changed. Reuse of an old executable/gate is deliberately
    # not accepted; the coherent build's full guarded suite runs first.
    for gate in GATES:
        subprocess.run([sys.executable, str(BASE / 'run-frozen-root-gate.py'),
                        '--prefix', PREFIX, gate], check=True)
    subprocess.run([sys.executable, str(BASE / 'summarize-frozen-root-gates.py'),
                    '--prefix', PREFIX, '--gates', *GATES], check=True)


if __name__ == '__main__':
    main()
