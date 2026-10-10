#!/usr/bin/env python3
"""Run fresh coherent native and RustFlow gates after capture; never invokes Cargo."""
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
PREFIX = 'terminal-promotion'
GATES = ['native', 'lib', 'flow-boundary', 'massless-sources',
         'source-fingerprint', 'python', 'capacity', 'polynomial-sources']


def main():
    # Fresh native guarded suite includes ten memo regressions; the root
    # finite_density filter includes three wrapper/final-audit regressions.
    for gate in GATES:
        subprocess.run([sys.executable, str(BASE / 'run-frozen-root-gate.py'),
                        '--prefix', PREFIX, gate], check=True)
    subprocess.run([sys.executable, str(BASE / 'summarize-frozen-root-gates.py'),
                    '--prefix', PREFIX, '--gates', *GATES], check=True)


if __name__ == '__main__':
    main()
