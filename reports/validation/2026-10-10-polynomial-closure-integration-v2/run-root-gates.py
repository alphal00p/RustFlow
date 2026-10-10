#!/usr/bin/env python3
"""Run the assigned prebuilt gates after build capture; never invokes Cargo.

This first authenticates reuse of the unchanged prior native gate. If source,
feature or binary identity differs, it stops for a separately recorded native
rerun decision. Full numerical acceptance is a separate recorded run.
"""
import subprocess
import sys
from pathlib import Path

BASE=Path(__file__).resolve().parent
PREFIX='polynomial-closure-v2'
GATES=['lib','flow-boundary','massless-sources','source-fingerprint','python','capacity','polynomial-sources']

def main():
    subprocess.run([sys.executable,str(BASE/'reuse-native-gate.py'),'--prefix',PREFIX],check=True)
    for gate in GATES:
        subprocess.run([sys.executable,str(BASE/'run-frozen-root-gate.py'),'--prefix',PREFIX,gate],check=True)
    subprocess.run([sys.executable,str(BASE/'summarize-frozen-root-gates.py'),'--prefix',PREFIX,'--gates',*GATES],check=True)

if __name__=='__main__':main()
