#!/usr/bin/env python3
"""Capture the completed shared build and its extra exact source-equivalence test."""
import json
import re
import subprocess
import sys
from pathlib import Path

BASE=Path(__file__).resolve().parent
PREFIX='polynomial-closure'

def main():
    log=(BASE/(PREFIX+'-build-resources.log')).read_text()
    resources=json.loads((BASE/(PREFIX+'-build-resources.json')).read_text())
    assert resources['exit_code']==0 and re.search(r'Finished .*release',log)
    emitted=re.findall(r'^\s*Executable .+ \(([^()]+)\)$',log,re.MULTILINE)
    extra=[name for name in emitted if Path(name).name.startswith('finite_density_polynomial_sources-')]
    assert len(extra)==1,'expected exactly one generated-source equivalence executable'
    subprocess.run([sys.executable,str(BASE/'capture-frozen-build.py'),'--prefix',PREFIX,
                    '--extra-artifact',extra[0]],check=True)

if __name__=='__main__':main()
