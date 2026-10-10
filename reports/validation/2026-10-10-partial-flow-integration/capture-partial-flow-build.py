"""Capture the completed thirteen-artifact shared build; no Cargo invocation."""
import subprocess
import sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
if __name__=='__main__':
    subprocess.run([sys.executable,str(BASE/'capture-frozen-build.py'),'--prefix','partial-flow'],check=True)
