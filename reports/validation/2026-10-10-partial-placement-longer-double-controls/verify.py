#!/usr/bin/env python3
"""Verify report and original capsule hashes, resolving the moved executable."""
import hashlib
import json
from pathlib import Path
BASE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((BASE/'artifact-manifest.json').read_text())
for path,record in manifest['files'].items():
    p=BASE/path
    assert p.stat().st_size==record['bytes'] and sha(p)==record['sha256'],path
move=json.loads((BASE/'executable-relocation.json').read_text())
for record in json.loads((BASE/'capsule-manifest.json').read_text())['files']:
    p=Path(move['current_path']) if record['capsule']==move['relative_original_path'] else BASE/record['capsule']
    assert p.stat().st_size==record['bytes'] and sha(p)==record['sha256'],record['capsule']
print(json.dumps({'report_files':len(manifest['files']),'capsule_including_relocated_static_binary':'PASS'}))
