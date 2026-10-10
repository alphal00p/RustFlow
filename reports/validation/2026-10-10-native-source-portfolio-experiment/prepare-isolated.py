#!/usr/bin/env python3
"""Reproduce the saved isolated native source copy; never edit production or run Cargo."""
import argparse,hashlib,json,shutil
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('destination',type=Path)
args=parser.parse_args()
report=Path(__file__).resolve().parent;root=report.parents[2]
record=json.loads((report/'isolated-source-map.json').read_text());base=root/'vendor/rustred'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert not args.destination.exists(),'refuse to overwrite an existing copy'
for name,info in record['files'].items():
 original=base/name;saved=report/info['saved_source']
 assert sha(saved)==info['after_sha256']
 assert (sha(original) if original.exists() else None)==info['before_sha256']
shutil.copytree(base,args.destination)
for name,info in record['files'].items():
 destination=args.destination/name;destination.parent.mkdir(parents=True,exist_ok=True)
 shutil.copyfile(report/info['saved_source'],destination)
 assert sha(destination)==info['after_sha256']
print(json.dumps({'status':'isolated_copy_reproduced','destination':str(args.destination),'modified_files':len(record['files']),'production_modified':False,'cargo_invoked':False}))
