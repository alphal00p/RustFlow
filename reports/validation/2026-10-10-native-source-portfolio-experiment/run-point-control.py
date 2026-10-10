#!/usr/bin/env python3
"""Run a bounded native source-only point diagnostic; no physical/reference values."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('tag',choices=['baseline','proposal']);p.add_argument('--mode',choices=['original'],default='original');p.add_argument('--witness',action='store_true');a=p.parse_args()
base=Path(__file__).resolve().parent;root=base.parents[2];name=a.tag+'-'+a.mode+('-odd-witness' if a.witness else '-19-points');output=base/(name+'.json');resources=base/(name+'-resources.json');binding_path=base/(name+'-binding.json')
assert not output.exists() and not resources.exists() and not binding_path.exists()
build_path=base/(a.tag+'-probe-build.json');build=json.loads(build_path.read_text());assert build['exit_code']==0
exe=Path(build['executable']['path']);fixture=base/('odd-zero-point.json' if a.witness else 'points.json');corpus=root/'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-singleton-source-diagnostic/round-006-provisional.bin'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
assert sha(exe)==build['executable']['sha256']
command=['timeout','600',str(exe),str(corpus),str(fixture),str(output),a.mode]
binding={'scope':'Bounded original-source discovery and proof replay; Point-only full-source portfolio control; no physical values or closure claim.','build_binding_sha256':sha(build_path),'executable_sha256':sha(exe),'source_corpus_sha256':sha(corpus),'fixture_sha256':sha(fixture),'launcher_sha256':sha(Path(__file__)),'command':command}
binding_path.write_text(json.dumps(binding,indent=2)+'\n')
r=subprocess.run([sys.executable,str(root/'reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py'),str(resources),*command],cwd=root)
assert sha(exe)==binding['executable_sha256'] and sha(corpus)==binding['source_corpus_sha256'] and sha(fixture)==binding['fixture_sha256']
binding['exit_code']=r.returncode;binding['post_run_inputs_unchanged']=True
if output.exists():binding['output_sha256']=sha(output)
binding_path.write_text(json.dumps(binding,indent=2)+'\n');raise SystemExit(r.returncode)
