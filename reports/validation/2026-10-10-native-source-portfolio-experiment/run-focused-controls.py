#!/usr/bin/env python3
"""Run the frozen cfg-only assertion controls against the exact isolated library."""
import hashlib,json,subprocess,sys
from pathlib import Path
base=Path(__file__).resolve().parent;root=base.parents[2]
build_path=base/'focused-runner-build.json';build=json.loads(build_path.read_text());assert build['exit_code']==0
native_path=base/'native-build/build-binding.json';native=json.loads(native_path.read_text());assert native['exit_code']==0
output=base/'focused-controls.json';resources=base/'focused-controls-resources.json';binding_path=base/'focused-controls-binding.json'
assert not output.exists() and not resources.exists() and not binding_path.exists()
exe=Path(build['executable']['path']);lib=Path(native['output'])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(exe)==build['executable']['sha256'] and sha(lib)==native['output_sha256']
source_map=json.loads((base/'isolated-source-map.json').read_text())
for name,meta in source_map['files'].items():assert sha(Path(source_map['workspace_copy'])/name)==meta['after_sha256']
command=['timeout','120',str(exe),str(output)]
binding={'scope':'19 cfg-only native assertion controls; no physical periods or whole closure.','native_build_binding_sha256':sha(native_path),'runner_build_binding_sha256':sha(build_path),'library_sha256':sha(lib),'executable_sha256':sha(exe),'source_map_sha256':sha(base/'isolated-source-map.json'),'command':command,'launcher_sha256':sha(Path(__file__))}
binding_path.write_text(json.dumps(binding,indent=2)+'\n')
r=subprocess.run([sys.executable,str(root/'reports/validation/2026-10-10-native-zero-domain-projection/run-resource-command.py'),str(resources),*command],cwd=root)
assert sha(exe)==binding['executable_sha256'] and sha(lib)==binding['library_sha256']
for name,meta in source_map['files'].items():assert sha(Path(source_map['workspace_copy'])/name)==meta['after_sha256']
binding['exit_code']=r.returncode;binding['post_run_sources_and_binaries_unchanged']=True
if output.exists():binding['output_sha256']=sha(output)
binding_path.write_text(json.dumps(binding,indent=2)+'\n');raise SystemExit(r.returncode)
