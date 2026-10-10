#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path.cwd();base=pathlib.Path(__file__).resolve().parent;mode=sys.argv[1];assert mode in ['original','extended-energy-Ward']
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
exe=pathlib.Path('/tmp/rustflow-singleton-energy-ward-probe');build=json.loads((base/'build-binding.json').read_text());assert h(exe)==build['executable']['sha256'];out=base/mode;assert not out.exists()
inputs={str(p):h(p)for p in [exe,base/'input.bin',base/'points.json',base/'density-input.json',base/'probe.rs']}
cmd=['timeout','300',str(exe),str(base/'input.bin'),str(base/'points.json'),str(out),mode,'3',str(base/'density-input.json')]
binding={'scope':'Fresh source-only polynomial energy Ward comparison; historical replay only validates actual NoRule selection. No numerical periods or imported rules.','inputs':inputs,'build_binding_sha256':h(base/'build-binding.json'),'command':cmd}
p=base/f'{mode}-binding.json';p.write_text(json.dumps(binding,indent=2)+'\n')
r=subprocess.run([sys.executable,str(root/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(base/f'{mode}-resources.json'),*cmd]);binding.update(exit_code=r.returncode,inputs_unchanged=all(h(pathlib.Path(p))==v for p,v in inputs.items()));p.write_text(json.dumps(binding,indent=2)+'\n');assert binding['inputs_unchanged'];raise SystemExit(r.returncode)
