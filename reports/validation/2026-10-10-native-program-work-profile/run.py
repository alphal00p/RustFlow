#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path.cwd();base=pathlib.Path(__file__).resolve().parent;mode=sys.argv[1];assert mode in ['single3','double'];h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
exe=pathlib.Path('/tmp/rustflow-native-program-work-profile');build=json.loads((base/'build-binding.json').read_text());assert h(exe)==build['executable']['sha256'];out=base/mode/'profile';assert not out.exists()
inputs={str(p):h(p)for p in [exe,base/mode/'input.bin',base/mode/'input.json',base/mode/'points.json',base/'probe.rs']};cmd=['timeout','600',str(exe),str(base/mode/'input.bin'),str(base/mode/'points.json'),str(out)]
binding={'scope':'Read-only native phase timing; no discovery or changed program. Shared host; run profile sectors sequentially.','command':cmd,'inputs':inputs,'build_binding_sha256':h(base/'build-binding.json'),'clock_ticks_per_second':int(subprocess.check_output(['getconf','CLK_TCK'],text=True))}
p=base/mode/'run-binding.json';assert not p.exists();p.write_text(json.dumps(binding,indent=2)+'\n')
r=subprocess.run([sys.executable,str(root/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(base/mode/'resources.json'),*cmd]);binding.update(exit_code=r.returncode,inputs_unchanged=all(h(pathlib.Path(p))==v for p,v in inputs.items()));p.write_text(json.dumps(binding,indent=2)+'\n');assert binding['inputs_unchanged'];raise SystemExit(r.returncode)
