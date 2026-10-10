import hashlib,json,subprocess,sys,os
from pathlib import Path
ROOT=Path.cwd();HERE=Path(__file__).resolve().parent
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
build=json.loads((HERE/'build-binding.json').read_text());assert build['exit_code']==0
exe=Path(build['executable']['path']);assert h(exe)==build['executable']['sha256']
assert not (HERE/'run-binding.json').exists() and not (HERE/'profile').exists()
inputs=[HERE/n for n in ['probe.rs','build.py','build-binding.json','input-bindings.json','input.bin','input.json','points.json','run.py']]+[exe]
hashes={str(p):h(p)for p in inputs}
command=['timeout','600',str(exe),str(HERE/'input.bin'),str(HERE/'points.json'),str(HERE/'profile')]
record={'scope':'Bounded standalone profile of one copied round18 provisional program; full cold replay followed by42 selected actual frontier/history label reductions. All libraries are statically linked and runtime evidence is independent of later root rebuilds. No closure or whole-workflow claim.', 'command':command,'inputs_sha256':hashes,'clock_ticks_per_second':os.sysconf('SC_CLK_TCK')}
(HERE/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n')
r=subprocess.run([sys.executable,str(ROOT/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(HERE/'resources.json'),*command])
record.update(exit_code=r.returncode,post_run_inputs_unchanged=all(h(Path(p))==sha for p,sha in hashes.items()))
(HERE/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n');assert record['post_run_inputs_unchanged'];raise SystemExit(r.returncode)
