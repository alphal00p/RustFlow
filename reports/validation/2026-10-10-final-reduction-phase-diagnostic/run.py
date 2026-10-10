from pathlib import Path
import os,json,hashlib,sys,subprocess
R=Path.cwd();B=Path(__file__).resolve().parent;S=R/'reports/validation/2026-10-10-final-history-work-profile';h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
j=json.loads((B/'build-binding.json').read_text());assert j['exit_code']==0;exe=Path(j['executable']['path']);assert h(exe)==j['executable']['sha256']
assert not (B/'run-binding.json').exists()
inputs=[B/'probe.rs',B/'build.py',B/'build-binding.json',B/'run.py',B/'input-bindings.json',exe]+[S/n for n in ['input.bin','input.json','points.json','profile/result.json']]
hashes={str(p):h(p) for p in inputs}
cmd=['timeout','600',str(exe),str(S/'input.bin'),str(S/'points.json'),str(B/'profile'),str(S/'profile/result.json')]
record={'scope':'Native full cold replay, seven terminal applicability checks, instrumented native-apply recursive scheduling on point24 control and point25 raised-target component only; no source discovery or numerical period values','command':cmd,'inputs_sha256':hashes,'budget':'600 wall seconds overall;480 seconds checked between native operations on one reduction; original100000 native applications/pending bounds retained','instrumentation':'Exact native apply API, current native harder-first pending comparator and coefficient operations; timing instrumentation adds I/O overhead; output checked against saved uninstrumented native result'}
(B/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n')
r=subprocess.run([sys.executable,str(R/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(B/'resources.json'),*cmd]);record.update(exit_code=r.returncode,inputs_unchanged=all(h(Path(p))==v for p,v in hashes.items()));(B/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n');assert record['inputs_unchanged'];raise SystemExit(r.returncode)
