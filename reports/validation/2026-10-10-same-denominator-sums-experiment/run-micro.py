from pathlib import Path
import json,hashlib,subprocess,sys
R=Path.cwd();B=Path(__file__).resolve().parent;h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert json.loads((B/'actual-run-binding.json').read_text())['exit_code']==0
build=json.loads((B/'micro-attempt2-build-binding.json').read_text());assert build['exit_code']==0
exe=Path(build['executable']['path']);assert h(exe)==build['executable']['sha256'];assert not(B/'micro-run-binding.json').exists()
inputs=[B/'micro-corrected.rs',B/'build-micro-corrected.py',B/'micro-attempt2-build-binding.json',B/'run-micro.py',B/'actual-run-binding.json',B/'actual/samples/samples.json',B/'actual/samples/symbolica-state.bin',B/'actual/samples/coefficients.bin',exe]
cas=Path('/home/codex-2/.cargo/git/checkouts/symbolica-db3dbb7e8d40efb5/ed2374f/src/domains/rational_polynomial.rs');inputs.append(cas)
hashes={str(p):h(p)for p in inputs};cmd=['timeout','300',str(exe),str(B/'actual/samples'),str(B/'micro-results.json')]
record={'scope':'three exact addition methods, three alternating repetitions on each frozen real repeated-denominator buffer; no native rule discovery or target reduction','command':cmd,'inputs_sha256':hashes,'maximum_samples':8,'maximum_coefficients_per_sample':32,'wall_limit_seconds':300,'symbolica_owner_anchors':{'integer_from_num_den':[423,459],'rational_add':[1075,1135]}};(B/'micro-run-binding.json').write_text(json.dumps(record,indent=2)+'\n')
r=subprocess.run([sys.executable,str(R/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(B/'micro-resources.json'),*cmd]);record.update(exit_code=r.returncode,inputs_unchanged=all(h(Path(p))==v for p,v in hashes.items()));(B/'micro-run-binding.json').write_text(json.dumps(record,indent=2)+'\n');assert record['inputs_unchanged'];raise SystemExit(r.returncode)
