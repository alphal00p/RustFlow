from pathlib import Path
import json,hashlib,gzip,shutil,sys,subprocess
R=Path.cwd();B=Path(__file__).resolve().parent;S=R/'reports/validation/2026-10-10-same-denominator-sums-experiment';T=Path('/tmp/rustflow-raw-lcm-samples-20261010');h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert not (B/'run-binding.json').exists();T.mkdir(exist_ok=False)
archive=json.loads((S/'archive-map.json').read_text())['files'][0];assert archive['original_path']=='actual/samples/coefficients.bin';compressed=S/archive['archive_path'];assert h(compressed)==archive['archive_sha256'];data=gzip.decompress(compressed.read_bytes());assert hashlib.sha256(data).hexdigest()==archive['original_sha256'];(T/'coefficients.bin').write_bytes(data)
for n in ['samples.json','symbolica-state.bin']:shutil.copyfile(S/'actual/samples'/n,T/n)
build=json.loads((B/'actual-attempt1-build-binding.json').read_text());assert build['exit_code']==0;exe=Path(build['executable']['path']);assert h(exe)==build['executable']['sha256']
cas=Path('/home/codex-2/.cargo/git/checkouts/symbolica-db3dbb7e8d40efb5/ed2374f/src/domains/rational_polynomial.rs')
assert json.loads((B/'tests-attempt2-resources.json').read_text())['exit_code']==0
inputs=[B/'tests-attempt2-resources.json',B/'probe.rs',B/'build.py',B/'run.py',B/'actual-attempt1-build-binding.json',S/'artifact-manifest.json',S/'archive-map.json',compressed,cas,exe]+list(T.iterdir());hashes={str(p):h(p)for p in inputs};cmd=['timeout','300',str(exe),str(T),str(B/'results.json')]
record={'scope':'only eight restored frozen coefficient samples; canonical raw-LCM with preflight-only fallback versus native sequential/balanced sums','command':cmd,'inputs_sha256':hashes,'restored_archive':archive,'wall_limit_seconds':300,'maximum_samples':8,'maximum_coefficients':32,'symbolica_revision':'ed2374f','owner_lines':{'integer_from_num_den':[423,459],'rational_add':[1075,1135]}};(B/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n')
r=subprocess.run([sys.executable,str(R/'reports/validation/2026-10-10-requested-ray-point-integration/run-resource-command.py'),str(B/'resources.json'),*cmd]);record.update(exit_code=r.returncode,inputs_unchanged=all(h(Path(p))==v for p,v in hashes.items()));(B/'run-binding.json').write_text(json.dumps(record,indent=2)+'\n');assert record['inputs_unchanged'];raise SystemExit(r.returncode)
