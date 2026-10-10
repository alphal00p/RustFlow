from pathlib import Path
import hashlib,json
ROOT=Path('/common/dev/rustflow_fermi'); REPORT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((REPORT/'artifact-manifest.json').read_text())
for name,row in m['files'].items():
 p=REPORT/name;assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],name
binding=json.loads((REPORT/'source-and-proof-binding.json').read_text())
for name,row in binding['proof_inputs'].items():assert sha(ROOT/name)==row['sha256'],name
build=json.loads((REPORT/'build-01.json').read_text())
for name,digest in build['source_files'].items():assert sha(REPORT/'source'/name)==digest,name
review=json.loads((REPORT/'independent-review.json').read_text())
for name,digest in review['files_sha256'].items():assert sha(REPORT/name)==digest,name
print(json.dumps({'status':'passed','manifest_files':len(m['files']),'proof_bindings':len(binding['proof_inputs']),'source_build_and_independent_review_bindings':True,'no_build_or_numerics':True}))
