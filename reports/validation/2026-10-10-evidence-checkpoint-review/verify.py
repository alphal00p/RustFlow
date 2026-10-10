"""Verify complete local manifests before checkpoint staging; no numerical claim."""
from pathlib import Path
import hashlib,json
ROOT=Path.cwd();HERE=Path(__file__).resolve().parent
names=['final-history-work-profile','final-reduction-phase-diagnostic','historical-candidate-collection-assessment','historical-candidate-limit-draft','longer-singleton-control','partial-boundary-certificate-final-review','partial-boundary-certificate-review','partial-boundary-certificate-validation','partial-placement-admission-design','partial-placement-boundary-assessment','partial-placement-closure-pilot','partial-placement-longer-double-controls','partial-placement-support-proof','partial-singleton-endpoint-assessment']
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for name in names:
 base=ROOT/('reports/validation/2026-10-10-'+name);manifest=base/'artifact-manifest.json'
 if not manifest.exists():manifest=base/'manifest.json'
 data=json.loads(manifest.read_text());files=data['files']
 entries=[{'path':key,**value}for key,value in files.items()]if isinstance(files,dict)else files
 for row in entries:
  p=(ROOT if row['path'].startswith('reports/')else base)/row['path'];assert p.resolve().is_relative_to(base.resolve())and p.is_file(),p
  assert h(p)==row['sha256'],p
  if 'bytes'in row:assert p.stat().st_size==row['bytes'],p
 rows.append({'report':str(base.relative_to(ROOT)),'manifest':manifest.name,'manifest_sha256':h(manifest),'files_verified':len(entries)})
report={'scope':'Hash/byte verification of every local file listed by the14 completed report manifests. This does not rerun computations, restore relocated executables, assert historical root artifacts still exist, or change the explicitly limited mathematical/numerical claims.','status':'passed','reports':rows,'total_files_verified':sum(r['files_verified']for r in rows),'script_sha256':h(Path(__file__))}
(HERE/'review.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'passed','reports':len(rows),'files_verified':report['total_files_verified']}))
