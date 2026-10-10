#!/usr/bin/env python3
"""Freeze completed guard-point validation evidence after reader-released archival."""
import datetime,gzip,hashlib,json
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
out=BASE/'guard-point-artifacts.json'
if out.exists():raise ValueError('preserve existing frozen manifest')
snapshot_path=BASE/'guard-point-source-hashes.json';snapshot=verify_snapshot(snapshot_path)
status=json.loads((BASE/'guard-point-validation-status.json').read_text())
assert status['checkpoint_regression_scope_validated'] and not status['full_feature_validated']
archive_path=BASE/'guard-point-native-archives.json'
assert archive_path.exists(),'native archival must finish before capture'
archive=json.loads(archive_path.read_text())
archive_check=json.loads((BASE/'guard-point-native-archive-verification.json').read_text())
assert archive_check['status'] in ['pass','passed']
records=archive.get('records',archive.get('files'))
assert isinstance(records,list)
for row in records:
 ap=ROOT/row['archive_path'];assert digest(ap)==row['archive_sha256']
 h=hashlib.sha256();n=0
 with gzip.open(ap,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b);n+=len(b)
 assert h.hexdigest()==row['original_sha256'] and n==row['original_bytes']
 assert not (ROOT/row['original_path']).exists()
# Only fresh stage evidence; old failed stage is linked without rewriting its bytes.
files={p for p in BASE.rglob('*') if p.is_file() and p.relative_to(BASE).parts[0].startswith('guard-point') and p!=out}
helpers=['capture-guard-point-artifacts.py','summarize-guard-point-checkpoint.py','launch-guard-point-build.py','run-guard-point-gate.py','run-guard-point-massive-regression.py','compare-guard-point-default-path.py','frozen_sources.py','capture-source-snapshot.py','capture-frozen-build.py','run-frozen-root-gate.py','summarize-frozen-root-gates.py','run-resource-command.py']
for name in helpers:
 p=BASE/name;assert p.is_file();files.add(p)
# Capture the exact archival implementation and verification report when present.
for p in BASE.glob('*guard-point*'):
 if p.is_file() and p!=out:files.add(p)
for p in files:
 assert p.suffix!='.tmp' and '__pycache__' not in p.parts
entries={str(p.relative_to(ROOT)):{'sha256':digest(p),'size_bytes':p.stat().st_size} for p in sorted(files)}
external={}
def add_external(p):
 p=ROOT/p if not Path(p).is_absolute() else Path(p)
 if str(p.relative_to(ROOT)) not in entries:external[str(p.relative_to(ROOT))]={'sha256':digest(p),'size_bytes':p.stat().st_size}
for name in ['pre-refinement-artifacts.json','released-native-archives.json','released-native-archives-verification.json']:
 p=BASE/name
 if p.exists():add_external(p)
for name in ['docs/finite-density-conditional-points.md','docs/finite-density-plan.md']:
 add_external(name)
for path in [BASE/'guard-point-two-loop-fixed-comparison.json',BASE/'guard-point-two-loop-laurent-comparison.json',BASE/'guard-point-massive-regression-comparison.json',BASE/'guard-point-massive-default-path-comparison.json']:
 report=json.loads(path.read_text())
 for name,sha in report.get('inputs_sha256',{}).items():
  assert digest(ROOT/name)==sha
  add_external(name)
build=json.loads((BASE/'guard-point-build-provenance.json').read_text())
for name,record in build['files'].items():
 assert digest(ROOT/name)==record['sha256'],name
# The build provenance retains external executable hashes; generated target artifacts are not committed.
verify_snapshot(snapshot_path)
report={'schema':1,'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Completed guard-point regression checkpoint only; no full-feature or three/four-loop numerical acceptance. Historical failed zero-projection runs remain unchanged and linked through their frozen manifest and archive resolution.','status':'captured_and_verified','source_snapshot_sha256':digest(snapshot_path),'source_files_verified_unchanged':len(snapshot['files']),'build_provenance_sha256':digest(BASE/'guard-point-build-provenance.json'),'build_artifacts_verified_unchanged':len(build['files']),'validation_status_sha256':digest(BASE/'guard-point-validation-status.json'),'archive_manifest_sha256':digest(archive_path),'archives_roundtrip_verified':len(records),'file_count':len(entries),'size_bytes':sum(r['size_bytes'] for r in entries.values()),'files':entries,'external_references':external,'full_feature_validated':False,'three_loop_full_amplitude_validated':False,'native_four_loop_acceptance':False,'supplied_oracle_numerical_records_compared':0}
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','source_files_verified_unchanged','build_artifacts_verified_unchanged','archives_roundtrip_verified','file_count','size_bytes']}))
