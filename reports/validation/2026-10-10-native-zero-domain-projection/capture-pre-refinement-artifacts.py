"""Freeze the completed failed projection experiment before a future source fix.

All computations are complete. Source and report capture only; no test execution.
"""
import gzip,json
from datetime import datetime,timezone
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
output=BASE/'pre-refinement-artifacts.json'
if output.exists():raise ValueError('preserve completed checkpoint manifest')
source=BASE/'zero-projection-source-hashes.json';snapshot=verify_snapshot(source)
archive=json.loads((BASE/'zero-projection-three-loop-archives.json').read_text());assert archive['complete']
import hashlib
for row in archive['records']:
 p=ROOT/row['archive_path'];assert digest(p)==row['archive_sha256'];h=hashlib.sha256();length=0
 with gzip.open(p,'rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block);length+=len(block)
 assert h.hexdigest()==row['original_sha256'] and length==row['original_bytes']
 assert not (ROOT/row['original_path']).exists()
def inventory():
 return sorted(p for p in BASE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=output and not p.name.endswith('.tmp'))
paths=inventory();files={str(p.relative_to(ROOT)):{'sha256':digest(p),'size_bytes':p.stat().st_size} for p in paths}
external={'docs/finite-density-zero-source-projection.md'}
comparison=json.loads((BASE/'massive-regression-comparison.json').read_text())
for name,expected in comparison['inputs_sha256'].items():
 assert digest(ROOT/name)==expected
 if not (ROOT/name).is_relative_to(BASE):external.add(name)
old=json.loads((BASE/'priority-three-loop-fixed-old-build-binding.json').read_text())
for key in ('compiled_build','source_snapshot','baseline_command'):external.add(old[key]['path'])
external_files={name:{'sha256':digest(ROOT/name),'size_bytes':(ROOT/name).stat().st_size} for name in sorted(external)}
assert inventory()==paths
for p in paths:assert digest(p)==files[str(p.relative_to(ROOT))]['sha256']
verify_snapshot(source)
report={'schema':1,'captured_utc':datetime.now(timezone.utc).isoformat(),'scope':'Completed pre-point-refinement projection source/build/regression/numerical/diagnostic artifacts. No production point refinement is included. Future fixes require a new prefix and preserve these bytes.','status':'blocked_by_massless_condition_guard_regression','release_validated':False,'source_snapshot_sha256':digest(source),'source_files_verified_unchanged':len(snapshot['files']),'source_reconstruction':'pre-refinement-source-reconstruction.json','file_count':len(files),'total_bytes':sum(x['size_bytes'] for x in files.values()),'files':files,'external_artifacts':external_files,'archive_original_hashes_verified':archive['file_count'],'excluded':['Python __pycache__ bytecode','this self-referential manifest'],'passing_regression_tests':176,'ignored_regression_tests':1,'massive_reference_comparisons':10,'massive_historical_comparisons':10,'massless_fixed_variants_failed':2,'massless_laurent_not_run':True,'three_loop_numerical_predictions':0,'native_four_loop_acceptance':False,'source_capsule_reconstruction_verified':True}
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','file_count','total_bytes','archive_original_hashes_verified','source_files_verified_unchanged']}))
