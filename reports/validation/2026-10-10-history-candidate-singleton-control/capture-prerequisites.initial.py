"""Capture completed lower-loop comparisons only as launch evidence."""
import hashlib,json,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'2026-10-10-history-candidate-limit-integration'
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
cap=json.loads((HERE/'capsule-binding.json').read_text())
exe=next(x['sha256']for x in cap['static_files']if Path(x['path']).name=='runtime')
rows=[]
for name,refs,refines in [('two-loop-fixed',40,30),('two-loop-laurent',30,24),('massive-regression',10,0)]:
 stem='history-candidate-limit-'+name
 cp=PARENT/(stem+'-comparison.json');bp=PARENT/(stem+'-comparison-binding.json');rp=PARENT/(stem+'-binding.json')
 c=json.loads(cp.read_text());b=json.loads(bp.read_text());r=json.loads(rp.read_text())
 assert c['status']=='passed'and b['status']=='passed'and b['exit_code']==0 and b['comparison_sha256']==h(cp)
 assert c['reference_comparison_count']==refs and c['independent_refinement_comparison_count']==refines
 assert r['exit_code']==0 and r['post_run_source_and_executable_unchanged'] is True
 assert r['source_snapshot_sha256']==cap['source_snapshot_sha256'] and r['build_provenance_sha256']==cap['build_provenance_sha256'] and r['executable_sha256']==exe
 assert b['native_provenance']['source_snapshot_sha256']==cap['source_snapshot_sha256'] and b['native_provenance']['executable_sha256']==exe
 if name=='massive-regression':assert c['historical_same_profile_comparison_count']==10 and len(c['assembly_checks'])==2
 for p in [cp,bp,rp]:
  q=HERE/'evidence'/p.name;assert not q.exists();shutil.copy2(p,q);assert h(q)==h(p)
  rows.append({'source_path':str(p),'copy_path':str(q.relative_to(HERE)),'sha256':h(q),'bytes':q.stat().st_size})
record={'scope':'Completed lower-loop comparisons are launch prerequisites only; their values are not passed to the singleton prediction process.','captured_files':rows,'source_snapshot_sha256':cap['source_snapshot_sha256'],'build_provenance_sha256':cap['build_provenance_sha256'],'executable_sha256':exe,'capture_script_sha256':h(Path(__file__)),'all_three_comparisons_passed':True}
(HERE/'numerical-prerequisites.json').write_text(json.dumps(record,indent=2)+'\n')
print(h(HERE/'numerical-prerequisites.json'))
