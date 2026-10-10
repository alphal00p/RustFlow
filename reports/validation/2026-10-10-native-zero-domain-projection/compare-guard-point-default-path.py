"""Compare saved matching massive profiles before/after the guard-point change.

No integration or reference generation. Equality of decimal strings is reported;
acceptance also uses the existing fixed-dimension regression tolerances.
"""
import json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
sys.path.insert(0,str(ROOT/'tools/finite_density'))
from compare_complete_massive_reference import compare,require
from compare_massive_reference import complex_decimal,digest,read
old=BASE/'massive-regression';new=BASE/'guard-point-massive-regression'
paths=[old/'input.json',new/'input.json',old/'configuration.json',new/'configuration.json',old/'prediction-28-80-12.json',new/'prediction-28-80-12.json',Path(__file__).resolve()]
require(read(paths[0])==read(paths[1]),'changed massive input')
a,b=read(paths[2]),read(paths[3])
for value in (a,b):
 value['native_closure']=re.sub(r'checkpoints: Some\("[^\"]+"\)','checkpoints: Some("<report-local>")',value['native_closure'])
require(a==b,'default massive configuration changed beyond report path')
a,b=read(paths[4]),read(paths[5])
require(a['full_amplitude'] and b['full_amplitude'],'incomplete amplitude')
require(a['independent_reference_comparisons']==b['independent_reference_comparisons']==0,'predictions must be saved before comparison')
require(a['normalization']==b['normalization'],'normalization differs')
expected=[[],[0],[1],[0,1]]
require([x['cut_slots'] for x in a['contributions']]==[x['cut_slots'] for x in b['contributions']]==expected,'cut coverage changed')
old_values=[x['values'] for x in a['contributions']]+[a['values']]
new_values=[x['values'] for x in b['contributions']]+[b['values']]
rows=[]
for cut,previous,current in zip(expected+['full'],old_values,new_values):
 require(len(previous)==len(current)==2,'target coverage changed')
 for target,(before,after) in enumerate(zip(previous,current)):
  rows.append({'cut_slots':cut,'target_index':target,'before':before,'after':after,'identical_decimal_string':before==after,**compare(complex_decimal(after),complex_decimal(before),'sample')})
report={'schema':1,'status':'passed' if all(x['passed'] for x in rows) else 'failed','scope':'One same-profile massive default-closure regression across the guard-point source change, excluding only report-local checkpoint directory. This adds no precision-refinement claim.','numerical_comparisons':len(rows),'identical_decimal_strings':sum(x['identical_decimal_string'] for x in rows),'configuration_equal_after_report_path_normalization':True,'comparisons':rows,'inputs_sha256':{str(p.relative_to(ROOT)):digest(p) for p in paths}}
output=BASE/'guard-point-massive-default-path-comparison.json'
require(not output.exists(),'preserve earlier comparison')
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','numerical_comparisons','identical_decimal_strings']}))
if report['status']!='passed':raise SystemExit(1)
