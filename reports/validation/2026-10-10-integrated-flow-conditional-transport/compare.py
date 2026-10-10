from pathlib import Path
from decimal import Decimal as D,localcontext
import json,hashlib
R=Path(__file__).resolve().parent
with localcontext() as ctx:
 ctx.prec=100
 finite=json.loads((R/'finite-result.json').read_text());rows=finite['records'];errs=[r for r in rows if 'error'in r]
 final={r['end_eta']:r for r in rows if r.get('digits')==50 and r.get('boundary_terms')==64 and r.get('start_eta')==32 and 'state'in r}
 diffs=[]
 for r in rows:
  if 'state'not in r:continue
  b=final[r['end_eta']];x=D(r['state'][0]['re']);y=D(b['state'][0]['re']);error=abs(x-y);rel=error/max(abs(y),D('1e-100'))
  diffs.append({'digits':r['digits'],'terms':r['boundary_terms'],'start':r['start_eta'],'eta':r['end_eta'],'absolute':str(error),'relative':str(rel)})
 record={'scope':'Candidate-internal refinements only, no physical reference values','identity_certified':False,'finite_profiles':len(rows)//5,'finite_records':len(rows),'errors':errs,'max_relative_refinement':str(max(D(r['relative'])for r in diffs)),'all_relative_below_1e_minus12':all(D(r['relative'])<D('1e-12')for r in diffs),'finite_differences':diffs}
 if (R/'endpoint-result.json').exists():
  end=json.loads((R/'endpoint-result.json').read_text()); er=end.get('records',[]);success=[r for r in er if 'endpoint'in r.get('result',{})]
  best=[r for r in success if r['digits']==50 and r['boundary_terms']==64 and r['start_eta']==32 and r['match_eta']=='1/8']
  record['endpoint_status']=end['status'];record['endpoint_errors']=[r for r in er if 'error'in r.get('result',{}) or 'error'in r];record['endpoint_records']=len(er);record['endpoint_successes']=len(success)
  if best:
   y=D(best[0]['result']['endpoint'][0]['re']);ed=[]
   for r in success:
    x=D(r['result']['endpoint'][0]['re']);error=abs(x-y);rel=error/max(abs(y),D('1e-100'));ed.append({'digits':r['digits'],'terms':r['boundary_terms'],'start':r['start_eta'],'match':r['match_eta'],'absolute':str(error),'relative':str(rel)})
   record['candidate_endpoint_normalized']=str(y);record['endpoint_max_relative_refinement']=str(max(D(r['relative'])for r in ed));record['endpoint_differences']=ed
 (R/'comparison.json').write_text(json.dumps(record,indent=2)+'\n')
 print({k:v for k,v in record.items()if not k.endswith('differences')})
