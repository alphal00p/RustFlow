"""Post-save diagnostic comparisons; no equation or coefficient generation."""
from pathlib import Path
import hashlib,json,itertools,re
import mpmath as mp
B=Path(__file__).resolve().parent

def bind(p):p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def cval(x):return mp.mpc(x['re'],x['im'])
def legacy(x):
 m=re.fullmatch(r'\(([-+0-9.eE]+)\+0i\)',x);assert m;return mp.mpc(m.group(1))
def main():
 assert not (B/'comparison.json').exists();plan=json.loads((B/'plan.json').read_text());run=json.loads((B/'prediction-run.json').read_text());assert run['exit_code']==0 and run['inputs_unchanged']
 assert bind(B/'predictions.json')==run['prediction'];prediction=json.loads((B/'predictions.json').read_text());assert prediction['reference_reads']==0
 for v in prediction['inputs'].values():assert bind(v['path'])==v
 inputs={n:bind(B/n) for n in ['plan.json','predictions.json','prediction-run.json','compare.py']}
 rows={};grid=set()
 for target,eta,digits in itertools.product(plan['producer']['targets'],plan['profiles']['eta'],plan['profiles']['digits']):
  for n in plan['profiles']['partial_sum_terms']:grid.add(('partial_sum',target,eta,digits,str(n)))
  for l,m in plan['profiles']['pade_degrees']:grid.add(('pade',target,eta,digits,f'{l}/{m}'))
 for row in prediction['rows']:
  degree=str(row['terms']) if row['method']=='partial_sum' else '/'.join(map(str,row['degrees']))
  key=(row['method'],row['target'],row['eta'],row['digits'],degree);assert key not in rows;rows[key]=row
 assert set(rows)<=grid
 if not prediction['failures']:assert set(rows)==grid
 mp.mp.dps=100;refs={}
 for eta in ['8','2']:
  record=plan['references']['eta'+eta];assert bind(record['path'])==record;inputs['eta'+eta]=record
  ref=json.loads(Path(record['path']).read_text());assert ref['profile']=={'dimension':'6.5','mu':'1','eta':eta,'nodes':32,'digits':80}
  assert ref['values']=='raw Euclidean two-cut contribution'
  refs[0,eta]=mp.mpc(ref['scalar']);refs[1,eta]=mp.mpc(ref['raised'])
 record=plan['references']['eta0'];assert bind(record['path'])==record;inputs['eta0']=record;ref=json.loads(Path(record['path']).read_text())
 assert ref['normalization']=='raw Euclidean' and type(ref['native_amf_predictions_read']) is int and ref['native_amf_predictions_read']==0 and type(ref['oracle_answers_read']) is int and ref['oracle_answers_read']==0
 samples=[x for x in ref['samples'] if x['epsilon']=='-5/4' and x['digits']==50];assert len(samples)==1;sample=samples[0]
 assert abs(legacy(sample['dimension'])-mp.mpf('6.5'))==0
 assert {tuple(x['cuts']) for x in sample['contributions']}=={(),(0,),(3,),(0,3)}
 values=[x['values'] for x in sample['contributions'] if x['cuts']==[0,3]][0]
 for target in plan['producer']['targets']:refs[target,'0']=legacy(values[target])
 criterion=mp.mpf(plan['prospective_numerical_criterion']['relative']);zero=mp.mpf(plan['prospective_numerical_criterion']['absolute_zero'])
 def check(a,b):
  error=abs(a-b);limit=criterion*abs(b) if b else zero
  return {'relative_error':mp.nstr(error/abs(b),40) if b else None,'absolute_error':mp.nstr(error,40),'passed':bool(error<=limit),'limit':mp.nstr(limit,40)}
 checks=[];precision=[];refinement=[]
 for key,row in sorted(rows.items()):
  method,target,eta,digits,degree=key;checks.append({'method':method,'target':target,'eta':eta,'digits':digits,'profile':degree,**check(cval(row['value']),refs[target,eta])})
 for method,target,eta,_,degree in sorted(k for k in rows if k[3]==30):
  if (method,target,eta,50,degree) in rows:precision.append({'method':method,'target':target,'eta':eta,'profile':degree,'digits':[30,50],**check(cval(rows[method,target,eta,30,degree]['value']),cval(rows[method,target,eta,50,degree]['value']))})
 for method,profile in [('partial_sum',['32','48','64']),('pade',['16/16','24/24','31/32'])]:
  for target,eta in itertools.product(plan['producer']['targets'],plan['profiles']['eta']):
   for lo,hi in zip(profile,profile[1:]):
    a=(method,target,eta,50,lo);b=(method,target,eta,50,hi)
    if a in rows and b in rows:refinement.append({'method':method,'target':target,'eta':eta,'profiles':[lo,hi],**check(cval(rows[a]['value']),cval(rows[b]['value']))})
 result={'schema':1,'status':'completed diagnostic comparison; no endpoint certification','prediction_binding':run['prediction'],'inputs':inputs,'reference_checks':checks,'precision_checks':precision,'profile_refinement_checks':refinement,'failed_pade_attempts':prediction['failures'],'endpoint_acceptance':False,'ode_identity_asserted':False,'full_three_or_four_loop_acceptance':False,'error_scope':'Empirical discrepancies and stability only; source-only scalar tail bound is a separate artifact'}
 assert all(bind(x['path'])==x for x in inputs.values())
 with (B/'comparison.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'final_profiles':[x for x in checks if x['digits']==50 and x['profile'] in ['64','31/32']],'precision_passed':sum(x['passed'] for x in precision),'precision_checks':len(precision),'failures':len(prediction['failures'])}))
if __name__=='__main__':main()
