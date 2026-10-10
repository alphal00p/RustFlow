"""Normalize and save frozen conditional ODE predictions before reading references."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,re,shutil,sys
import mpmath as mp
B=Path(__file__).resolve().parent

def bind(p):p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def load(record):assert bind(record['path'])==record;return json.loads(Path(record['path']).read_text())
def cval(x):return mp.mpc(x['re'],x['im'])
def printable(x):return {'re':mp.nstr(mp.re(x),mp.mp.dps),'im':mp.nstr(mp.im(x),mp.mp.dps)}
def legacy(s):m=re.fullmatch(r'\(([-+0-9.eE]+)\+0i\)',s);assert m;return mp.mpc(m.group(1))

def main():
 assert not (B/'comparison.json').exists();plan=json.loads((B/'plan.json').read_text());inputs={'plan':bind(B/'plan.json'),'script':bind(Path(__file__))}
 for name,item in plan['producer'].items():assert bind(item['path'])==item;inputs[name]=item
 source={name:json.loads(Path(item['path']).read_text()) for name,item in plan['producer'].items() if name.endswith('.json')}
 for run in ['finite-run.json','endpoint-run.json']:assert source[run]['exit_code']==0
 normal=load(plan['normalization']);normal_run=load(plan['normalization_run']);normal_build=load(plan['normalizer_build']);coef=load(plan['moment_coefficients']);contract=load(plan['series_contract'])
 assert normal_run['exit_code']==0 and normal_build['exit_code']==0 and normal_run['output']==plan['normalization']
 assert normal['cut_slots']==[0,3] and normal['dimension']=='13/2' and normal['input_identity']==coef['input_identity']==contract['input_identity']
 assert contract['beta']=='-1/4' and contract['tail']=='unknown' and contract['reference_values_used'] is False
 frozen_coefs=source['scalar-coefficients.json'];assert frozen_coefs['beta']=='-1/4'
 scalar=sorted([x for x in coef['coefficients'] if x['target']==0],key=lambda x:x['n'])
 assert [x['n'] for x in scalar]==list(range(64)) and len({x['master'] for x in scalar})==1
 assert list(map(Q,frozen_coefs['coefficients']))==[Q(x['coefficient']) for x in scalar]
 master=scalar[0]['master'];companion=source['companion.json'];assert companion['identity_certified'] is False and companion['target_index']==0 and companion['beta']=='-1/4'
 seen=set()
 for build_name,source_name in [('build.json','transport.rs'),('endpoint-build.json','endpoint.rs')]:
  build=source[build_name];assert build['exit_code']==0 and build['source_sha256']==plan['producer'][source_name]['sha256']
  for name,item in build['dependencies'].items():assert item['sha256']==normal_build['dependencies'][name]['sha256']
  for item in [build['executable'],*build['dependencies'].values()]:
   if item['path'] not in seen:assert bind(item['path'])['sha256']==item['sha256'];seen.add(item['path'])
 finite=source['finite-result.json'];endpoint=source['endpoint-result.json']
 assert finite['identity_certified'] is False and finite['reference_values_read'] is False and finite['target_index']==0
 assert endpoint['identity_certified'] is False and endpoint['physical_endpoint_accepted'] is False
 expected=set(itertools.product([30,50],[32,48,64],[16,32],[16,8,4,2,1]));found=set();rows=[]
 profiles={x['digits']:x for x in normal['profiles']};mp.mp.dps=90
 def factor(digits):
  p=profiles[digits];assert len(p['masters'])==1 and p['masters'][0]['channel']==master
  return cval(p['native_measure_to_euclidean'])*cval(p['compact_normalization'])*cval(p['masters'][0]['value'])
 for row in finite['records']:
  key=(row['digits'],row['boundary_terms'],row['start_eta'],row['end_eta']);assert key not in found;found.add(key)
  assert row['working_digits']==row['digits']+20 and len(row['state'])==5
  rows.append({'kind':'finite','digits':row['digits'],'boundary_terms':row['boundary_terms'],'start_eta':row['start_eta'],'eta':str(row['end_eta']),'value':printable(factor(row['digits'])*cval(row['state'][0]))})
 assert found==expected
 expected=set(itertools.product([30,50],[32,48,64],[16,32]));found=set();failures=[]
 for row in endpoint['records']:
  assert row['working_digits']==row['digits']+20 and row['frobenius_terms']==row['boundary_terms']
  if row['match_eta']=='1/8':
   key=(row['digits'],row['boundary_terms'],row['start_eta']);assert key not in found;found.add(key)
   assert 'error' not in row['result'] and len(row['result']['endpoint'])==5
   rows.append({'kind':'endpoint','digits':row['digits'],'boundary_terms':row['boundary_terms'],'start_eta':row['start_eta'],'frobenius_terms':row['frobenius_terms'],'match_eta':row['match_eta'],'eta':'0','value':printable(factor(row['digits'])*cval(row['result']['endpoint'][0]))})
  else:
   assert row['match_eta']=='1/16' and row['digits']==50 and row['boundary_terms']==64 and row['start_eta']==32 and 'error' in row['result'];failures.append(row)
 assert found==expected and len(failures)==1
 for k in ['normalization','normalization_run','normalizer_build','moment_coefficients','series_contract']:inputs[k]=plan[k]
 physical={'schema':1,'scope':'Normalized saved scalar candidate predictions; equation remains unproved','inputs':inputs,'rows':rows,'failed_extra_endpoint_match':failures,'reference_reads':0,'identity_certified':False}
 with (B/'physical-predictions.json').open('x') as f:f.write(json.dumps(physical,indent=2)+'\n')
 prediction_binding=bind(B/'physical-predictions.json')
 # Numerical references are opened only after physical predictions have been saved.
 eta2=load(plan['references']['eta2']);eta0=load(plan['references']['eta0'])
 assert eta2['profile']=={'dimension':'6.5','mu':'1','eta':'2','nodes':32,'digits':80} and eta2['values']=='raw Euclidean two-cut contribution'
 assert eta0['normalization']=='raw Euclidean' and type(eta0['native_amf_predictions_read']) is int and eta0['native_amf_predictions_read']==0 and eta0['oracle_answers_read']==0
 sample=[x for x in eta0['samples'] if x['epsilon']=='-5/4' and x['digits']==50];assert len(sample)==1
 cut=[x for x in sample[0]['contributions'] if x['cuts']==[0,3]];assert len(cut)==1
 refs={'2':mp.mpc(eta2['scalar']),'0':legacy(cut[0]['values'][0])}
 rel=mp.mpf(plan['criterion']['relative']);zero=mp.mpf(plan['criterion']['absolute_zero']);checks=[]
 for row in rows:
  if row['eta'] not in refs:continue
  expected=refs[row['eta']];actual=cval(row['value']);error=abs(actual-expected);limit=rel*abs(expected) if expected else zero
  checks.append({k:v for k,v in row.items() if k!='value'}|{'absolute_error':mp.nstr(error,45),'relative_error':mp.nstr(error/abs(expected),45),'passed':bool(error<=limit),'criterion_limit':mp.nstr(limit,45)})
 result={'schema':1,'status':'all conditional scalar comparisons passed' if all(x['passed'] for x in checks) else 'conditional scalar comparison failure','physical_predictions':prediction_binding,'reference_bindings':plan['references'],'checks':checks,'failed_extra_endpoint_match':failures,'scope':'Numerical evidence conditional on an unproved candidate equation; finite-prefix fit is not an all-order identity','identity_certified':False,'physical_endpoint_certified':False,'raised_result':False,'full_amplitude_acceptance':False}
 assert bind(B/'physical-predictions.json')==prediction_binding
 for item in inputs.values():assert bind(item['path'])==item
 with (B/'comparison.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 copies=B/'producer-freeze';copies.mkdir()
 for name,item in plan['producer'].items():shutil.copyfile(item['path'],copies/name)
 print(json.dumps({'status':result['status'],'checks':len(checks),'passed':sum(x['passed'] for x in checks),'maximum_relative_error':str(max(mp.mpf(x['relative_error']) for x in checks)),'failed_extra_endpoint_preserved':len(failures),'identity_certified':False}))
if __name__=='__main__':main()
