"""Post-save comparison only; this tool never writes or changes producer coefficients."""
from pathlib import Path
import hashlib,json,itertools
import mpmath as mp
BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2];REF=ROOT/'reports/validation/2026-10-10-integrated-moment-validation'
def bind(p):
 p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def val(v):return mp.mpc(v['re'],v['im'])
def main():
 prediction=BASE/'series-predictions.json';planpath=BASE/'plan.json';out=BASE/'comparison.json';assert not out.exists()
 # Authenticate complete predictions and their source bindings before opening references.
 pred=json.loads(prediction.read_text());plan=json.loads(planpath.read_text())
 assert pred['numerical_reference_reads']==0 and pred['unknown_tail'] is True
 assert pred['cut_slots']==plan['cut_slots'] and pred['dimension']==plan['dimension']
 for b in pred['inputs'].values():assert bind(b['path'])==b
 grid=set(itertools.product(plan['original_targets'],plan['profiles']['eta'],plan['profiles']['known_terms'],plan['profiles']['digits']))
 rows={}
 for r in pred['rows']:
  key=(r['target'],r['eta'],r['terms'],r['digits']);assert key not in rows;rows[key]=r
 assert set(rows)==grid and pred['known_prefix_length']>=max(plan['profiles']['known_terms'])
 inputs={'prediction':bind(prediction),'plan':bind(planpath),'comparator':bind(Path(__file__))}
 manifest_path=REF/'artifact-manifest.json';assert bind(manifest_path)==plan['reference_manifest']
 manifest=json.loads(manifest_path.read_text());reference_files={x['path']:x for x in manifest['files']}
 def saved(name):
  p=REF/name;b=bind(p);expected=reference_files[name]
  assert b['sha256']==expected['sha256'] and b['bytes']==expected['bytes'];inputs[name]=b;return json.loads(p.read_text())
 mp.mp.dps=100;reference={};assessments={}
 for eta in plan['profiles']['eta']:
  assessments[eta]=saved(f'eta{eta}-refinement.json');assert assessments[eta]['status']=='empirical reference refinement passed'
  assert all(x['passed'] for x in assessments[eta]['refinement_checks']+assessments[eta]['independent_parameter_checks'])
  r=saved(f'runs/eta{eta}-n32-d80/reference.json');assert r['profile']=={'dimension':'6.5','mu':'1','eta':eta,'nodes':32,'digits':80}
  assert r['values']=='raw Euclidean two-cut contribution' and r['empirical_accuracy_only'] is True
  reference[(0,eta)]=mp.mpc(r['scalar']);reference[(1,eta)]=mp.mpc(r['raised'])
 rel=mp.mpf(plan['numerical_checks']['relative_criterion']);zero=mp.mpf(plan['numerical_checks']['absolute_zero_criterion'])
 def check(actual,expected):
  error=abs(actual-expected);denom=abs(expected);limit=zero if not denom else rel*denom
  return {'absolute_error':mp.nstr(error,40),'relative_error':mp.nstr(error/denom,40) if denom else None,'criterion_limit':mp.nstr(limit,40),'passed':bool(error<=limit)}
 refs=[];precisions=[];truncations=[]
 for key,r in sorted(rows.items()):
  target,eta,terms,digits=key;refs.append({'target':target,'eta':eta,'terms':terms,'digits':digits,'scope':'eta 2 raw sums are diagnostic outside the initial disk' if eta=='2' else 'empirical finite-eta comparison',**check(val(r['value']),reference[(target,eta)])})
 for target,eta,terms in itertools.product(plan['original_targets'],plan['profiles']['eta'],plan['profiles']['known_terms']):
  lo,hi=plan['profiles']['digits'];precisions.append({'target':target,'eta':eta,'terms':terms,'digits':[lo,hi],**check(val(rows[(target,eta,terms,lo)]['value']),val(rows[(target,eta,terms,hi)]['value']))})
 terms=plan['profiles']['known_terms'];digits=max(plan['profiles']['digits'])
 for target,eta in itertools.product(plan['original_targets'],plan['profiles']['eta']):
  for lo,hi in zip(terms,terms[1:]):truncations.append({'target':target,'eta':eta,'terms':[lo,hi],'digits':digits,**check(val(rows[(target,eta,lo,digits)]['value']),val(rows[(target,eta,hi,digits)]['value']))})
 final=[r for r in refs if r['eta']=='8' and r['terms']==max(terms) and r['digits']==digits]
 final_refinement=[r for r in truncations if r['eta']=='8' and r['terms']==terms[-2:]]
 result={'schema':1,'status':'eta8 final partial sums match empirical reference and final term refinement' if all(x['passed'] for x in final+final_refinement) else 'eta8 final partial sums did not establish the prospective empirical criterion','reference_checks':refs,'precision_checks':precisions,'truncation_checks':truncations,'inputs':inputs,'reference_manifest':bind(manifest_path),'numerical_acceptance':{'scope':'finite eta8 empirical consistency of saved source-derived coefficients only','eta2':'diagnostic raw-series extrapolation, never continuation acceptance','rigorous_error_bound':False,'endpoint':False,'ode_identity':False,'three_or_four_loop_full_amplitude_acceptance':False},'source_germ_coverage':'Separate producer proof required; this comparator does not infer analyticity from numeric agreement'}
 assert all(bind(v['path'])==v for v in inputs.values())
 with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'eta8_final':final,'eta8_final_refinement':final_refinement,'precision_passed':sum(x['passed'] for x in precisions),'precision_checks':len(precisions),'eta2_raw_reference_passed':sum(x['passed'] for x in refs if x['eta']=='2')}))
if __name__=='__main__':main()
