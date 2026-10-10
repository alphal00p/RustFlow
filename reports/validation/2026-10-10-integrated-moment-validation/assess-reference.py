"""Assess saved reference refinements only; never read producer values."""
import argparse, hashlib, json
from decimal import Decimal, localcontext
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def values(d):return {'scalar':d['scalar'],'raised':d['raised'],**d['raised_components']}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('runs',nargs=4,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 assert not a.output.exists()
 by={};inputs=[]
 for run in a.runs:
  result=json.loads((run/'reference.json').read_text());res=json.loads((run/'resources.json').read_text())
  assert res['exit_code']==0 and res['bound_inputs_unchanged'] and res['within_individual_budget']
  assert res['reference_sha256']==sha(run/'reference.json')
  assert result['receipt_sha256']==sha(run/'receipt.json')
  receipt=json.loads((run/'receipt.json').read_text());assert receipt['profile']==result['profile']
  profile=result['profile'];key=(profile['nodes'],profile['digits']);assert key not in by;by[key]=result
  inputs.append({'run':str(run.resolve()),'reference_sha256':sha(run/'reference.json'),'resources_sha256':sha(run/'resources.json'),'receipt_sha256':sha(run/'receipt.json')})
 assert set(by)=={(16,50),(24,50),(32,50),(32,80)}
 profiles=list(by.values());assert len({(x['profile']['dimension'],x['profile']['mu'],x['profile']['eta']) for x in profiles})==1
 comparisons=[]
 with localcontext() as ctx:
  ctx.prec=100
  for left,right,kind in [((16,50),(24,50),'node refinement 16 to 24'),((24,50),(32,50),'node refinement 24 to 32'),((32,50),(32,80),'precision refinement 50 to 80')]:
   av,bv=values(by[left]),values(by[right]);assert set(av)==set(bv)
   for name in av:
    aa,bb=Decimal(av[name]),Decimal(bv[name]);assert aa.is_finite() and bb.is_finite()
    error=abs(aa-bb);scale=max(abs(aa),abs(bb))
    small=scale<Decimal('1e-20');limit=Decimal('1e-27') if small else Decimal('1e-14')*scale
    comparisons.append({'kind':kind,'component':name,'absolute_error':str(error),'relative_error':str(error/scale) if scale else '0','limit':str(limit),'passed':error<=limit})
  kernel=[]
  for sample in by[(32,80)]['parameter_cross_checks']:
   for name in ['value_relative_error','derivative_relative_error']:
    error=Decimal(sample[name]);kernel.append({'h':sample['h'],'kind':name,'error':str(error),'passed':error<=Decimal('1e-14')})
  result={'schema':1,'status':'empirical reference refinement passed' if all(x['passed'] for x in comparisons+kernel) else 'reference accuracy unresolved',
   'profile':by[(32,80)]['profile'],'inputs':inputs,'assessor_sha256':sha(__file__),
   'refinement_checks':comparisons,'independent_parameter_checks':kernel,
   'rigorous_error_bound':False,'native_or_producer_comparison_performed':False,
   'criterion':'14 relative digits; absolute 1e-27 when both magnitudes below 1e-20; independently varied node count and arithmetic precision',
   'scope':'scalar, raised total and all four raised components; two-cut raw Euclidean reference only'}
  a.output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
