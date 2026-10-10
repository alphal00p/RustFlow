"""Compare saved one-sided original-mass differences to separate raised pieces."""
import argparse,hashlib,json
from decimal import Decimal,localcontext
from pathlib import Path
BASE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(tag):
 folder=BASE/'runs'/tag
 result=json.loads((folder/'reference.json').read_text());resources=json.loads((folder/'resources.json').read_text())
 assert resources['exit_code']==0 and resources['bound_inputs_unchanged'] and resources['within_individual_budget']
 assert resources['reference_sha256']==sha(folder/'reference.json')
 assert result['receipt_sha256']==sha(folder/'receipt.json')
 return result,{'tag':tag,'reference_sha256':sha(folder/'reference.json'),'resources_sha256':sha(folder/'resources.json')}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args();assert not args.output.exists()
 plan=json.loads((BASE/'prospective-mass-check.json').read_text())
 analytic_path=BASE/'runs/eta2-n32-d80/reference.json';assert sha(analytic_path)==plan['reference_baseline_sha256']
 analytic=json.loads(analytic_path.read_text());inputs=[];rows=[]
 with localcontext() as ctx:
  ctx.prec=100
  components={k:Decimal(v) for k,v in analytic['raised_components'].items()}
  exact={'bulk':sum(v for k,v in components.items() if k!='upper_surface'),'upper_surface':components['upper_surface'],'total':Decimal(analytic['raised'])}
  computed={}
  for n in [24,32]:
   zero,binding=load(f'mass-zero-n{n}-d50');inputs.append(binding)
   J0=Decimal(zero['fixed_radius_unraised_original_numerator']);assert J0==Decimal(zero['moving_radius_unraised_original_numerator'])
   for denominator in ([256,1024,4096] if n==24 else [4096]):
    value,binding=load(f'mass-a{denominator}-n{n}-d50');inputs.append(binding)
    assert value['profile']['physical_mass_squared']==f'1/{denominator}'
    fixed=Decimal(value['fixed_radius_unraised_original_numerator']);moving=Decimal(value['moving_radius_unraised_original_numerator'])
    values={'bulk':(J0-fixed)*denominator,'upper_surface':(fixed-moving)*denominator,'total':(J0-moving)*denominator}
    computed[n,denominator]=values
    for name,v in values.items():
     assert exact[name]!=0
     relative=abs(v/exact[name]-1)
     rows.append({'nodes':n,'mass_squared':f'1/{denominator}','component':name,'finite_difference':str(v),'analytic_mass_jet':str(exact[name]),'relative_error':str(relative)})
  checks=[]
  for name in exact:
   errors=[abs(computed[24,d][name]/exact[name]-1) for d in [256,1024,4096]]
   checks.append({'component':name,'kind':'error decreases as positive mass decreases','errors':[str(e) for e in errors],'passed':errors[2]<errors[1]<errors[0]})
   final=abs(computed[32,4096][name]/exact[name]-1)
   checks.append({'component':name,'kind':'final relative mass error','error':str(final),'limit':plan['criteria']['final_mass_relative_error_max'],'passed':final<=Decimal(plan['criteria']['final_mass_relative_error_max'])})
   refinement=abs((computed[32,4096][name]-computed[24,4096][name])/exact[name])
   checks.append({'component':name,'kind':'independent node refinement at smallest mass','error':str(refinement),'limit':plan['criteria']['smallest_mass_node24_to32_relative_change_max'],'passed':refinement<=Decimal(plan['criteria']['smallest_mass_node24_to32_relative_change_max'])})
  report={'schema':1,'status':'physical mass and Fermi-surface diagnostic passed' if all(c['passed'] for c in checks) else 'physical mass diagnostic unresolved','scope':'Independent finite-difference derivative/sign check only; not a 14-digit amplitude gate','prospective_plan_sha256':sha(BASE/'prospective-mass-check.json'),'analytic_reference_sha256':sha(analytic_path),'source_sha256':sha(__file__),'inputs':inputs,'comparisons':rows,'checks':checks,'producer_coefficients_read':False,'integer_Richardson_assumption':False,'rigorous_error_bound':False}
  args.output.write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
