"""Meaningful saved-data contract checks, without generating physical series values."""
from pathlib import Path
import copy,hashlib,importlib.util,json
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('saved_prefix',BASE/'sum-saved-prefix.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def bind(p):return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
plan={'producer_schema':'integrated-gaussian-rational-prefix.v1','dimension':'13/2','cut_slots':[0,3],'original_targets':[0,1],'profiles':{'known_terms':[2],'digits':[30]}}
data={'schema':plan['producer_schema'],'q':1,'variable':'z=1/eta','tail':'unknown','dimension':'13/2','input_identity':'fixture','known_prefix_length':2,'audit':[{'target':t,'full_off_shell_replay':True,'eta_leading_exponent_at_D':'1/4'} for t in [0,1]],'conditions':['eta'],'coefficients':[{'target':t,'n':n,'master':'channel','coefficient':'0' if n==1 else '-3/7'} for t in [0,1] for n in [0,1]]}
normal={'dimension':'13/2','input_identity':'fixture','cut_slots':[0,3],'profiles':[{'digits':30,'masters':[{'channel':'channel'}],'checked_conditions':[{'condition':'eta','value':{'re':'8','im':'0'}}]}]}
checks=[]
def case(name,change=None):
 d=copy.deepcopy(data);n=copy.deepcopy(normal);p=copy.deepcopy(plan)
 if change:change(d,n,p)
 try:m.validate(d,n,p);ok=change is None
 except (AssertionError,ValueError,ZeroDivisionError):ok=change is not None
 checks.append({'name':name,'passed':ok});assert ok,name
case('explicit zero coefficient admitted')
case('insufficient prefix rejected',lambda d,n,p:p['profiles'].update(known_terms=[3]))
case('unknown tail cannot be filled',lambda d,n,p:d['coefficients'].pop())
case('duplicate channel rejected',lambda d,n,p:d['coefficients'].append(copy.deepcopy(d['coefficients'][0])))
case('different physical input rejected',lambda d,n,p:n.update(input_identity='other'))
case('different cut set rejected',lambda d,n,p:n.update(cut_slots=[0]))
case('unbound master rejected',lambda d,n,p:n['profiles'][0]['masters'][0].update(channel='other'))
case('missing retained condition rejected',lambda d,n,p:n['profiles'][0].update(checked_conditions=[]))
case('vanished retained condition rejected',lambda d,n,p:n['profiles'][0]['checked_conditions'][0]['value'].update(re='0'))
case('mismatched leading branches rejected',lambda d,n,p:d['audit'].append({'target':0,'full_off_shell_replay':True,'eta_leading_exponent_at_D':'3/4'}))
case('unreplayed offshell factor rejected',lambda d,n,p:d['audit'][0].update(full_off_shell_replay=False))
case('expression coefficient rejected',lambda d,n,p:d['coefficients'][0].update(coefficient='1+epsilon'))
actual=json.loads((BASE/'normalization-24/coefficients.json').read_text());native=json.loads((BASE/'normalization-24/normalization.json').read_text());prospective=json.loads((BASE/'plan.json').read_text())
try:m.validate(actual,native,prospective);raise RuntimeError('24 coefficients incorrectly supplied planned 64')
except AssertionError as e:checks.append({'name':'actual frozen 24-prefix refuses planned 32/48/64','passed':'insufficient saved prefix' in str(e)});assert checks[-1]['passed']
result={'schema':1,'status':'passed','checks':checks,'physical_partial_sums_evaluated':False,'inputs':[bind(BASE/'sum-saved-prefix.py'),bind(Path(__file__)),bind(BASE/'plan.json'),bind(BASE/'normalization-24/coefficients.json'),bind(BASE/'normalization-24/normalization.json')]}
with (BASE/'adapter-controls.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':len(checks),'physical_partial_sums_evaluated':False}))
