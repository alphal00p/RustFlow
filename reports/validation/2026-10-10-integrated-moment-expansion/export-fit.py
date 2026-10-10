from pathlib import Path
from fractions import Fraction
import hashlib,json
R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=R/'gaussian-wide-prefix-64/coefficients.json';d=json.loads(p.read_text());old=json.loads((R/'gaussian-wide-prefix-24/coefficients.json').read_text())
assert d['known_prefix_length']==64 and d['dimension']=='13/2' and d['q']==1
assert [x for x in d['coefficients'] if x['n']<24]==old['coefficients']
contract={'schema':'integrated-moment-series-contract.v1','scope':'actual cut[0,3], all uncut shifted, original fixed-index scalar and raised mixed-medium targets','input_identity':d['input_identity'],'dimension':'13/2','q':1,'variable':'z=1/eta','beta':'-1/4','log_power':0,'normalization':'A(D)^2 * native one-loop unit-mass Minkowski tadpole I[1], A=1/(sqrt(pi)*Gamma((D-1)/2)), mu_i=1; native physical target coefficients included','all_regions':'census enumerates two; soft virtual polynomial is ordinary native ProvedZero, all-hard Gaussian channel remains','original_virtual_numerator_degree_max':1,'known_prefix_length':64,'training_count':48,'holdout_count':16,'tail':'unknown','reference_values_used':False,'sources':{str(x.relative_to(R)):sha(x) for x in [p,R/'census.json',R/'grading.md',R/'gaussian-method.md',R/'gaussian-wide-build.json',R/'gaussian-four-prefix-comparison.json',R/'wide-24-comparison.json']}}
out=R/'fit-data';out.mkdir(exist_ok=False);(out/'series-contract.json').write_text(json.dumps(contract,indent=2)+'\n');proof=sha(R/'grading.md')
config={'schema':'integrated-per-target-fit-plan.v1','targets':[]}
for ti in [0,1]:
 rows=[x for x in d['coefficients'] if x['target']==ti];assert len(rows)==64
 masters=sorted({x['master'] for x in rows});assert len(masters)==1
 series=hashlib.sha256((sha(out/'series-contract.json')+':target='+str(ti)).encode()).hexdigest()
 for name,start,end in [('training',0,48),('holdout',48,64)]:
  channels=[]
  for master in masters:
   by={x['n']:str(Fraction(x['coefficient'])) for x in rows if x['master']==master};assert set(by)==set(range(64))
   channels.append({'id':hashlib.sha256(master.encode()).hexdigest(),'beta':'-1/4','log_power':0,'coefficients':[by[n] for n in range(start,end)]})
  x={'schema':'integrated-rational-prefix.v1','dimension':'13/2','q':1,'series_identity':series,'grading_proof_identity':proof,'target_index':ti,'first_index':start,'channels':channels}
  (out/f'target-{ti}-{name}.json').write_text(json.dumps(x,indent=2)+'\n')
 config['targets'].append({'id':ti,'training':str(out/f'target-{ti}-training.json'),'holdout':str(out/f'target-{ti}-holdout.json')})
(out/'fit-plan.json').write_text(json.dumps(config,indent=2)+'\n')
manifest={'scope':'Producer generated all64 before training; fitter must not open holdout before both candidates frozen','files':{x.name:sha(x) for x in sorted(out.glob('*.json'))}}
(out/'producer-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'status':'PASS','coefficients':128,'first24_exact_match':True,'fit_plan':str(out/'fit-plan.json'),'series_contract_sha256':sha(out/'series-contract.json')}))
