"""Exact audit arithmetic and source binding; no fits or numerical profiles."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib, json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'reports/validation'
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def put(name,value): (HERE/name).write_text(json.dumps(value,indent=2)+'\n')
def choose(a,n):
    out=Q(1)
    for j in range(n): out*= (a-j)/(j+1)
    return out
m=lambda e:1/(Q(9,2)+e)
checks=[]
def check(name,actual,expected):
    assert actual==expected,(name,actual,expected)
    checks.append({'name':name,'actual':str(actual),'expected':str(expected),'pass':True})
bulk=Q(5,2)*(m(0)**2+m(-1)*m(1))
surface=(m(0)+3*m(1))/2
derivative=4*(m(2)*m(0)+m(1)**2)
check('bulk_absolute',bulk,Q(1580,6237))
check('surface_absolute',surface,Q(38,99))
check('kernel_derivative_absolute',derivative,Q(3808,14157))
check('B0_multiplier',Q(9,8)*(bulk+surface),Q(1987,2772))
check('B1_multiplier',Q(9,8)*derivative,Q(476,1573))
check('signed_leading_raised',Q(9,8)*(Q(3,2)*m(0)**2-m(0)/2+m(1)),Q(43,264))
# Triangular identity on a deterministic synthetic sequence; no holdout opened.
a=Q(1,4);R=Q(4)
g=[Q((-1)**n*(n+1),n+2) for n in range(16)]
h=[sum((g[n]/R**n*(-1)**(k-n)*choose(a-n,k-n) for n in range(k+1)),Q(0)) for k in range(16)]
for k in range(16):
    check('composition_inverse_'+str(k),R**k*sum((h[n]*choose(a-n,k-n) for n in range(k+1)),Q(0)),g[k])
endpoint=BASE/'2026-10-10-integrated-flow-conditional-transport/endpoint-result.json'
e=json.loads(endpoint.read_text())
assert e['identity_certified'] is False and e['physical_endpoint_accepted'] is False
assert e['exponents']==['0','1/4','5/2','7/2','19/4']
records=e['records']; success=[r for r in records if 'endpoint' in r.get('result',{})]
assert len(success)==12
assert all(r['match_eta']=='1/8' for r in success)
failed=[r for r in records if r not in success]
assert len(failed)==1 and failed[0]['match_eta']=='1/16'
put('exact-checks.json',{'schema':1,'scope':'rational arithmetic and frozen metadata only; no new coefficient generation, fit, reference or numerical profile','rational_checks':checks,'conditional_endpoint_metadata':{'successful_profiles':len(success),'failed_profiles':len(failed),'physical_endpoint_accepted':False,'identity_certified':False}})
sources={
'producer':['artifact-manifest.json','scalar-source-kernel.json','grading.md','convergent-germ.md','sharp-germ-bound.md','generation-summary.json','gaussian-wide-source/generate.rs','gaussian-wide-source/compact_exact.rs'],
'continuation':['plan.json','predict.py','predictions.json','transformed-coefficients.json','prediction-run.json','bound-plan.json'],
'candidate':['target-0/recurrence-candidate.json','certificate-pole-check.py','certificate-pole-check.json'],
'conditional':['derive.py','companion.json','transport.rs','finite-result.json','endpoint.rs','endpoint-result.json'],
'normalization':['normalization-64/normalization.json']}
dirs={'producer':'integrated-moment-expansion','continuation':'integrated-moment-continuation','candidate':'integrated-flow-physical-fit','conditional':'integrated-flow-conditional-transport','normalization':'integrated-moment-series-comparison'}
bindings=[]
for group,names in sources.items():
    for name in names:
        path=BASE/('2026-10-10-'+dirs[group])/name
        assert path.is_file(),path
        bindings.append({'group':group,'path':str(path.relative_to(ROOT)),'sha256':sha(path),'bytes':path.stat().st_size})
put('source-bindings.json',{'schema':1,'status':'exact-read-only-bindings','files':bindings,'no_reference_values_used_for_bound_derivation':True,'no_raised_holdout_parsed':True})
files=[{'path':str(f.relative_to(HERE)),'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(HERE.rglob('*')) if f.is_file() and f.name!='artifact-manifest.json']
put('artifact-manifest.json',{'schema':1,'status':'frozen','scope':'independent source/sign/normalization/conditional-ODE and finite-circle majorant audit; no physical eta0 acceptance','files':files})
print(json.dumps({'rational_checks':len(checks),'files':len(files),'manifest_sha256':sha(HERE/'artifact-manifest.json'),'raised_norm_sha256':sha(HERE/'raised-norm.md')}))
