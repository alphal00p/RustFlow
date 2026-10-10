"""Exact finite-bound and retained native-polynomial checks; not a limit proof."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import hashlib,json
BASE=Path(__file__).resolve().parent
ORIGIN=BASE.parent/'2026-10-10-partial-origin-capability-prototype'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def evaluate(poly,alpha):
    return sum((Q(t['coefficient'])*__import__('functools').reduce(lambda a,b:a*b,(a**e for a,e in zip(alpha,t['exponents'])),Q(1)) for t in poly),Q(0))
sigma=Q(1,2); witnesses=[]
for h,k,n,ell,s,p,r,p0 in product([1,2],[1,2,3],[1,3],[0,4],[0,3],[1,6],[0,5],[0,3]):
    ji=n-1+max(s-1,0)+max(ell-1,0);j=k*ji;power=Q(p0+j)+sigma
    bounds=[Q(2*(p+r+j+3)),2*(power+1),2+2*ji+power+sigma,2*n+ell+2*max(s-1,0)+power+sigma,Q(k+2)]
    bound=max(bounds);integer=(bound.numerator+bound.denominator-1)//bound.denominator;d=Q(integer)+Q(1,5)
    count=0 if h==1 else (d/2+r).__floor__()+1
    margins={'virtual':h*d/2-p-j-count-sigma,'angular':d/2-1-power,'radial':d-2-2*ji-power-sigma,'lower':d-2*n-ell-2*max(s-1,0)-power-sigma,'positive_angle_measure':d-k-2}
    assert all(v>0 for v in margins.values())
    if h==2:assert d/2+r<count<d/2+r+2
    witnesses.append({'h':h,'compact':k,'n':n,'lower':ell,'upper':s,'P':p,'R':r,'P0':p0,'J':j,'D':str(d),'N_terms':count,'positive_margins':{x:str(y)for x,y in margins.items()}})
records=[]
for path in sorted((ORIGIN/'final-evidence').glob('cuts-*.json')):
    cap=json.loads(path.read_text());checks=zero_faces=0
    for support in cap['complete_supports']:
        full=support['full']
        if full is None:continue
        for alpha in product([Q(0),Q(1)],repeat=len(support['active'])):
            u=evaluate(full['u'],alpha)
            if u<=0:continue
            # At physical null shells and eta=epsilon=0 the retained native
            # polynomial reduces to the exact nonnegative pair sum.
            for angle in [Q(0),Q(1,2),Q(1)]:
                real=Q(0)
                for pair,poly in full['pairs']:
                    c=evaluate(poly,alpha);assert c>=0
                    ri=Q(pair[0]+1);rj=Q(pair[1]+1)
                    real+=c*4*ri*rj*angle
                assert real>=0
                zero_faces+=int(real==0);checks+=1
    records.append({'source':str(path.relative_to(BASE.parent)),'sha256':sha(path),'eta_zero_face_angle_checks':checks,'vanishing_F0_checks':zero_faces})
out={'scope':'exact arithmetic witnesses and replay of already native-audited coefficient records; not fresh HEPKit, limit, or numerical integration','sigma':str(sigma),'finite_witnesses':len(witnesses),'witnesses':witnesses,'native_records':records,'strict_gap_all_pass':True,'cone_zero_counterexample':{'exponent':'5/4','jet_order':2,'function_limit_at_zero':0,'second_derivative_exponent':'-3/4','meaning':'pointwise cone zero alone cannot justify mass jets'}}
(BASE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'finite_witnesses':len(witnesses),'eta_zero_checks':sum(x['eta_zero_face_angle_checks']for x in records),'vanishing_F0_checks':sum(x['vanishing_F0_checks']for x in records),'output_sha256':sha(BASE/'checks.json')}))
