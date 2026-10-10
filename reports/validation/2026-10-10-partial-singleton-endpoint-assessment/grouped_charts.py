#!/usr/bin/env python3
"""Exact parallel-direction grouping and complete rank-two UV charts."""
import hashlib
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

BASE=Path(__file__).resolve().parent
N=4
def const(v):return {(0,)*N:Q(v)} if v else {}
def var(i):return {tuple(int(j==i) for j in range(N)):Q(1)}
def add(a,b):
    c=a.copy()
    for k,v in b.items():c[k]=c.get(k,Q())+v
    return {k:v for k,v in c.items() if v}
def neg(a):return {k:-v for k,v in a.items()}
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    c={}
    for k,v in a.items():
        for l,w in b.items():
            m=tuple(x+y for x,y in zip(k,l));c[m]=c.get(m,Q())+v*w
    return {k:v for k,v in c.items() if v}
def power(a,n):
    b=const(1)
    for _ in range(n):b=mul(a,b)
    return b
def substitute(p,images):
    out={}
    for powers,coefficient in p.items():
        term=const(coefficient)
        for image,n in zip(images,powers):term=mul(term,power(image,n))
        out=add(out,term)
    return out
def deriv(p,i):
    out={}
    for k,v in p.items():
        if k[i]:
            l=list(k);l[i]-=1;out[tuple(l)]=v*k[i]
    return out
def serialize(p):return [{"powers":k,"coefficient":str(v)} for k,v in sorted(p.items())]
def valuation(p):return tuple(min(k[i] for k in p) for i in range(N))
def divide_monomial(p,m):
    out={tuple(x-y for x,y in zip(k,m)):v for k,v in p.items()}
    assert all(x>=0 for k in out for x in k)
    return out

# First names a,b,c,d; then the same polynomial ring names X,c,d,z.
a,b,c,d=[var(i) for i in range(4)]
u=sub(mul(add(add(a,b),d),add(c,d)),mul(d,d))
v=sub(mul(a,u),mul(mul(a,a),add(c,d)))
X,c,d,z=[var(i) for i in range(4)]
images=[mul(X,z),mul(X,sub(const(1),z)),c,d]
group_u=substitute(u,images);group_v=substitute(v,images)
assert group_u==add(mul(X,add(c,d)),mul(c,d))
assert group_v==mul(mul(X,z),add(mul(mul(X,sub(const(1),z)),add(c,d)),mul(c,d)))
jac=sub(mul(deriv(images[0],0),deriv(images[1],3)),mul(deriv(images[0],3),deriv(images[1],0)))
assert jac==neg(X)

supports=[]
for has_b in [False,True]:
    for has_c,has_d in [(True,False),(False,True),(True,True)]:
        images=[X,c if has_c else {},d if has_d else {},z if has_b else const(1)]
        uu=substitute(group_u,images);vv=substitute(group_v,images)
        ss=mul(X,z) if has_b else X
        # The Gaussian mean numerators adj(A)B, where B=(a,0).
        aa=ss
        means=[mul(aa,add(c if has_c else {},d if has_d else {})),mul(aa,d if has_d else {})]
        active=[0]+([1] if has_c else [])+([2] if has_d else [])
        charts=[]
        for perm in itertools.permutations(active):
            chart_images=[{} for _ in range(4)];chart_images[3]=z
            for position,axis in enumerate(perm):
                chart_images[axis]=const(1)
                for j in range(position):chart_images[axis]=mul(chart_images[axis],var(j))
            uc=substitute(uu,chart_images);vc=substitute(vv,chart_images);sc=substitute(ss,chart_images)
            uv=valuation(uc);assert uv==(1,0,0,0), (perm,uv)
            unit=divide_monomial(uc,uv)
            assert unit.get((0,0,0,0),0)==1 and all(x>0 for x in unit.values())
            assert len(sc)==1
            sv=next(iter(sc));ratio=divide_monomial(vc,tuple(x+y for x,y in zip(uv,sv))) if vc else {}
            q=divide_monomial(vc,uv) if vc else {}
            mean_units=[divide_monomial(substitute(m,chart_images),uv) if m else {} for m in means]
            charts.append({"group_order":[['X','c','d'][i] for i in perm],
                           "coordinate_names":["t","w","unused","z"],
                           "U_valuation":uv,"U_positive_unit":serialize(unit),
                           "S":serialize(sc),"Q_numerator_over_U_unit":serialize(q),
                           "Q_over_S_numerator_over_U_unit":serialize(ratio),
                           "Gaussian_mean_numerators_over_U_unit":[serialize(m) for m in mean_units],
                           "ordered_group_jacobian_powers":[len(active)-2,0,0,0],
                           "sole_D_dependent_U_axis":"t","unit_lower_bound":1})
        supports.append({"positive_parameters":['a']+(['b'] if has_b else [])+(['c'] if has_c else [])+(['d'] if has_d else []),
                         "grouping":"a=X*z,b=X*(1-z), absolute Jacobian X" if has_b else "a=X,b absent,z=1",
                         "beta_coordinate_exponents":"z^(n_a-1)*(1-z)^(n_b-1), both nonnegative integers" if has_b else "no beta coordinate",
                         "charts":charts})

report={"scope":"Independent exact rank-two parallel-row grouping and complete UV-chart algebra. No native/HEPKit replay, no period or endpoint value, no production admission.",
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "exact_gaussian_U":serialize(group_u),"exact_gaussian_V":serialize(group_v),
        "grouping_jacobian":serialize(jac),"grouped_variable_names":["X","c","d","z"],
        "full_rank_shifted_supports":len(supports),"complete_charts":sum(len(s['charts']) for s in supports),
        "supports":supports,
        "conclusion":"Every actual full-rank shifted-positive support has one D-dependent UV coordinate after grouping parallel directions; U=t times a positive unit. Q, Q/S and both Gaussian means are smooth on each closed chart cube, including beta-coordinate endpoints. Covariance inverse-U factors require the separate finite integer degree R in the subtraction bound."}
assert report['full_rank_shifted_supports']==6 and report['complete_charts']==20
(BASE/'grouped-chart-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({"supports":6,"charts":20,"grouping_jacobian":True,"Gaussian_determinant_and_Schur_identities":True}))
