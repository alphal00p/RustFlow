#!/usr/bin/env python3
"""Independent exact Gaussian/mass-placement algebra, not native admission."""
import hashlib
import importlib.util
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
helper=HERE.parent/'2026-10-10-occupied-coordinate-control/check_routing.py'
spec=importlib.util.spec_from_file_location('polys',helper)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def value(poly, assignment):
    return sum((c*prod(assignment[v] for v in key) for key,c in poly.items()),Q(0))


def prod(values):
    result=Q(1)
    for x in values:result*=x
    return result


def main():
    out=HERE/'checks.json';assert not out.exists()
    u=m.poly('(a+b)*(c+d)+c*d')
    gaussian=m.add(m.mul(m.poly('a+b+d'),m.poly('c+d')),m.poly('d*d'),-1)
    assert u==gaussian
    v=m.add(m.mul(m.poly('a'),u),m.poly('a*a*(c+d)'),-1)
    assert v==m.poly('a*(b*(c+d)+c*d)')
    remainder=m.add(m.mul(m.poly('a'),u),v,-1)
    assert remainder==m.poly('a*a*(c+d)') and all(c>0 for c in remainder.values())
    bubble_u=m.poly('x+y');bubble_v=m.poly('x*y')
    assert m.add(m.mul(m.poly('x'),bubble_u),bubble_v,-1)==m.poly('x*x')
    samples=[]
    for a,b,c,d in itertools.product([Q(0),Q(1),Q(3,2)],repeat=4):
        at=dict(zip('abcd',[a,b,c,d]));den=a*value(u,at)
        if den>0:
            ratio=value(v,at)/den;assert 0<=ratio<=1
            samples.append({'parameters':[str(x) for x in [a,b,c,d]],'ratio':str(ratio)})
    unsupported=[]
    for b in [Q(1,100),Q(1,10000)]:
        at=dict(a=Q(1),b=b,c=Q(1),d=Q(1));ratio=value(v,at)/(b*value(u,at))
        assert ratio>1;unsupported.append({'b':str(b),'V_over_bU':str(ratio)})
    path=ROOT/'examples/finite_density/massless_three_loop_chain.json';original=json.loads(path.read_text())
    rows=[list(map(Q,e['routing'])) for e in original['edges']]
    # P1=q+t, P2=l, P3=t.
    inverse=[[Q(1),Q(1),Q(0)],[Q(0),Q(0),Q(1)],[Q(0),Q(1),Q(0)]]
    routed=[[sum(row[i]*inverse[i][j] for i in range(3)) for j in range(3)] for row in rows]
    assert routed==[[1,1,0],[0,0,1],[0,1,0],[1,0,0],[0,-1,1]]
    assert all(routed[s][0]==0 for s in [1,2,4])
    report={'scope':'Exact independent Gaussian polynomial/routing assessment only; no native HEPKit replay, contour permit, source-zero installation or values',
            'singleton_polynomials':{'U':'(a+b)(c+d)+cd','V':'a[b(c+d)+cd]','partial_mass':'eta*a*U','positive_difference':'aU-V=a²(c+d)'},
            'double_polynomials':{'U':'x+y','V':'xy','partial_mass':'eta*x*U','positive_difference':'xU-V=x²'},
            'singleton_positive_domain_checks':len(samples),'samples':samples,'nonuniform_slot2_examples':unsupported,
            'singleton_routed_physical_rows':[[str(x) for x in row] for row in routed],
            'shift0_leaves_unshifted_factors_independent_of_compact_q':True,
            'inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,helper,Path(__file__)]}}
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'positive_ratio_checks':len(samples),'exact_polynomial_and_routing_checks':'passed'}))


if __name__=='__main__':main()
