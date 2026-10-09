from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json,sys
sys.path.insert(0,'/tmp')
from rustflow_free_virtual_flats import rank,inverse,mul

def add(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,Q(0))+c
    return {m:c for m,c in out.items() if c}
def times(p,q):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            m=tuple(i+j for i,j in zip(a,b));out[m]=out.get(m,Q(0))+x*y
    return {m:c for m,c in out.items() if c}
def scale(p,c):return {m:x*c for m,x in p.items() if x*c}
def rows(cuts):
    data=json.loads(Path('examples/finite_density/triangular_prism.json').read_text())
    original=[[Q(x) for x in e['routing']] for e in data['edges']];L=data['loops']
    occupied=[[x*(1 if data['edges'][i]['charges'][0]>0 else -1) for x in original[i]] for i in cuts]
    basis=occupied[:]
    for j in range(L):
        unit=[Q(i==j) for i in range(L)]
        if rank(basis+[unit],L)>len(basis):basis.append(unit)
    return mul(original,inverse(basis))
def compute(cuts):
    routed=rows(cuts)
    slots=[s for s,r in enumerate(routed) if s not in cuts and any(r[2:])]
    external=[(s,r[:2]) for s,r in enumerate(routed) if s not in cuts and not any(r[2:])]
    assert all(b[0]+b[1]==0 and b[0] for _,b in external),external
    n=len(slots);assert n==6
    x=[{tuple(int(i==j) for i in range(n)):Q(1)} for j in range(n)]
    A=[[add(*(scale(x[j],routed[s][2+i]*routed[s][2+k]) for j,s in enumerate(slots))) for k in range(2)] for i in range(2)]
    B=[[add(*(scale(x[j],routed[s][2+i]*routed[s][k]) for j,s in enumerate(slots))) for k in range(2)] for i in range(2)]
    C12=add(*(scale(x[j],routed[s][0]*routed[s][1]) for j,s in enumerate(slots)))
    U=add(times(A[0][0],A[1][1]),scale(times(A[0][1],A[1][0]),-1))
    cross=add(times(A[1][1],times(B[0][0],B[0][1])),scale(times(A[0][1],add(times(B[0][0],B[1][1]),times(B[1][0],B[0][1]))),-1),times(A[0][0],times(B[1][0],B[1][1])))
    V=add(cross,scale(times(U,C12),-1))
    assert U and V and all(c>0 for c in U.values()) and all(c>0 for c in V.values())
    UX=times(U,add(*x)); points=sorted(set([(*m,0) for m in U]+[(*m,0) for m in V]+[(*m,1) for m in UX]))
    assert rank([[Q(a-b) for a,b in zip(v,points[0])] for v in points[1:]],n+1)==n+1
    candidates={}
    # Discovery only, not claimed exhaustive. Integer normals in this bounded
    # grid are followed by exact supporting-face rank verification.
    for r in product(range(-2,3),repeat=n):
        values=[sum(a*b for a,b in zip(r,pt[:n]))+pt[n] for pt in points]
        m=min(values); face=[pt for pt,val in zip(points,values) if val==m]
        if not any(pt[n] for pt in face) or len(face)<n+1:continue
        if rank([[Q(a-b) for a,b in zip(v,face[0])] for v in face[1:]],n+1)!=n:continue
        candidates[tuple(r)]={'parameter_scalings':list(r),'leading_G_power':m,'twice_D_coefficient':-m,'face_points':len(face)}
    def encode(p):return [{'powers':list(m),'coefficient':str(c)} for m,c in sorted(p.items())]
    return {'cuts':cuts,'parameter_slots':slots,'pure_transfer_slots':[s for s,_ in external],'routed_rows':[[str(c) for c in row] for row in routed], 'U':encode(U),'V':encode(V),'positive_coefficients_verified':True,'G_points':len(points),'candidate_eta_facets':list(candidates.values()),'facet_search_complete':False,'admission':False}
if __name__=='__main__':
    out=[compute(c) for c in [[1,5],[1,7],[5,7]]]
    Path('/tmp/rustflow_prism_parametric.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out:print(json.dumps({k:r[k] for k in ['cuts','G_points','candidate_eta_facets','facet_search_complete']}))
