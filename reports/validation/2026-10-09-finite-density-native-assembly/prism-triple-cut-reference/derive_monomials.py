from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json
# Laurent polynomial variables h,b,c,r1,r2,r3,x0,x1,x2. ν and Γν live in kernel tag.
N=9
def term(coefficient=1,**kw):
 names=['h','b','c','r1','r2','r3','x0','x1','x2'];return {tuple(kw.get(k,0)for k in names):F(coefficient)}
def add(*xs):
 r=defaultdict(F)
 for x in xs:
  for m,c in x.items():r[m]+=c
 return {m:c for m,c in r.items()if c}
def mul(a,b):
 r=defaultdict(F)
 for m,c in a.items():
  for n,d in b.items():r[tuple(x+y for x,y in zip(m,n))]+=c*d
 return {m:c for m,c in r.items()if c}
def scale(a,s):return {m:c*s for m,c in a.items()}
Na=add(term(F(1,4),b=2),term(F(1,4),h=2,x1=1),term(F(1,4),h=1,b=1,x1=1),term(F(-1,4),h=1,c=1,x1=1),term(F(1,2),h=1,b=1,x2=1))
hp=add(term(-1),term(1,r1=-1,r2=1));bp=add(term(-1),term(1,r1=-1,r3=1))
# Complete derivative of N(a), including explicit -a in two places.
Nap=add(mul(term(F(1,2),b=1),add(bp,term(-1))),mul(add(hp,term(-1)),add(term(F(1,4),h=1,x1=1),term(F(1,4),b=1,x1=1),term(F(-1,4),c=1,x1=1),term(F(1,2),b=1,x2=1))),mul(term(1,h=1),add(mul(term(F(1,4),x1=1),add(hp,bp)),mul(term(F(1,2),x2=1),bp))))
base=term(1,h=-1,b=-1,c=-1)
bulk0=mul(base,add(Nap,scale(mul(Na,add(mul(term(1,h=-1),hp),mul(term(1,b=-1),bp),term(F(1,2),r1=-2))),-1)))
bulk1=scale(mul(base,mul(Na,add(mul(term(1,x0=1,x1=1),hp),mul(term(1,x0=1,x2=1),bp)))),-1)
surface=scale(mul(base,Na),F(-1,2)) # r1=mu; overall Ad^3/8 measure r2,r3 and mu^(D−4)
rows=[]
for sector,shift,poly in [('bulk',0,bulk0),('bulk',1,bulk1),('surface',0,surface)]:
 for m,c in sorted(poly.items()):
  A,B,C,*rest=m;k=rest[:3];x=rest[3:]
  assert 2*(A+B+C-shift)+sum(k)==-4 if sector=='bulk' else 2*(A+B+C-shift)+sum(k)==-2
  rows.append({'sector':sector,'nu_shift':shift,'coefficient':str(c),'h_powers':[A,B,C],'energy_powers':k,'simplex_powers':x})
obj={'schema':1,'scope':'definition-only exact original raised prism triple-cut polynomial expansion; no values','dimension':'4-2*epsilon','nu':'3-D/2','common_factor':'A_d^3/(8*(4*pi)^(D/2))','bulk_kernel':'coefficient*Gamma(nu+nu_shift)*F^(-nu-nu_shift)*product h_ij^h_powers*product r_i^energy_powers*product x_i^simplex_powers','surface_measure':'r1=mu; replace its dr1*r1^(D-3) by mu^(D-4), with coefficient including -1/2; remaining two radial measures unchanged','rows':rows}
Path('/tmp/prism-triple-cut-monomials.json').write_text(json.dumps(obj,indent=2)+'\n');print('rows',len(rows),'bulk0',len(bulk0),'bulk1',len(bulk1),'surface',len(surface))
