from pathlib import Path
from fractions import Fraction
import itertools,json,time,sys
N=6;names='abcdef';start=time.monotonic();LIMIT=40
U=['ab','ad','af','bc','cd','cf','be','de','ef','bf','df'];Phi=['abc','acd','acf','bcf','abd','bcd','bde','bdf','adf']
def hepp(s,order):
 e=tuple(s.count(x)for x in names);return tuple(sum(e[i]for i in order[j:])for j in range(1,6))+(0,)
def change(v,S,p):
 r=list(v);r[p]=sum(v[i]for i in S);return tuple(r)
def factor(terms):
 m=tuple(min(t[j]for t in terms)for j in range(N));r=[tuple(t[j]-m[j]for j in range(N))for t in terms];return m,r,not any(m)and tuple([0]*N)in terms or m in terms
count=0;slopes={};bad=[];zero=[];maxdepth=0;rows=[]
def rec(u,psi,jac,eta,trail,order):
 global count,maxdepth
 if time.monotonic()-start>LIMIT:raise TimeoutError
 if len(trail)>24:raise RuntimeError(('depth',trail))
 f,r,unit=factor(psi)
 if unit:
  um,_,uu=factor(u);assert uu
  count+=1;maxdepth=max(maxdepth,len(trail));data=[]
  for i,m in enumerate(eta):
   if m:
    slope=Fraction(2*f[i]-3*um[i],2*m);constant=Fraction(jac[i]+1+6*(um[i]-f[i]),m)
    slopes[str(slope)]=slopes.get(str(slope),0)+1;d={'axis':i,'slope':str(slope),'constant':str(constant),'m':m}
    data.append(d)
    if slope<0:bad.append({'order':order,'trail':trail,**d})
    if slope==0 and constant<0:zero.append({'order':order,'trail':trail,**d})
  rows.append({'order':order,'trail':trail,'U_monomial':um,'F_eta_monomial':f,'jacobian':jac,'eta_exponents':eta,'pole_exponents':data});return
 supports=[set(i for i,n in enumerate(t)if n)for t in r];options=[]
 for k in range(2,N+1):
  options=[s for s in itertools.combinations(range(N),k)if all(set(s)&a for a in supports)]
  if options:break
 S=min(options,key=lambda s:sum(sum(t[i]for i in s)for t in r))
 for p in S:
  newu=tuple(set(change(t,S,p)for t in u));newpsi=tuple(set(change(t,S,p)for t in psi));newjac=list(change(jac,S,p));newjac[p]+=len(S)-1
  rec(newu,newpsi,tuple(newjac),change(eta,S,p),trail+[{'set':S,'pivot':p}],order)
completed=0;status='complete'
try:
 for order in itertools.permutations(range(6)):
  u=tuple(hepp(s,order)for s in U);phi=[hepp(s,order)for s in Phi];alpha=[hepp(s,order)for s in names];psi=phi+[tuple(a+b+(1 if i==5 else 0)for i,(a,b)in enumerate(zip(v,w)))for v in u for w in alpha]
  rec(u,tuple(set(psi)),(4,3,2,1,0,-1),(0,0,0,0,0,1),[],''.join(names[i]for i in order));completed+=1
except (TimeoutError,RuntimeError)as e:status=repr(e)
report={'scope':'eta Mellin scalar prism chart diagnostic, not a proof until all tensor/jet exponents and UV subtraction remainders checked','status':status,'completed_hepp_charts':completed,'sector_count':count,'max_depth':maxdepth,'seconds':time.monotonic()-start,'D_slopes':slopes,'negative_D_slopes':bad,'zero_D_negative_constants':zero,'sectors':rows};Path('/tmp/prism-joint-eta-sector-probe.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items()if k not in ['sectors','negative_D_slopes','zero_D_negative_constants']});print('negative slope',len(bad),'zero slope negative constants',len(zero));print('first bad',bad[:2],zero[:2])
