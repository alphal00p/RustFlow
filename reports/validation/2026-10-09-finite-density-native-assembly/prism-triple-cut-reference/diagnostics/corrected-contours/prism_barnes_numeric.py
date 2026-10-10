"""Binary64 continuation diagnostic only; not a production/reference backend."""
import cmath,math,json,sys,time
from fractions import Fraction
from pathlib import Path
CO=[.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7]
def lgamma(z):
 z=complex(z)
 if z.real<.5:return math.log(math.pi)-cmath.log(cmath.sin(math.pi*z))-lgamma(1-z)
 z-=1;s=CO[0]
 for k,c in enumerate(CO[1:],1):s+=c/(z+k)
 t=z+7.5
 return .5*math.log(2*math.pi)+(z+.5)*cmath.log(t)-t+cmath.log(s)
def f(x):return float(Fraction(x))
def group_kernel(items,y):return cmath.exp(sum(sign*lgamma(c+1j*s*y)for c,s,sign in items))
def leaf_value(leaf,step,trunc,t=1.):
 axes=sorted(int(x)for x in leaf['contour']);cont={int(k):f(v)for k,v in leaf['contour'].items()};groups={};constant=0j
 for name,sign in [('num',1),('den',-1)]:
  for encoded in leaf[name]:
   g=list(map(f,encoded));c=g[0]+g[1]*t+sum(g[i]*cont[i]for i in axes);ks=tuple(g[i]for i in axes)
   if all(x==0 for x in ks):
    if c<=0 and c==int(c):
     if sign==-1:return 0j
     raise ValueError(('uncancelled numerator pole',leaf,g))
    constant+=sign*lgamma(c);continue
   nonzero=next(k for k in ks if k);key=tuple(k/nonzero for k in ks)
   groups.setdefault(key,[]).append((c,nonzero,sign))
 fac=f(leaf['coefficient'])*cmath.exp(constant)
 n=math.ceil(trunc/step);ys=[step*i for i in range(-n,n+1)]
 if not axes:return fac
 if len(axes)==1:
  items=groups.get((1.,),[])
  vals=[group_kernel(items,y) for y in ys]
  return fac*complex(math.fsum(v.real for v in vals),math.fsum(v.imag for v in vals))*step/(2*math.pi)
 assert set(groups)<=set([(1.,0.),(0.,1.),(1.,1.)]),groups
 a=[group_kernel(groups.get((1.,0.),[]),y)for y in ys]
 b=[group_kernel(groups.get((0.,1.),[]),y)for y in ys]
 c=[group_kernel(groups.get((1.,1.),[]),i*step)for i in range(-2*n,2*n+1)]
 vals=[]
 for i,av in enumerate(a):
  row=[av*bv*c[i+j]for j,bv in enumerate(b)]
  vals.append(complex(math.fsum(x.real for x in row),math.fsum(x.imag for x in row)))
 return fac*complex(math.fsum(v.real for v in vals),math.fsum(v.imag for v in vals))*(step/(2*math.pi))**2

def outside(row):
 v=lambda g:f(g[0])+f(g[1])
 log=sum(lgamma(v(g))for g in row['outside_num'])-sum(lgamma(v(g))for g in row['outside_den'])+math.log(4)*v(row['four_exponent'])
 return cmath.exp(log)
def test():
 # Γ(a)2^-a = MB Γ(-z)Γ(a+z); initiala1, contour-.4, finala=-1/3.
 leaf={'num':[['0','0','-1','0'],['1','-4/3','1','0']],'den':[],'contour':{'2':'-2/5'},'coefficient':'1'}
 exact=cmath.exp(lgamma(-1/3))*2**(1/3);residue=cmath.exp(lgamma(-1/3));value=leaf_value(leaf,.02,20)+residue
 assert abs(value/exact-1)<1e-12,(value,exact)
 print('Barnes crossed-pole sign toy',value,exact,flush=True)
if __name__=='__main__':
 test();data=json.loads(Path(sys.argv[1]).read_text());step=float(sys.argv[2]);trunc=float(sys.argv[3]);start=time.monotonic();rows=[]
 for i,row in enumerate(data['rows']):
  vals=[leaf_value(l,step,trunc)for l in row['leaves']];raw=complex(math.fsum(v.real for v in vals),math.fsum(v.imag for v in vals));v=outside(row)*raw
  rows.append(v);print(i,v,'dim',len(row['leaves']),flush=True)
 total=complex(math.fsum(v.real for v in rows),math.fsum(v.imag for v in rows));D=f(data['rows'][0]['D'][0])+f(data['rows'][0]['D'][1]);d=D-1;Ad=2**(1-d)*math.pi**(-d/2)/math.gamma(d/2);total*=Ad**3/(8*(4*math.pi)**(D/2))
 report={'status':'binary64 regulated Barnes diagnostic only, not a reference','epsilon':data['epsilon'],'regulator':data.get('analytic_regulator'),'step':step,'truncation':trunc,'value':[total.real,total.imag],'wall_seconds':time.monotonic()-start,'rows':[[x.real,x.imag]for x in rows]};Path(sys.argv[4]).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
