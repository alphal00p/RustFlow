"""Independent binary64 Cartesian-angle quadrature, validation only."""
import math,json,time
from pathlib import Path

def nodes(n):
 out=[]
 for i in range((n+1)//2):
  x=math.cos(math.pi*(i+.75)/(n+.5))
  for _ in range(30):
   a,b=1.,x
   for k in range(2,n+1):a,b=b,((2*k-1)*x*b-(k-1)*a)/k
   der=n*(x*b-a)/(x*x-1);delta=b/der;x-=delta
   if abs(delta)<2e-16:break
  w=2/((1-x*x)*der*der);out.append((x,w))
  if i!=n-1-i:out.append((-x,w))
 return sorted(out)

def direct(p,a,b,c,n):
 # x=n1.n2,y=n1.n3,z=relative transverse cosine.
 # x,y densities (1-x²)^((p-2)/2); z density (1-z²)^((p-3)/2).
 q=nodes(n);xy=[(x,w*(1-x*x)**((p-2)/2),math.sqrt(1-x*x))for x,w in q]
 zz=[(z,w*(1-z*z)**((p-3)/2))for z,w in q]
 normxy=math.gamma((p+1)/2)/(math.sqrt(math.pi)*math.gamma(p/2))
 normz=math.gamma(p/2)/(math.sqrt(math.pi)*math.gamma((p-1)/2))
 rows=[]
 for x,wx,sx in xy:
  for y,wy,sy in xy:
   v=math.fsum(wz*((1-x*y-sx*sy*z)/2)**c for z,wz in zz)
   rows.append(wx*wy*((1-x)/2)**a*((1-y)/2)**b*v)
 return math.fsum(rows)*normxy**2*normz

def exact(p,a,b,c):
 al=p/2
 return math.exp(2*math.lgamma(p)+sum(math.lgamma(al+s)for s in[a,b,c])+math.lgamma(p+a+b+c)-3*math.lgamma(al)-sum(math.lgamma(p+s)for s in[a+b,a+c,b+c]))
start=time.monotonic();results=[]
for p,a,b,c in[(9,1/3,2/3,4/3),(10,-1/3,2/3,-2/3),(8,1,1,1)]:
 expected=exact(p,a,b,c);ref=[]
 for n in[24,40,64,96]:
  value=direct(p,a,b,c,n);ref.append({'order':n,'value':value,'relative_gamma_difference':abs(value/expected-1)})
  print(p,a,b,c,n,ref[-1],flush=True)
 results.append({'sphere_dimension':p,'powers':[a,b,c],'gamma_formula':expected,'quadratures':ref})
Path('/tmp/prism-angular-quadrature.json').write_text(json.dumps({'status':'binary64 independent diagnostic; not a Laurent reference','method':'direct three relative cosines with normalized Gegenbauer weights','results':results,'wall_seconds':time.monotonic()-start},indent=2)+'\n')
