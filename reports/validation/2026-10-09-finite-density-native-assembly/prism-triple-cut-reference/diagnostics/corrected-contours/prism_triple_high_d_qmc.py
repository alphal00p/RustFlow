"""Independent direct original three-cut high-D check, binary64 diagnostics."""
import math,json,time,sys
from pathlib import Path
from fractions import Fraction as Q
D=float(Q(sys.argv[1]));rho=float(Q(sys.argv[2]));count=int(sys.argv[3]);beta=[1+rho/x for x in[7,11,13]];delta=[rho/x for x in[17,19,23]];kappa=[rho/x for x in[29,31,37]];nu=sum(beta)-D/2;p=D-2
A=2**(2-D)*math.pi**(-(D-1)/2)/math.gamma((D-1)/2)
common=A**3/(8*(4*math.pi)**(D/2))*math.gamma(nu)/math.prod(math.gamma(x)for x in beta)
normalxy=math.gamma((p+1)/2)/(math.sqrt(math.pi)*math.gamma(p/2));normalz=math.gamma(p/2)/(math.sqrt(math.pi)*math.gamma((p-1)/2))
def radical(n,b):
 out=0.;fac=1./b
 while n:
  n,d=divmod(n,b);out+=fac*d;fac/=b
 return out

def integrand(us,surface):
 r1=1. if surface else us[0]**(1/p);r2=us[1]**(1/p);r3=us[2]**(1/p)
 x=2*us[3]-1;y=2*us[4]-1;z=2*us[5]-1
 cos23=x*y+math.sqrt((1-x*x)*(1-y*y))*z
 h=2*r1*r2*(1-x);b=2*r1*r3*(1-y);c=2*r2*r3*(1-cos23)
 x0=us[6];x1=(1-x0)*us[7];x2=(1-x0)*(1-us[7]);F=x0*x1*h+x0*x2*b+x1*x2*c
 N=b*b/4+x1*h*h/4+x1*h*b/4-x1*h*c/4+x2*h*b/2
 if surface:core=-.5*N
 else:
  hp=-1+r2/r1;bp=-1+r3/r1
  Na=b*(bp-1)/2+(hp-1)*(x1*(h+b-c)/4+x2*b/2)+h*(x1*(hp+bp)/4+x2*bp/2)
  Fp=x0*x1*hp+x0*x2*bp
  core=Na-N*(hp/h+bp/b+1/(2*r1*r1)+nu*Fp/F)
 weight=(1-x*x)**((p-2)/2)*(1-y*y)**((p-2)/2)*(1-z*z)**((p-3)/2)*8*normalxy**2*normalz
 weight*=math.prod(xx**(bb-1)for xx,bb in zip([x0,x1,x2],beta))*(1-x0)
 weight*=h**delta[0]*b**delta[1]*c**delta[2]*r1**kappa[0]*r2**kappa[1]*r3**kappa[2]
 weight/=p**(2 if surface else 3)
 return weight*core/(h*b*c)*F**(-nu)
start=time.monotonic();bulk=0.;surface=0.;results=[];primes=[2,3,5,7,11,13,17,19]
for i in range(1,count+1):
 us=[radical(i,b)for b in primes];bulk+=integrand(us,False);surface+=integrand(us,True)
 if i in[count//4,count//2,count]:
  results.append({'points':i,'bulk':bulk/i*common,'surface':surface/i*common,'total':(bulk+surface)/i*common});print(results[-1],flush=True)
Path(sys.argv[4]).write_text(json.dumps({'status':'independent direct high-D QMC diagnostic, not an accuracy-certified reference','D':str(Q(sys.argv[1])),'regulator':str(Q(sys.argv[2])),'results':results,'wall_seconds':time.monotonic()-start},indent=2)+'\n')
