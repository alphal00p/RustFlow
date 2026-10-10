from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import itertools,json,sys
# Linear form [constant,t,z,w], endpoint t=1. Generic rational epsilon avoids coincidences.
ACTIVE_RHO=F(1,100003)
# c + t*ct + z*cz + w*cw + rho*cr + rho*t*crt
def v(c=0,t=0,z=0,w=0,r=0,rt=0):return (F(c),F(t),F(z),F(w),F(r),F(rt))
def add(*args):return tuple(sum(a[i]for a in args)for i in range(6))
def scale(a,c):return tuple(F(c)*x for x in a)
def constant(a):return not any(a[1:])
def substitute(a,axis,expr):
 b=list(a);factor=b[axis];b[axis]=F(0);return add(tuple(b),scale(expr,factor))
def encode(a):return [str(x)for x in a]
def clean(num,den):
 n=Counter(num);d=Counter(den)
 for x in list(n):q=min(n[x],d[x]);n[x]-=q;d[x]-=q
 return sorted(n.elements()),sorted(d.elements())
def value(a,t,contour):return a[0]+a[1]*t+ACTIVE_RHO*(a[4]+a[5]*t)+sum(a[i]*contour.get(i,F(0))for i in [2,3])
def initial(row,eps,rho=F(1,100003)):
 D=v(7,-3-2*eps);alpha=scale(add(D,v(-2)),F(1,2));p=scale(alpha,2);beta=[v(2,-1,rt=F(1,q))for q in [7,11,13]];nu=add(*beta,scale(D,F(-1,2)),v(row['nu_shift']));delta=[v(3,-3,rt=F(1,q))for q in [17,19,23]];z=v(z=1);w=v(w=1)
 A,B,C=row['h_powers'];x=row['simplex_powers'];k=row['energy_powers'];a=add(v(A),delta[0],z);b=add(v(B),delta[1],w);c=add(v(C),delta[2],scale(nu,-1),scale(z,-1),scale(w,-1));S=add(a,b,c)
 num=[scale(z,-1),scale(w,-1),add(nu,z,w),add(beta[0],v(x[0]),z,w),add(beta[1],v(x[1]),scale(nu,-1),scale(w,-1)),add(beta[2],v(x[2]),scale(nu,-1),scale(z,-1)),add(alpha,a),add(alpha,b),add(alpha,c)]
 den=[add(p,a,b),add(p,a,c),add(p,b,c)]
 radii=[add(p,a,b,v(k[0])),add(p,a,c,v(k[1])),add(p,b,c,v(k[2]))]
 for i,r in enumerate(radii):
  r=add(r,v(rt=F(1,[29,31,37][i])))
  if row['sector']=='surface'and i==0:continue
  num.append(r);den.append(add(r,v(1)))
 # Outside factors do not affect contour residues; evaluate at endpoint later.
 outside_num=[p,p,add(p,S)]
 outside_den=[alpha,alpha,alpha,add(*beta,v(sum(x)),scale(nu,-2)),*beta]
 num,den=clean(num,den)
 candidates=[]
 floatnum=[[float(g[0]+ACTIVE_RHO*g[4]),float(g[1]+ACTIVE_RHO*g[5]),float(g[2]),float(g[3])]for g in num]
 for i,j in itertools.product(range(1,80),repeat=2):
  cz=-i/23;cw=-j/29;margin=min(g[0]+g[2]*cz+g[3]*cw for g in floatnum)
  if margin>0:candidates.append((margin,i,j))
 assert candidates,('no starting contour',row)
 _,i,j=max(candidates);contour={2:F(-i,23),3:F(-j,29)}
 assert all(value(g,F(0),contour)>0 for g in num)
 return num,den,contour,outside_num,outside_den,S,D

def crossings(num,den,start,contour):
 events=[]
 for index,g in enumerate(num):
  if not any(g[i]for i in contour):continue
  a=value(g,start,contour);b=value(g,F(1),contour)
  if a==b:continue
  lo,hi=sorted([a,b])
  for n in range(max(0,int(-hi)-1),max(1,int(-lo)+2)):
   pole=F(-n);t=start+(pole-a)*(1-start)/(b-a)
   if start<t<1:events.append((t,index,n,b<a))
 return sorted(events)

def continue_branch(num,den,start,contour,coef,history,leaves,problems,forced_axis=None):
 num,den=clean(num,den)
 events=crossings(num,den,start,contour)
 # A fixed straight-contour integral is retained, and every crossing adds its residue.
 leaves.append({'num':num,'den':den,'contour':contour,'coefficient':coef,'history':history})
 for t,index,n,decreasing in events:
  g=num[index];axis=forced_axis if forced_axis in contour and g[forced_axis] else next(i for i in sorted(contour)if g[i]);other_num=num[:index]+num[index+1:];expr=scale(add(g,v(n)),F(-1,g[axis]));expr=list(expr);expr[axis]=0;expr=tuple(expr)
  newnum=[substitute(x,axis,expr)for x in other_num];newden=[substitute(x,axis,expr)for x in den];newnum,newden=clean(newnum,newden);newcontour={i:c for i,c in contour.items()if i!=axis}
  singular=[x for x in newnum if constant(x)and x[0].denominator==1 and x[0]<=0]
  if singular:
   problems.append({'kind':'coincident_residue','history':history,'crossing':(str(t),index,n),'gamma':encode(g),'singular':[encode(x)for x in singular]});continue
  from math import factorial
  newcoef=coef*F((-1)**n,factorial(n))*F(1,abs(g[axis]))*(1 if decreasing else -1)
  continue_branch(newnum,newden,t,newcontour,newcoef,history+[{'t':str(t),'gamma':encode(g),'pole':n,'axis':axis}],leaves,problems,forced_axis)

if __name__=='__main__':
 eps=F(sys.argv[1])if len(sys.argv)>1 else F(1,101);rho=F(sys.argv[2])if len(sys.argv)>2 else F(1,100003);ACTIVE_RHO=rho;rows=json.loads(Path('/tmp/prism-triple-cut-monomials.json').read_text())['rows'];report=[];problems=[]
 for i,row in enumerate(rows):
  num,den,c,on,od,S,D=initial(row,eps,rho);leaves=[];continue_branch(num,den,F(0),c,F(row['coefficient']),[],leaves,problems)
  report.append({'row':row,'outside_num':[encode(a)for a in on],'outside_den':[encode(a)for a in od],'four_exponent':encode(S),'D':encode(D),'leaves':[{'num':[encode(g)for g in l['num']],'den':[encode(g)for g in l['den']],'contour':{str(k):str(x)for k,x in l['contour'].items()},'coefficient':str(l['coefficient']),'history':l['history']}for l in leaves]})
 print('rows',len(rows),'leaves',sum(len(r['leaves'])for r in report),'dimensions',dict(Counter(len(l['contour'])for r in report for l in r['leaves'])),'problems',len(problems));print(problems[:2]);Path(sys.argv[3] if len(sys.argv)>3 else '/tmp/prism-barnes-continuation.json').write_text(json.dumps({'epsilon':str(eps),'analytic_regulator':'rho','continuation_seed_rho':str(rho),'status':'experimental contour continuation requiring independent audit and regulator removal','rows':report,'problems':problems},indent=2)+'\n')
