from fractions import Fraction as Q
from pathlib import Path
import hashlib,json

def dot(x,y):return x[0]*y[0]-sum(a*b for a,b in zip(x[1:],y[1:]))
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def denom(k,l,p,q):return [dot(x,x) for x in [k,l,sub(k,p),sub(l,p),sub(k,q),sub(k,l)]]
def original(v):
    p1,p2,p3,p4=v
    return dot(p1,p2)**2+dot(p1,p3)*dot(p2,p4)
q=tuple(map(Q,['1/2','1/2','0','0']))
p=tuple(map(Q,['0','1','0','0']))
k=tuple(map(Q,['2','1/3','2/5','-1/7']))
l=tuple(map(Q,['3/2','-2/7','1/5','2/3']))
d0,d1,d2,d3,d4,d5=denom(k,l,p,q)
na=((d0-d4)**2-(d0+d1-d5))/4
nb=(1+(d1-1-d3)*(d0-d4))/4
nc=((d0+d4)**2+(d0-d1+d5)*(d4+d2))/4
assert original([k,q,l,p])==na
assert original([p,q,l,k])==nb
qc=(-q[0],q[1],q[2],q[3])
c0,c1,c2,c3,c4,c5=denom(k,l,p,qc)
nc_physical=((c0+c4)**2+(c0-c1+c5)*(c4+c2))/4
assert original([k,sub(k,qc),sub(k,l),sub(k,p)])==nc_physical
assert sub(p,qc)[0]>0 and -qc[0]>0
v=sub(p,q)
assert 2*dot(q,v)==-1 and 2*dot(sub(q,p),v)==0
assert -2*dot(sub(k,q),v)==d2-d4
nap=(d0-d4)*(d4-d2-1)/2-(d0+d1-d5)/4
nbp=Q(1,2)+(d1-1-d3)*(d4-d2-1)/4
assert 2*dot(k,q)*dot(k,v)+dot(k,l)*dot(p,v)==nap
assert 2*dot(p,q)*dot(p,v)+dot(p,l)*dot(k,v)==nbp

# Exact degree-four formal jets in y=sqrt(t), not finite differences.
class Jet:
    def __init__(self,a):self.a=[Q(a)]+[Q(0)]*4 if not isinstance(a,list) else a
    def __add__(self,b):
        b=b if isinstance(b,Jet) else Jet(b)
        return Jet([x+y for x,y in zip(self.a,b.a)])
    __radd__=__add__
    def __neg__(self):return Jet([-x for x in self.a])
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        b=b if isinstance(b,Jet) else Jet(b)
        return Jet([sum(self.a[j]*b.a[i-j] for j in range(i+1)) for i in range(5)])
    __rmul__=__mul__
    def inverse(self):
        assert self.a[0]
        out=[1/self.a[0]]
        for i in range(1,5):out.append(-sum(self.a[j]*out[i-j] for j in range(1,i+1))/self.a[0])
        return Jet(out)
    def __truediv__(self,b):return self*(b if isinstance(b,Jet) else Jet(b)).inverse()
    def __rtruediv__(self,b):return self.inverse()*b
    def __pow__(self,n):
        out=Jet(1)
        for _ in range(n):out=out*self
        return out
r1,r2,r3=Q(1,3),Q(1,2),Q(2,5)
z,w=Q(4,5),Q(3,5)
x0,x1,x2=Q(1,2),Q(1,3),Q(1,6)
ell=4*r1*r2;b=2*r1*r3*(1-z)
c0=2*r2*r3*(1-z);c1=4*r2*r3*z
sigma=-4*r2*r3*Q(3,5)*w
dh=-1+r2/r1;db=-1+r3/r1
fa=x0*x1*dh+x0*x2*db
y=Jet([Q(0),Q(1),Q(0),Q(0),Q(0)])
t=y*y
checks=[]
for nu in [1,2]:
    def f(sigma):
        h=ell*t;c=Jet(c0)+sigma*y+c1*t-sigma/2*y**3
        ff=x0*x1*h+x0*x2*b+x1*x2*c
        nn=b*b/4+h*x1*(h+b-c)/4+h*x2*b/2
        nna=b*(db-1)/2+(dh-1)*(x1*(h+b-c)/4+x2*b/2)+h*(x1*(dh+db)/4+x2*db/2)
        gg=nn/(b*c*ff**nu)
        hh=(nna-nn*(db/b+nu*fa/ff+1/(2*r1*r1)))/(b*c*ff**nu)
        return -dh/ell**2*gg+t/ell*hh
    averaged=(f(sigma)+f(-sigma))/2
    f0=x0*x2*b+x1*x2*c0;n0=b*b/4;g0=n0/(b*c0*f0**nu)
    nh=x1*(b-c0)/4+x2*b/2
    gh=(nh-nu*n0*x0*x1/f0)/(b*c0*f0**nu)
    gc=-g0*(1/c0+nu*x1*x2/f0)
    gcc=g0*(2/c0**2+2*nu*x1*x2/(c0*f0)+nu*(nu+1)*(x1*x2)**2/f0**2)
    na0=b*(db-1)/2+(dh-1)*(x1*(b-c0)/4+x2*b/2)
    h0=(na0-n0*(db/b+nu*fa/f0+1/(2*r1*r1)))/(b*c0*f0**nu)
    expected0=-dh/ell**2*g0
    expected1=-dh/ell**2*(ell*gh+c1*gc+sigma**2/2*gcc)+h0/ell
    assert averaged.a[0]==expected0 and averaged.a[2]==expected1
    assert averaged.a[1]==0 and averaged.a[3]==0
    checks.append({'nu':nu,'f0':str(expected0),'f1_paired_fixed_w':str(expected1),'odd_sqrt_t_jets':'zero'})
out={'status':'passed','scope':'Exact definition-only numerator/kinematic-derivative and local paired-collinear jet identities; no integral values','two_cut_numerator_maps_checked':3,'off_shell_jet_numerator_maps_checked':2,'collinear_jet_checks':checks,'complete_prism_reference':False,'numerical_oracle_records_compared':0,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path('/tmp/finite_density_prism_subtractions_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
