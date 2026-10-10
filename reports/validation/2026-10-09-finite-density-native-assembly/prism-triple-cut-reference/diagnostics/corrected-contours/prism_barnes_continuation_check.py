from prism_barnes_continuation import *
from prism_barnes_numeric import leaf_value,lgamma
import cmath,math,json

def barnes(a,b,c,d):return cmath.exp(sum(lgamma(x)for x in[a+c,a+d,b+c,b+d])-lgamma(a+b+c+d))
# u=z+w. The coupled two-dimensional integral factorizes after u substitution.
# Its continuation can independently be checked by twice Barnes's first lemma.
start=[F(1),F(6,5),F(13,10),F(7,5),F(11,10),F(9,10),F(8,5),F(17,10)]
end=[F(-11,10),F(2,5),F(3,5),F(7,5),F(-7,10),F(3,10),F(4,5),F(8,5)]
num=[]
for i,(a,b)in enumerate(zip(start,end)):
 if i<4:sgn=1 if i<2 else -1;num.append(v(a,b-a,z=sgn,w=sgn))
 else:sgn=1 if i<6 else -1;num.append(v(a,b-a,w=sgn))
leaves=[];problems=[];continue_branch(num,[],F(0),{2:F(-1,7),3:F(-1,11)},F(1),[],leaves,problems)
assert not problems,problems
leaves=[{'num':[encode(g)for g in l['num']],'den':[encode(g)for g in l['den']],'contour':{str(k):str(x)for k,x in l['contour'].items()},'coefficient':str(l['coefficient']),'history':l['history']}for l in leaves]
expected=barnes(*map(float,end[:4]))*barnes(*map(float,end[4:]));result=[]
for step in[.025,.0125,.00625]:
 vals=[leaf_value(l,step,12)for l in leaves];value=complex(math.fsum(x.real for x in vals),math.fsum(x.imag for x in vals));error=abs(value/expected-1)
 print(step,value,expected,error,flush=True);result.append({'step':step,'value':[value.real,value.imag],'expected':[expected.real,expected.imag],'relative_error':error})
assert result[-1]['relative_error']<1e-10
Path('/tmp/prism-barnes-continuation-check.json').write_text(json.dumps({'status':'independent exact Barnes-first-lemma continuation diagnostic','leaves':len(leaves),'results':result},indent=2)+'\n')
