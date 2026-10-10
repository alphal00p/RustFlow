from prism_barnes_continuation import *
from prism_barnes_numeric import leaf_value,lgamma
import cmath,math
m0=F(5,2);m1=F(-7,4);num=[v(z=-1),v(w=-1),v(m0,m1-m0,1,1)]
leaves=[];problems=[];continue_branch(num,[],F(0),{2:F(-7,10),3:F(-9,10)},F(1),[],leaves,problems)
print('leaves',len(leaves),'problems',problems)
encoded=[{'num':[encode(g)for g in l['num']],'den':[encode(g)for g in l['den']],'contour':{str(k):str(x)for k,x in l['contour'].items()},'coefficient':str(l['coefficient']),'history':l['history']}for l in leaves]
expected=cmath.exp(lgamma(float(m1)))*3**(-float(m1))
for step in[.1,.05,.025]:
 vals=[leaf_value(l,step,16)for l in encoded];value=sum(vals);print(step,value,expected,abs(value/expected-1))
