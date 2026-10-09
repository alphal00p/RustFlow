"""Validation-only exact Laurent audit of saved native matrix coefficients.

No integral elimination, source generation, or feature evaluation occurs here.
"""
import ast,json,pathlib,argparse
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parent
LOW,HIGH=-16,24
class S:
 def __init__(self,d): self.d={k:Q(v) for k,v in d.items() if v and LOW<=k<=HIGH}
 def __add__(self,b):
  z=self.d.copy()
  for k,v in b.d.items():z[k]=z.get(k,Q(0))+v
  return S(z)
 def __neg__(self):return S({k:-v for k,v in self.d.items()})
 def __sub__(self,b):return self+-b
 def __mul__(self,b):
  z={}
  for k,v in self.d.items():
   for l,w in b.d.items():
    if LOW<=k+l<=HIGH:z[k+l]=z.get(k+l,Q(0))+v*w
  return S(z)
 def inv(self):
  k=min(self.d); c=[1/self.d[k]]
  for n in range(1,HIGH+k+1):c.append(-sum((self.d.get(k+j,Q(0))*c[n-j] for j in range(1,n+1)),Q(0))/self.d[k])
  return S({i-k:v for i,v in enumerate(c)})
 def __pow__(self,p):
  if p<0:return self.inv()**(-p)
  r=S({0:1});a=self
  while p:
   if p&1:r=r*a
   a=a*a;p//=2
  return r

def parse(n):
 if isinstance(n,ast.Constant):return S({0:n.value})
 if isinstance(n,ast.Name):return {'x':S({1:1}),'e':S({0:Q(4,5)})}[n.id]
 if isinstance(n,ast.UnaryOp):return -parse(n.operand) if isinstance(n.op,ast.USub) else parse(n.operand)
 a,b=parse(n.left),parse(n.right)
 if isinstance(n.op,ast.Add):return a+b
 if isinstance(n.op,ast.Sub):return a-b
 if isinstance(n.op,ast.Mult):return a*b
 if isinstance(n.op,ast.Div):return a*b.inv()
 if isinstance(n.op,ast.Pow):
  power=b.d.get(0,Q(0));assert set(b.d)<={0} and power.denominator==1
  return a**int(power)
 raise ValueError(ast.dump(n))

arguments=argparse.ArgumentParser()
arguments.add_argument('--input',default='guard-refinement/legacy-3passes/sunset-N9/result.json')
arguments.add_argument('--output',default='n9-endpoint-matrix-audit.json')
options=arguments.parse_args()
d=json.loads((ROOT/options.input).read_text())['detail']
z=[[parse(ast.parse(v.replace('sunset_closure::{}::eta','x').replace('sunset_closure::{}::epsilon','e').replace('^','**'),mode='eval').body) for v in row] for row in d['matrix']]
rows=[]
for i,row in enumerate(z):
 pp={p:{str(j):str(v.d[p]) for j,v in enumerate(row) if p in v.d} for p in range(-8,0)}
 pp={str(p):v for p,v in pp.items() if v}
 if pp:rows.append({'row':i,'basis':d['basis'][i],'principal_part_at_epsilon_4_5':pp})
report={'scope':'Read-only exact rational Laurent expansion of the saved closed matrix, not an elimination backend or feature prediction. Expansion coefficients use Fraction arithmetic.','epsilon':'4/5','rows':rows}
(ROOT/options.output).write_text(json.dumps(report,indent=2)+'\n')
print('pole orders',[(r['row'],list(r['principal_part_at_epsilon_4_5'])) for r in rows])
for i,p in [(23,-2),(25,-4)]:
 if i<len(z):print('row',i,'power',p,{j:str(v.d[p]) for j,v in enumerate(z[i]) if p in v.d})
