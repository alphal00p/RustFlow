"""Exact Cauchy contour translation after a regulated continuation."""
from prism_barnes_continuation import *
import sys

def optimize(leaf):
 num=[tuple(F(x)for x in g)for g in leaf['num']];den=[tuple(F(x)for x in g)for g in leaf['den']];cont={int(k):F(x)for k,x in leaf['contour'].items()}
 if not cont:return [leaf],[]
 # Set physical parameters to their endpoint before translating straight contours.
 num=[v(g[0]+g[1],0,g[2],g[3])for g in num];den=[v(g[0]+g[1],0,g[2],g[3])for g in den]
 axes=sorted(cont);candidates=[]
 floating=[[float(x)for x in g]for g in num if any(g[a]for a in axes)]
 for ij in itertools.product(range(-35,36),repeat=len(axes)):
  new={axis:i/(19 if axis==2 else 23)for axis,i in zip(axes,ij)}
  if any(x==int(x)for x in new.values()):continue
  vals=[g[0]+sum(g[a]*new[a]for a in axes)for g in floating]
  def distance(x):return x if x>0 else min(x-math.floor(x),math.ceil(x)-x)
  margin=min(map(distance,vals));movement=sum(abs(new[a]-float(cont[a]))for a in axes)
  candidates.append((margin,-movement,ij))
 _,_,ij=max(candidates);new={axis:F(i,19 if axis==2 else 23)for axis,i in zip(axes,ij)}
 def translate(g):return (g[0],sum(g[a]*(new[a]-cont[a])for a in axes),g[2],g[3])
 result=[];problems=[];continue_branch([translate(g)for g in num],[translate(g)for g in den],F(0),cont,F(leaf['coefficient']),[],result,problems)
 encoded=[]
 for item in result:
  encoded.append({'num':[encode(g)for g in item['num']],'den':[encode(g)for g in item['den']],'contour':{str(k):str(x)for k,x in item['contour'].items()},'coefficient':str(item['coefficient']),'history':leaf['history'],'contour_translation':item['history']})
 return encoded,problems
if __name__=='__main__':
 import math
 data=json.loads(Path(sys.argv[1]).read_text());problems=[]
 for i,row in enumerate(data['rows']):
  leaves=[]
  for old in row['leaves']:
   new,p=optimize(old);leaves.extend(new);problems.extend(p)
  print(i,len(row['leaves']),len(leaves),flush=True);row['leaves']=leaves
 data['contour_translation_problems']=problems
 Path(sys.argv[2]).write_text(json.dumps(data,indent=2)+'\n');print('problems',len(problems))
