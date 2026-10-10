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
 # Translate one coordinate at a time. Residues of that translation no
 # longer depend on its interpolation parameter; this is ordinary Cauchy.
 result=[{'num':num,'den':den,'contour':cont,'coefficient':F(leaf['coefficient']),'history':[]}];problems=[]
 for axis in axes:
  next_result=[]
  for item in result:
   if axis not in item['contour']:next_result.append(item);continue
   oldcont=item['contour'];delta=new[axis]-oldcont[axis]
   def translate(g):return (g[0]+g[1],g[axis]*delta,g[2],g[3])
   pieces=[]
   continue_branch([translate(g)for g in item['num']],[translate(g)for g in item['den']],F(0),oldcont,item['coefficient'],item['history'],pieces,problems,forced_axis=axis)
   for piece in pieces:
    remains=axis in piece['contour']
    def finish(g):return (g[0]+g[1]-(g[axis]*delta if remains else 0),F(0),g[2],g[3])
    piece['num']=[finish(g)for g in piece['num']];piece['den']=[finish(g)for g in piece['den']]
    if remains:piece['contour']={**piece['contour'],axis:new[axis]}
    next_result.append(piece)
  result=next_result
 encoded=[]
 for item in result:
  encoded.append({'num':[encode(g)for g in item['num']],'den':[encode(g)for g in item['den']],'contour':{str(k):str(x)for k,x in item['contour'].items()},'coefficient':str(item['coefficient']),'history':leaf['history'],'contour_translation':item['history']})
 if len(axes)==2:
  refined=[]
  for item in encoded:
   if len(item['contour'])==1:
    new_items,more_problems=optimize(item);refined.extend(new_items);problems.extend(more_problems)
   else:refined.append(item)
  encoded=refined
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
