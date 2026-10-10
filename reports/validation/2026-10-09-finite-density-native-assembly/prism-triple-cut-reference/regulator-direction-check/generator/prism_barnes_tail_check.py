from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,sys
j=json.load(open(sys.argv[1]));counts=Counter();bad=[];gaps=[]
for ri,row in enumerate(j['rows']):
 for li,leaf in enumerate(row['leaves']):
  axes=sorted(map(int,leaf['contour']));
  if not axes:continue
  slopes={}
  for name,sign in [('num',1),('den',-1)]:
   for g in leaf[name]:
    cs=tuple(F(g[a])for a in axes)
    if not any(cs):continue
    fac=next(c for c in cs if c);key=tuple(c/fac for c in cs);slopes[key]=slopes.get(key,F(0))+sign*abs(fac)
  if len(axes)==1:gap=slopes.get((F(1),),F(0))
  else:
   assert set(slopes)<=set([(F(1),F(0)),(F(0),F(1)),(F(1),F(1))])
   a=slopes.get((F(1),F(0)),F(0));b=slopes.get((F(0),F(1)),F(0));c=slopes.get((F(1),F(1)),F(0));gap=min(a+c,b+c,(a+b)/2)
  counts[(len(axes),str(gap))]+=1
  if gap<=0:bad.append([ri,li,str(gap)])
  gaps.append(gap)
print('counts',counts,'bad',len(bad));assert not bad
Path(sys.argv[2]).write_text(json.dumps({'status':'strict uniform Stirling exponential damping checked','epsilon':j['epsilon'],'criterion':'sum sign*abs(linear imaginary Gamma arguments) >= gap*L1(imaginary coordinates); exp(-pi*gap*L1/2) dominates any bounded-rho power','minimum_gap':str(min(gaps)),'counts':[{'dimension':d,'gap':g,'leaves':n}for(d,g),n in sorted(counts.items())],'bad':bad},indent=2)+'\n')
