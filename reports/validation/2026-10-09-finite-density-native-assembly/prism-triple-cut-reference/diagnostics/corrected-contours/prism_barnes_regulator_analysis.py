from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import json
j=json.load(open('/tmp/prism-barnes-symbolic-contours.json'));rho=F(j['continuation_seed_rho']);counts=Counter();margin=F(1);problems=[];zero=0
for ri,row in enumerate(j['rows']):
 for li,leaf in enumerate(row['leaves']):
  axes=sorted(map(int,leaf['contour']));cont={int(k):F(v)for k,v in leaf['contour'].items()};power=0;identically_zero=False
  for name,sign in[('num',1),('den',-1)]:
   for gs in leaf[name]:
    g=list(map(F,gs));real=g[0]+g[1]+sum(g[a]*cont[a]for a in axes);slope=g[4]+g[5];hasvar=any(g[a]for a in axes)
    pole=real<=0 and real.denominator==1
    if not hasvar:
     if pole:
      if slope==0:
       if sign>0:problems.append(('identical numerator pole',ri,li,gs))
       else:identically_zero=True
      power+=sign
    else:
     real_end=real+slope*rho
     lo,hi=sorted([real,real_end]);n0=(-hi).__ceil__();n1=(-lo).__floor__()
     if n1>=max(0,n0):problems.append(('rho contour crossing',ri,li,gs,str(real),str(real_end)))
     dist=real if real>0 else min(real-real.__floor__(),real.__ceil__()-real);margin=min(margin,dist)
  if identically_zero:zero+=1
  else:counts[(len(axes),power)]+=1
print('degree counts',counts,'identicallyzero',zero,'minmargin',float(margin),'problems',len(problems));print(problems[:2])
Path('/tmp/prism-barnes-symbolic-regulator-analysis.json').write_text(json.dumps({'counts':[{ 'dimension':a,'pole_order':b,'terms':n}for(a,b),n in sorted(counts.items())],'identically_zero_terms':zero,'minimum_nonconstant_gamma_pole_distance_at_rho_zero':str(margin),'problems':problems},indent=2)+'\n')
# Specialize symbolic regulator for the existing numeric reader.
for row in j['rows']:
 for name in['outside_num','outside_den']:
  row[name]=[[str(F(g[0])+rho*F(g[4])),str(F(g[1])+rho*F(g[5])),g[2],g[3]]for g in row[name]]
 for name in['D','four_exponent']:
  g=row[name];row[name]=[str(F(g[0])+rho*F(g[4])),str(F(g[1])+rho*F(g[5])),g[2],g[3]]
 for leaf in row['leaves']:
  for name in['num','den']:
   leaf[name]=[[str(F(g[0])+rho*F(g[4])),str(F(g[1])+rho*F(g[5])),g[2],g[3]]for g in leaf[name]]
j['analytic_regulator']=str(rho);Path('/tmp/prism-barnes-symbolic-specialized.json').write_text(json.dumps(j,indent=2)+'\n')
