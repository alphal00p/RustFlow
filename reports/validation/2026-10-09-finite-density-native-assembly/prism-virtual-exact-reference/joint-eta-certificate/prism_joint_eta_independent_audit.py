from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json
rows=json.loads(Path('/tmp/prism-joint-eta-sector-probe.json').read_text())['sectors']
names='abcdef';U=['ab','ad','af','bc','cd','cf','be','de','ef','bf','df'];Phi=['abc','acd','acf','bcf','abd','bcd','bde','bdf','adf']
def initial(s,order):return tuple(sum(s.count(names[k])for k in order[j:])for j in range(1,6))+(0,)
def change(v,S,p):return tuple(sum(v[k]for k in S)if i==p else n for i,n in enumerate(v))
def minterm(terms):return tuple(map(min,zip(*terms)))
vol=defaultdict(F); zero=[]; maxbound=F(0); witness=None; paths=defaultdict(set);positive_axes=0
for row in rows:
 order=tuple(names.index(c)for c in row['order']); alpha=[initial(c,order)for c in names];u=[initial(s,order)for s in U];f=[initial(s,order)for s in Phi]+[tuple(x+y+(i==5)for i,(x,y)in enumerate(zip(a,b)))for a in u for b in alpha];j=(4,3,2,1,0,-1);eta=(0,0,0,0,0,1);trail=[]
 for step in row['trail']:
  S=step['set'];p=step['pivot'];paths[(row['order'],tuple(trail))].add((tuple(S),p));trail.append((tuple(S),p))
  alpha=[change(x,S,p)for x in alpha];u=[change(x,S,p)for x in u];f=[change(x,S,p)for x in f];eta=change(eta,S,p);j=tuple(v+(len(S)-1 if i==p else 0)for i,v in enumerate(change(j,S,p)))
 um=minterm(u);fm=minterm(f)
 assert list(um)==row['U_monomial'] and list(fm)==row['F_eta_monomial']
 assert list(j)==row['jacobian'] and list(eta)==row['eta_exponents'] and um in u and fm in f
 volume=F(1)
 for ji,mi in zip(j,eta):assert ji+mi+1>0;volume/=ji+mi+1
 vol[row['order']]+=volume
 zs=[]
 for i,m in enumerate(eta):
  if not m:continue
  slope=F(2*fm[i]-3*um[i],2*m)
  # Rank <=6 and P<=7: worst Gaussian monomial U^(1-3D/2) F^(D-7).
  constant=F(j[i]+1+um[i]-7*fm[i],m)
  assert slope>=0
  if slope==0:
   assert um[i]==fm[i]==0 and j[i]+1==0 and m==1
   assert all(a[i]==0 for a in alpha)
   zs.append(i)
  else:
   positive_axes+=1;bound=-constant/slope
   if bound>maxbound:maxbound=bound;witness={'order':row['order'],'trail':row['trail'],'axis':i,'slope':str(slope),'constant':str(constant)}
 assert len(zs)<=1
 if zs:zero.append((row['order'],zs[0]))
assert len(vol)==720 and all(v==F(1,120)for v in vol.values())
for key,children in paths.items():
 sets={S for S,p in children};assert len(sets)==1
 S=sets.pop();assert {p for _,p in children}==set(S)
result={'scope':'Independent replay of every exact map, positive-unit factor, Jacobian and full branching/volume coverage. No integral values or production admission.','sector_count':len(rows),'chart_count':len(vol),'per_chart_original_cube_volume':'1/120','total_projective_max_gauge_volume':str(sum(vol.values())),'complete_branch_points':len(paths),'zero_D_axes':len(zero),'zero_axis_all_original_alpha_valuations_zero':True,'positive_D_axes':positive_axes,'rank6_P7_strict_D_lower_bound':str(maxbound),'bound_witness':witness}
Path('/tmp/prism-joint-eta-independent-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
