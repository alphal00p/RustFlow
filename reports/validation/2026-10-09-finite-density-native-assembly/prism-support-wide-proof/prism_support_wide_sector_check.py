from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,permutations
import json,time
source=Path('reports/validation/2026-10-09-finite-density-native-assembly/prism-endpoint-proof/rustflow_prism_parametric.json')
raw=json.loads(source.read_text());start=time.monotonic();results=[]
def valuation(poly):
 v=tuple(map(min,zip(*poly)));return v,v in poly
for graph in raw:
 fullu=[tuple(t['powers'])for t in graph['U']];fullv=[tuple(t['powers'])for t in graph['V']];N=len(graph['parameter_slots']);result={'cuts':graph['cuts'],'supports':[]};totalwork=totalleaves=totalnodes=0
 for mask in range(1<<N):
  keep=[i for i in range(N)if mask>>i&1];drop=set(range(N))-set(keep);n=len(keep)
  restrict=lambda poly:set(tuple(t[i]for i in keep)for t in poly if all(t[i]==0 for i in drop))
  u=restrict(fullu);v=restrict(fullv)
  if not u:result['supports'].append({'mask':mask,'rank_deficient':True});continue
  U={t+(0,)for t in u};F={t+(0,)for t in v}
  for t in u:
   for i in range(n):
    e=list(t)+( [1] );e[i]+=1;F.add(tuple(e))
  stats={'mask':mask,'V_zero':not v,'nodes':0,'leaves':0,'work':0,'regular':0,'min_D_slope':None,'worst_D_bound_Pn_rank4_J1':Q(0)}
  def change(t,S,p):return tuple(sum(t[i]for i in S)if j==p else x for j,x in enumerate(t))
  def rec(u,f,j,m,a):
   stats['nodes']+=1;stats['work']+=(len(u)+len(f)+len(a)+2)*n
   uv,uu=valuation(u);fv,ff=valuation(f)
   if uu and ff:
    stats['leaves']+=1;regular=0
    for i,mi in enumerate(m):
     if not mi:continue
     d=2*fv[i]-3*uv[i]
     assert d>=0,(mask,u,f,j,m)
     assert mi+uv[i]>=fv[i]
     if d==0:
      assert uv[i]==fv[i]==j[i]==0 and mi==1 and all(x[i]==0 for x in a)
      regular+=1
     else:
      slope=Q(d,2*mi);stats['min_D_slope']=slope if stats['min_D_slope']is None else min(stats['min_D_slope'],slope)
      c=j[i]+(n-4)*uv[i]-(n+1)*fv[i]
      stats['worst_D_bound_Pn_rank4_J1']=max(stats['worst_D_bound_Pn_rank4_J1'],Q(-2*c,d))
    assert regular<=1;stats['regular']+=regular;return
   poly,minimum=(f,fv)if not ff else(u,uv)
   stats['work']+=(1<<n)*len(poly)*n
   candidates=[]
   for k in range(2,n+1):
    for S in combinations(range(n),k):
     if all(any(t[i]>minimum[i]for i in S)for t in poly):candidates.append((sum(sum(t[i]-minimum[i]for i in S)for t in poly),S))
    if candidates:break
   _,S=min(candidates)
   for p in S:rec({change(t,S,p)for t in u},{change(t,S,p)for t in f},change(j,S,p),change(m,S,p),[change(t,S,p)for t in a])
  for order in permutations(range(n)):
   def hepp(t):return tuple(sum(t[i]for i in order[k:])for k in range(1,n))+(t[n],)
   a=[hepp(tuple(int(i==j)for i in range(n))+(0,))for j in range(n)]
   rec({hepp(t)for t in U},{hepp(t)for t in F},tuple(range(n-1,-1,-1)),(0,)*(n-1)+(1,),a)
  totalwork+=stats['work'];totalnodes+=stats['nodes'];totalleaves+=stats['leaves'];stats['min_D_slope']=str(stats['min_D_slope']);stats['worst_D_bound_Pn_rank4_J1']=str(stats['worst_D_bound_Pn_rank4_J1']);result['supports'].append(stats)
 result.update(total_work=totalwork,total_nodes=totalnodes,total_leaves=totalleaves,rank_deficient=sum(s.get('rank_deficient',False)for s in result['supports']),full_rank_V_zero=sum(s.get('V_zero',False)for s in result['supports']))
 print({k:v for k,v in result.items()if k!='supports'});results.append(result)
Path('/tmp/prism-support-wide-sector-check.json').write_text(json.dumps({'scope':'Independent exact restriction of routed U/V supports; no HEPKit constructor test or flow admission','seconds':time.monotonic()-start,'results':results},indent=2)+'\n')
