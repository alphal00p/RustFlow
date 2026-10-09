import json,itertools
from pathlib import Path
from fractions import Fraction
base=json.loads(Path('/tmp/prism-virtual-sector-probe.json').read_text())
N=5
def transform(v,S,p):
 out=list(v);out[p]=sum(v[i]for i in S);return out
def prepare(poly):
 m=[min(t[j]for t in poly)for j in range(N)];r=[[t[j]-m[j]for j in range(N)]for t in poly];return m,r,any(not any(t)for t in r)
def recurse(polys,jac,matrix,trail,out):
 assert len(trail)<15,(polys,trail)
 f,r,unit=prepare(polys['Phi'])
 if unit:
  fields={k:{'monomial':prepare(v)[0],'residual':prepare(v)[1]}for k,v in polys.items()}
  slopes=[Fraction(2*f[j]-3*fields['U']['monomial'][j],2)for j in range(N)]
  out.append({'fields':fields,'jacobian':jac,'matrix':matrix,'trail':trail,'D_slopes':[str(s)for s in slopes]});return
 supports=[set(i for i,v in enumerate(t)if v)for t in r]
 options=[S for size in range(2,N+1)for S in itertools.combinations(range(N),size)if all(set(S)&t for t in supports)]
 minsize=len(options[0]);options=[s for s in options if len(s)==minsize]
 S=min(options,key=lambda s:sum(sum(t[i]for i in s)for t in r))
 for p in S:
  new={k:[transform(t,S,p)for t in v]for k,v in polys.items()};newjac=transform(jac,S,p);newjac[p]+=len(S)-1
  newmatrix=[transform(row,S,p)for row in matrix]
  recurse(new,newjac,newmatrix,trail+[{'set':S,'pivot':p}],out)
allrows=[]
for chart in base['charts_data']:
 polys={name:[[a+b for a,b in zip(chart[name]['monomial'],t)]for t in chart[name]['residual']]for name in ['U','Phi','Psi']}
 out=[];recurse(polys,[4,3,2,1,0],[[int(i==j)for j in range(N)]for i in range(N)],[],out)
 for row in out:row['hepp_order']=chart['order']
 allrows+=out
bad=[];signs={}
for i,row in enumerate(allrows):
 u=row['fields']['U']['monomial'];f=row['fields']['Phi']['monomial']
 for j in range(N):
  if f[j]>u[j]:
   slope=Fraction(row['D_slopes'][j]);signs[str(slope)]=signs.get(str(slope),0)+1
   if slope<=0:bad.append({'sector':i,'axis':j,'u':u[j],'f':f[j],'jac':row['jacobian'][j],'slope':str(slope),'order':row['hepp_order'],'trail':row['trail']})
report={'scope':'exact positive sector charts for ordinary virtual prism Phi/U; not a contour proof or reference evaluation','sector_count':len(allrows),'max_blowup_depth':max(len(x['trail'])for x in allrows),'eta_sensitive_axis_D_slopes':signs,'nonpositive_slope_axes':bad,'sectors':allrows}
Path('/tmp/prism-virtual-resolved-sectors.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items()if k not in ['sectors','nonpositive_slope_axes']});print('nonpositive',len(bad),bad[:3])
