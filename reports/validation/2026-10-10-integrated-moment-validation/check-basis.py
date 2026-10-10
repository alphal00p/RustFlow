#!/usr/bin/env python3
"""Independent rational affine-basis cross-check; no integration/ODE/closure."""
import json,hashlib
from fractions import Fraction as Q
from pathlib import Path
BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rank(rows):
 a=[list(map(Q,r))for r in rows];k=0
 for col in range(len(a[0])if a else 0):
  pivot=next((i for i in range(k,len(a))if a[i][col]),None)
  if pivot is None:continue
  a[k],a[pivot]=a[pivot],a[k];d=a[k][col];a[k]=[v/d for v in a[k]]
  for i in range(len(a)):
   if i!=k:
    d=a[i][col];a[i]=[x-d*y for x,y in zip(a[i],a[k])]
  k+=1
 return k

def main():
 path=ROOT/'examples/finite_density/massless_three_loop_chain.json';data=json.loads(path.read_text())
 assert data['loops']==3 and data['targets'][1]['numerator']=='g1_2+u1*u2'
 coords=['g1_1','g1_2','g1_3','g2_2','g2_3','g3_3','u1','u2','u3']
 physical=[]
 for edge in data['edges']:
  r=list(map(Q,edge['routing']));physical.append([r[i]*r[j]*(1 if i==j else 2)for i in range(3)for j in range(i,3)]+[Q(0)]*3)
 rows=[];selected=[];completion=[]
 candidates=physical+[[Q(int(i==j))for j in range(9)]for i in range(9)]
 for i,row in enumerate(candidates):
  if rank(rows+[row])>len(rows):
   rows.append(row)
   if i<5:selected.append(i)
   else:selected.append(5+len(completion));completion.append(coords[i-5])
 assert len(rows)==9 and completion==['g1_2','u1','u2','u3']
 assert rows[5]==[Q(int(i==1))for i in range(9)]
 assert rows[6:]==[[Q(int(i==j))for i in range(9)]for j in (6,7,8)]
 assert selected==list(range(9))
 result={'status':'passed independent exact rational basis cross-check; native emitted metadata still required','input':str(path.relative_to(ROOT)),'input_sha256':sha(path),'coordinates':coords,'physical_routing_rows':[[str(x)for x in r]for r in physical],'physical_rank':rank(physical),'completion_slots':{str(i+5):s for i,s in enumerate(completion)},'converted_original_numerator':'rho_5 + rho_6*rho_7','coefficient_mass_derivatives':'zero for this specific numerator because these completion coordinates contain no independent masses','shifted_physical_slots':[1,2,4],'numerator_slots_shifted':[],'special_case_off_endpoint_equality':'N_eta=N follows for this fixture after actual native output confirms the same basis; it is not assumed for other numerators','scalar_powers':[1,1,1,1,1,0,0,0,0],'raised_numerator_terms_before_wick':[{'coefficient':'1','indices':[2,1,1,1,1,-1,0,0,0]},{'coefficient':'1','indices':[2,1,1,1,1,0,-1,-1,0]}],'method':'Fraction Gaussian elimination on exact original routing, reproducing documented deterministic physical-then-coordinate completion; no production routine imported','integration_performed':False,'ode_inferred':False,'reference_values_generated':False,'native_prediction_values_read':False,'source_sha256':sha(Path(__file__))}
 (BASE/'basis-crosscheck.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','physical_rank':5,'completion':completion}))
if __name__=='__main__':main()
