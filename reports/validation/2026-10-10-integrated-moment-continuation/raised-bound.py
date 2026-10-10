"""Source-only raised absolute C2 norm on the unchanged fixed grid."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,importlib.util,json
B=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('bounds',B/'scalar-bound.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
bind=mod.bind;root=mod.upper_fourth_root;display=mod.display
assert not (B/'raised-bound.json').exists()
plan=json.loads((B/'raised-bound-plan.json').read_text())
for k in ['source_derivation','base_plan','transformed_coefficients']:assert bind(plan[k]['path'])==plan[k]
base=json.loads(Path(plan['base_plan']['path']).read_text());data=json.loads(Path(plan['transformed_coefficients']['path']).read_text())
raised=[x for x in data['transformed'] if x['target']==1];assert len(raised)==1;coefficients=list(map(Q,raised[0]['coefficients']));assert len(coefficients)==64 and coefficients[0]==Q(43,264)
# Independently replay all positive radial moment constants and the signed n=0 check.
m=lambda e:Q(1)/(Q(9,2)+e)
bulk=Q(5,2)*(m(0)**2+m(-1)*m(1));surface=(m(0)+3*m(1))/2;derivative=4*(m(2)*m(0)+m(1)**2)
assert bulk==Q(1580,6237) and surface==Q(38,99) and derivative==Q(3808,14157)
a=Q(9,8)*(bulk+surface);b=Q(9,8)*derivative
assert a==Q(plan['constant_B0']) and b==Q(plan['constant_B1'])
assert Q(9,8)*(Q(3,2)*m(0)**2-m(0)/2+m(1))==coefficients[0]
rows=[];best=[];R=Q(base['R'])
for eta in base['eta']:
 w=R/(Q(eta)+R)
 if w==1:rows.append({'eta':eta,'status':'no admissible radius; endpoint not bounded'});continue
 for N in base['terms']:
  S=sum((coefficients[k]*w**k for k in range(N)),Q(0));candidates=[]
  for radius in base['radii']:
   r=Q(radius)
   if not w<r<1:continue
   power1=root(1+r);power5=root((1+r)**5)
   B0=power5/(1-r);B1=r/4*(Q(5,16)*power1/(1-r)+power5/(1-r)**2)
   M=a*B0+b*B1;T=M*(w/r)**N/(1-w/r)
   relative=T/(abs(S)-T) if abs(S)>T else None
   row={'eta':eta,'w':str(w),'terms':N,'radius':radius,'circle_norm_upper':str(M),'exact_partial_sum':str(S),'absolute_tail_upper':str(T),'absolute_tail_decimal':display(T),'relative_error_upper':str(relative) if relative is not None else None,'relative_error_decimal':display(relative) if relative is not None else None,'exact_upward_root_checks':True,'status':'relative bound established' if relative is not None else 'tail does not establish nonzero denominator'}
   rows.append(row)
   if relative is not None:candidates.append(row)
  best.append(min(candidates,key=lambda x:Q(x['relative_error_upper'])) if candidates else {'eta':eta,'terms':N,'status':'no relative bound on fixed radius grid'})
result={'schema':1,'status':'separate raised absolute component tail bound evaluated','plan':bind(B/'raised-bound-plan.json'),'source':bind(Path(__file__)),'root_enclosure_source':bind(B/'scalar-bound.py'),'rows':rows,'best_fixed_grid_bounds':best,'exact_radial_constants':{'bulk':str(bulk),'surface':str(surface),'derivative':str(derivative),'signed_leading_coefficient':str(coefficients[0])},'reference_values_read':False,'pade_values_read':False,'endpoint_bound':False,'physical_normalization_rounding_included':False}
with (B/'raised-bound.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
print(json.dumps({'best_bounds':[{'eta':x['eta'],'terms':x['terms'],'radius':x.get('radius'),'relative_bound':x.get('relative_error_decimal'),'status':x['status']} for x in best]}))
