"""Exact rational source-based Cauchy bounds; no reference or approximant read."""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,localcontext
from math import isqrt
import hashlib,json
B=Path(__file__).resolve().parent

def bind(p):p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def upper_fourth_root(x):
 assert x>0;scale=10**60;num=x.numerator*scale**4;den=x.denominator
 ceiling=(num+den-1)//den;root=isqrt(isqrt(ceiling))
 if root**4<ceiling:root+=1
 result=Q(root,scale);assert result**4>=x and (result-Q(1,scale))**4<x
 return result

def display(x):
 with localcontext() as ctx:ctx.prec=50;return str(Decimal(x.numerator)/Decimal(x.denominator))

def main():
 assert not (B/'scalar-bound.json').exists();plan=json.loads((B/'bound-plan.json').read_text())
 for item in plan['source_bindings'].values():assert bind(item['path'])==item
 data=json.loads(Path(plan['source_bindings']['gaussian-wide-prefix-64/coefficients.json']['path']).read_text())
 scalar_audits=[x for x in data['audit'] if x['target']==0];assert len(scalar_audits)==1 and scalar_audits[0]['indices'][:5]==[1,1,1,1,1] and all(i==0 for i in scalar_audits[0]['indices'][5:])
 h0=[x for x in data['coefficients'] if x['target']==0 and x['n']==0];assert len(h0)==1 and Q(h0[0]['coefficient'])>0
 R=Q(plan['R']);rows=[];best=[]
 for eta in plan['eta']:
  w=R/(Q(eta)+R)
  if w==1:rows.append({'eta':eta,'status':'no admissible radius; endpoint not bounded'});continue
  for N in plan['terms']:
   candidates=[]
   for rs in plan['radii']:
    r=Q(rs)
    if not w<r<1:continue
    circle_upper=upper_fourth_root((1+r)**5)
    relative_circle_upper=upper_fourth_root(((1+r)/(1-w))**5)
    geometric=(w/r)**N/((1-r)*(1-w/r))
    absolute=circle_upper*geometric;relative=relative_circle_upper*geometric
    row={'eta':eta,'w':str(w),'terms':N,'radius':rs,'circle_root_upper':str(circle_upper),'relative_root_upper':str(relative_circle_upper),'exact_fourth_power_enclosures_checked':True,'absolute_tail_over_abs_h0_upper':str(absolute),'absolute_tail_over_abs_h0_decimal':display(absolute),'relative_error_upper':str(relative),'relative_error_decimal':display(relative),'scope':'exact-series truncation only; excludes decimal arithmetic/master evaluation rounding'}
    rows.append(row);candidates.append(row)
   selected=min(candidates,key=lambda x:Q(x['relative_error_upper']))
   best.append(selected)
 result={'schema':1,'status':'source-based scalar truncation bound evaluated with exact rational upper enclosures','source_scope':'N=1 original scalar only; positive parameter/compact measure after the fixed-dimensional common Gamma and measure prefactor','proof':'Cauchy circle bound plus real-axis positive lower kernel bound, with exact upward fourth-root enclosure','plan':bind(B/'bound-plan.json'),'source':bind(Path(__file__)),'input_bindings':plan['source_bindings'],'rows':rows,'best_fixed_grid_bounds':best,'reference_values_read':False,'pade_values_read':False,'endpoint_bound':False,'raised_bound':False,'numerical_rounding_bound_included':False}
 with (B/'scalar-bound.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'best_bounds':[{'eta':x['eta'],'terms':x['terms'],'radius':x['radius'],'relative_bound':x['relative_error_decimal']} for x in best]}))
if __name__=='__main__':main()
