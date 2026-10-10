"""Exact metadata-only matching. No moment or numerical reference is evaluated."""
import argparse, ast, hashlib, json, re
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
INPUT=ROOT/'examples/finite_density/massless_three_loop_chain.json'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(a,b):
 out=dict(a)
 for k,v in b.items(): out[k]=out.get(k,Fraction(0))+v
 return {k:v for k,v in out.items() if v}
def mul(a,b):
 out={}
 for x,v in a.items():
  for y,w in b.items():
   z=tuple(sorted(x+y));out[z]=out.get(z,Fraction(0))+v*w
 return {k:v for k,v in out.items() if v}
def parse(s):
 s=re.sub(r'rustflow_occupied::\{\}::q_(\d+)',r'q\1',s)
 s=re.sub(r'rustflow_density::\{\}::line_mass_squared_(\d+)',r'm\1',s)
 s=s.replace('integrated_moment::{}::eta','eta').replace('^','**')
 def walk(n):
  if isinstance(n,ast.Constant) and type(n.value) is int: return {():Fraction(n.value)} if n.value else {}
  if isinstance(n,ast.Name) and re.fullmatch(r'q[0-8]|m[0-4]|eta',n.id): return {(n.id,):Fraction(1)}
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): return {k:-v for k,v in walk(n.operand).items()}
  if isinstance(n,ast.BinOp):
   a,b=walk(n.left),walk(n.right)
   if isinstance(n.op,ast.Add):return add(a,b)
   if isinstance(n.op,ast.Sub):return add(a,{k:-v for k,v in b.items()})
   if isinstance(n.op,ast.Mult):return mul(a,b)
  raise ValueError('not an admitted exact polynomial: '+ast.dump(n))
 return walk(ast.parse(s,mode='eval').body)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('metadata',type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 assert not a.output.exists()
 d=json.loads(a.metadata.read_text()); independent=json.loads((BASE/'basis-crosscheck.json').read_text())
 assert sha(INPUT)==independent['input_sha256']
 assert d['cut_slots']==[0,3] and d['shifted_slots']==[1,2,4]
 assert d['routing']==[['1','0','0'],['0','0','1'],['1','-1','0']] and d['routing_determinant']=='1'
 expected=['q0','q5','q0+q3-2*q1','q3','q0+q3+q5+2*q4-2*q1-2*q2','q2','q6','q8','q6-q7','1-q6','q6','1-q7','q7']
 for slot,value in enumerate(expected):
  assert parse(d['physical_factors'][slot])==parse(value),slot
  assert parse(d['independent_factors'][slot])==parse(value+(f'-m{slot}' if slot<5 else '')),slot
  assert parse(d['eta_factors'][slot])==parse(value+('-eta' if slot in [1,2,4] else '')),slot
 assert d['roles']==['RequiredCut','Ordinary','Ordinary','RequiredCut']+['Ordinary']*5+['Occupation']*4
 terms=[[(independent['scalar_powers'],'-1')],[(row['indices'],'-1') for row in independent['raised_numerator_terms_before_wick']]]
 target_maps=[]
 for target,rows in zip(d['physical_targets'],terms,strict=True):
  wanted={tuple(indices+[0]*4):coefficient for indices,coefficient in rows}
  actual={tuple(row['indices']):row['coefficient'] for row in target}
  assert len(actual)==len(target) and actual==wanted
  target_maps.append(actual)
 assert d['physical_targets']==d['independent_targets']
 assert d['shells']==[{'chemical_potential':'1','loop_index':0,'lower_slot':10,'mass_squared':'0','physical_slot':0,'upper_slot':9}, {'chemical_potential':'1','loop_index':1,'lower_slot':12,'mass_squared':'0','physical_slot':3,'upper_slot':11}]
 # Algebraic raised-jet and symmetry checks at exact rational sample points.
 checks=0
 for r in [Fraction(1,5),Fraction(2,3),Fraction(1)]:
  for s in [Fraction(1,7),Fraction(3,4)]:
   for t in [Fraction(1,9),Fraction(2,5),Fraction(7,8)]:
    N=-r*r/2+r*s*(Fraction(1,2)+t); Na=-1+s/(2*r); ha=-1+s/r
    assert N/(2*r*r)-Na==Fraction(3,4)+(s/r)*(t/2-Fraction(1,4))
    assert N*ha==r*r/2-r*s*(t+1)+s*s*(t+Fraction(1,2))
    Nswap=-s*s/2+s*r*(Fraction(1,2)+t)
    assert (N*ha+Nswap*(-1+r/s))/2==(t+1)*(r-s)**2/2
    checks+=3
 result={'schema':1,'status':'passed','scope':'exact actual prepared conversion and fixed-original numerator jet audit; no numerical integration',
  'metadata':str(a.metadata.resolve()),'metadata_sha256':sha(a.metadata),'input_sha256':sha(INPUT),
  'basis_crosscheck_sha256':sha(BASE/'basis-crosscheck.json'),'audit_source_sha256':sha(__file__),
  'actual_input_identity':d['input_identity'],'factor_checks':39,'target_terms':3,'rational_jet_checks':checks,
  'special_case_N_eta_equals_N':True,'generic_equality_claim':False,
  'mass_conversion_coefficient_derivatives':'zero for this fixture only',
  'wick_signs':'scalar (-1)^5=-1; raised (-1)^6 times degree-two numerator Wick sign -1=-1',
  'normalization_bridge':'Three ordinary Minkowski factors Wick-rotate with three minus signs; combined with scalar coefficient -1 the Euclidean scalar kernel is positive. The Minkowski C2 is +d/da C1, while the original Euclidean numerator is the negative of its continued polynomial, reproducing -d/da at fixed original N.',
  'upper_surface':'positive from -d/da Theta(mu-E1); C2H0 acts on support despite original upper indices zero',
  'producer_coefficients_read':False,'reference_values_generated':False}
 a.output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
