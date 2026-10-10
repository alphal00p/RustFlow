"""Saved rational series evaluation. This process never reads reference values."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,re,sys,time
import mpmath as mp
BASE=Path(__file__).resolve().parent
PATTERN=re.compile(r'-?\d+(?:/\d+)?\Z')
def rational(s):
 assert isinstance(s,str) and PATTERN.fullmatch(s),f'nonrational coefficient {s!r}'
 return Fraction(s)
def bind(p):
 p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def value(q):return mp.mpf(q.numerator)/q.denominator
def complex_value(v):return mp.mpc(mp.mpf(v['re']),mp.mpf(v['im']))
def printable(v):return {'re':mp.nstr(v.real,mp.mp.dps),'im':mp.nstr(v.imag,mp.mp.dps)}
def validate(data,normalization,plan):
 assert data['schema']==plan['producer_schema'] and data['q']==1 and data['variable']=='z=1/eta' and data['tail']=='unknown'
 assert data['dimension']==normalization['dimension']==plan['dimension']
 assert data['input_identity']==normalization['input_identity']
 assert normalization['cut_slots']==plan['cut_slots']
 assert int(data['known_prefix_length'])>=max(plan['profiles']['known_terms']),'insufficient saved prefix; unknown coefficients cannot be invented'
 exponents={}
 for audit in data['audit']:
  assert audit['full_off_shell_replay'] is True
  t=audit['target'];e=rational(audit['eta_leading_exponent_at_D'])
  if t in exponents:assert exponents[t]==e,'multiple exponents need separate channels'
  exponents[t]=e
 assert set(exponents)==set(plan['original_targets'])
 coefficients={};coverage={t:set() for t in exponents}
 for row in data['coefficients']:
  target=row['target'];n=row['n'];channel=row['master']
  assert target in exponents and isinstance(n,int) and 0<=n<data['known_prefix_length']
  key=(target,n,channel);assert key not in coefficients,'duplicate coefficient channel'
  coefficients[key]=rational(row['coefficient']);coverage[target].add(n)
 assert all(ns==set(range(data['known_prefix_length'])) for ns in coverage.values()),'missing coefficients must be explicit zeros'
 channels={key[2] for key in coefficients}
 assert {p['digits'] for p in normalization['profiles']}==set(plan['profiles']['digits'])
 for profile in normalization['profiles']:
  assert len(profile['masters'])==len(channels) and {x['channel'] for x in profile['masters']}==channels
  assert [x['condition'] for x in profile['checked_conditions']]==data['conditions']
  assert all(Fraction(x['value']['re'])!=0 or Fraction(x['value']['im'])!=0 for x in profile['checked_conditions'])
 return exponents,coefficients

def main():
 if len(sys.argv)!=3:raise SystemExit('sum-saved-prefix.py NORMALIZER_RUN_DIRECTORY OUTPUT_JSON')
 directory=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();assert not out.exists()
 source=directory/'coefficients.json';native=directory/'normalization.json';runpath=directory/'run.json';planpath=BASE/'plan.json'
 run=json.loads(runpath.read_text());assert run['exit_code']==0 and run['inputs_unchanged'] and run['executable_unchanged'] and run['dependencies_unchanged']
 assert bind(native)==run['output'];assert bind(source)['sha256']==run['inputs']['coefficients']['sha256']
 data=json.loads(source.read_text());normalization=json.loads(native.read_text());plan=json.loads(planpath.read_text())
 exponents,coefficients=validate(data,normalization,plan)
 inputs={name:bind(p) for name,p in [('saved_coefficients',source),('native_normalization',native),('normalization_run',runpath),('plan',planpath),('adapter',Path(__file__))]}
 started=time.monotonic();rows=[]
 for profile in normalization['profiles']:
  digits=profile['digits'];mp.mp.dps=digits
  native_factor=complex_value(profile['native_measure_to_euclidean'])*complex_value(profile['compact_normalization'])
  masters={x['channel']:complex_value(x['value']) for x in profile['masters']}
  for eta_string in plan['profiles']['eta']:
   eta=rational(eta_string);assert eta>0
   for terms in plan['profiles']['known_terms']:
    for target,exponent in sorted(exponents.items()):
     exact_channel_sums={channel:sum((coefficient/eta**n for (t,n,c),coefficient in coefficients.items() if t==target and n<terms and c==channel),Fraction(0)) for channel in masters}
     physical=native_factor*mp.power(value(eta),value(exponent))*mp.fsum(value(q)*masters[c] for c,q in exact_channel_sums.items())
     rows.append({'target':target,'eta':eta_string,'terms':terms,'digits':digits,'leading_eta_exponent':str(exponent),'value':printable(physical),'exact_channel_sums':{c:str(q) for c,q in exact_channel_sums.items()},'status':'saved-prefix partial sum; no reference read, no endpoint or ODE assertion'})
 result={'schema':1,'scope':'Exact saved coefficient channel sums with native Gaussian/measure normalization; no coefficient constructed from reference','input_identity':data['input_identity'],'cut_slots':plan['cut_slots'],'dimension':plan['dimension'],'known_prefix_length':data['known_prefix_length'],'unknown_tail':True,'inputs':inputs,'rows':rows,'wall_seconds':time.monotonic()-started,'numerical_reference_reads':0}
 assert all(bind(v['path'])==v for v in inputs.values())
 with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'saved':str(out),'rows':len(rows),'sha256':bind(out)['sha256'],'reference_reads':0}))
if __name__=='__main__':main()
