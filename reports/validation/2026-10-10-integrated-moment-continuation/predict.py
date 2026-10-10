"""Exact same-prefix composition and bounded numerical diagnostics; no reference reads."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,re,sys,time
import mpmath as mp
B=Path(__file__).resolve().parent
PATTERN=re.compile(r'-?\d+(?:/\d+)?\Z')
def q(s):assert isinstance(s,str) and PATTERN.fullmatch(s);return Q(s)
def bind(p):p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def val(x):return mp.mpf(x.numerator)/x.denominator
def cval(x):return mp.mpc(x['re'],x['im'])
def text(x):return {'re':mp.nstr(mp.re(x),mp.mp.dps),'im':mp.nstr(mp.im(x),mp.mp.dps)}
def binomial(a,n):
 result=Q(1)
 for j in range(n):result*= (a-j)/(j+1)
 return result

def main():
 started=time.monotonic();plan=json.loads((B/'plan.json').read_text());assert not (B/'predictions.json').exists()
 for item in plan['producer'].values():
  if isinstance(item,dict) and 'path' in item:assert bind(item['path'])==item
 data=json.loads(Path(plan['producer']['coefficients']['path']).read_text());norm=json.loads(Path(plan['producer']['normalization']['path']).read_text())
 assert data['schema']=='integrated-gaussian-rational-prefix.v1' and data['q']==1 and data['tail']=='unknown' and data['dimension']=='13/2'
 assert data['known_prefix_length']==plan['producer']['coefficient_count']==64 and data['input_identity']==norm['input_identity']
 alpha=q(plan['transformation']['alpha']);R=q(plan['transformation']['R']);targets=plan['producer']['targets'];nmax=64
 assert all(q(x['eta_leading_exponent_at_D'])==alpha and x['full_off_shell_replay'] is True for x in data['audit'])
 channels=sorted({x['master'] for x in data['coefficients']});coefs={}
 for row in data['coefficients']:
  key=(row['target'],row['master'],row['n']);assert key not in coefs;coefs[key]=q(row['coefficient'])
 assert set(coefs)==set(itertools.product(targets,channels,range(nmax)))
 h={};roundtrips=[]
 for target,channel in itertools.product(targets,channels):
  a=[coefs[target,channel,n] for n in range(nmax)]
  hs=[sum((a[n]*R**(-n)*(-1)**(k-n)*binomial(alpha-n,k-n) for n in range(k+1)),Q(0)) for k in range(nmax)]
  inverse=[R**k*sum((hs[m]*binomial(alpha-m,k-m) for m in range(k+1)),Q(0)) for k in range(nmax)]
  assert inverse==a;roundtrips.append({'target':target,'channel':channel,'exact_coefficients_checked':nmax,'passed':True});h[target,channel]=hs
 exact={'schema':1,'scope':'Exact finite lower-triangular composition of the same 64 coefficients; unknown tail remains unknown','plan':bind(B/'plan.json'),'producer':plan['producer']['coefficients'],'transformed':[{'target':t,'channel':c,'coefficients':[str(x) for x in hs]} for (t,c),hs in h.items()],'inverse_composition_checks':roundtrips}
 with (B/'transformed-coefficients.json').open('x') as f:f.write(json.dumps(exact,indent=2)+'\n')
 rows=[];pades=[];failures=[]
 for profile in norm['profiles']:
  digits=profile['digits'];assert digits in plan['profiles']['digits'];mp.mp.dps=digits
  masters={x['channel']:cval(x['value']) for x in profile['masters']};assert set(masters)==set(channels)
  factor=cval(profile['native_measure_to_euclidean'])*cval(profile['compact_normalization'])
  for target in targets:
   coeff=[mp.fsum(val(h[target,c][k])*masters[c] for c in channels) for k in range(nmax)]
   assert all(mp.im(x)==0 for x in coeff);coeff=[mp.re(x) for x in coeff]
   for eta_string in plan['profiles']['eta']:
    eta=q(eta_string);w=R/(eta+R);physical=factor*mp.power(val(R),val(alpha))*mp.power(val(w),-val(alpha))
    for terms in plan['profiles']['partial_sum_terms']:
     # Sum each channel in Q, then apply the independently evaluated hard period.
     hs={c:sum((h[target,c][k]*w**k for k in range(terms)),Q(0)) for c in channels}
     result=physical*mp.fsum(val(v)*masters[c] for c,v in hs.items())
     rows.append({'method':'partial_sum','target':target,'eta':eta_string,'w':str(w),'terms':terms,'digits':digits,'value':text(result)})
   for degree in plan['profiles']['pade_degrees']:
    L,M=degree;fit_terms=L+M+1
    try:
     numerator,denominator=mp.pade(coeff[:fit_terms],L,M)
     residuals=[]
     for k in range(nmax):
      terms=[denominator[j]*coeff[k-j] for j in range(min(M,k)+1)]
      actual=mp.fsum(terms)-(numerator[k] if k<=L else 0)
      scale=mp.fsum(abs(x) for x in terms)+(abs(numerator[k]) if k<=L else 0)
      residuals.append({'n':k,'normalized_residual':mp.nstr(abs(actual)/scale if scale else abs(actual),mp.mp.dps),'used_for_fit':k<fit_terms})
     record={'target':target,'digits':digits,'degrees':degree,'fit_terms':fit_terms,'remaining_prefix_terms':nmax-fit_terms,'numerator':[mp.nstr(x,mp.mp.dps) for x in numerator],'denominator':[mp.nstr(x,mp.mp.dps) for x in denominator],'prefix_residuals':residuals,'status':'numerical rational approximant; no all-order identity or pole-free path proof'}
     pades.append(record)
     for eta_string in plan['profiles']['eta']:
      eta=q(eta_string);w=R/(eta+R);wx=val(w);den=mp.polyval(list(reversed(denominator)),wx)
      assert den!=0,'Pade denominator vanishes at evaluation point'
      pvalue=mp.polyval(list(reversed(numerator)),wx)/den
      physical=factor*mp.power(val(R),val(alpha))*mp.power(wx,-val(alpha))
      denom_norm=mp.fsum(abs(x)*wx**k for k,x in enumerate(denominator))
      rows.append({'method':'pade','target':target,'eta':eta_string,'w':str(w),'degrees':degree,'digits':digits,'value':text(physical*pvalue),'denominator':mp.nstr(den,mp.mp.dps),'relative_denominator':mp.nstr(abs(den)/denom_norm,mp.mp.dps)})
    except (ValueError,ZeroDivisionError,AssertionError) as e:
     failures.append({'target':target,'digits':digits,'degrees':degree,'error':str(e),'status':'preserved failed rational approximant; no replacement profile'})
 result={'schema':1,'scope':'Both targets, same frozen 64-term prefix; no reference or ODE input','inputs':{'plan':bind(B/'plan.json'),'source':bind(Path(__file__)),'coefficients':plan['producer']['coefficients'],'normalization':plan['producer']['normalization'],'transformed':bind(B/'transformed-coefficients.json')},'rows':rows,'pade_models':pades,'failures':failures,'reference_reads':0,'unknown_tail':True,'wall_seconds':time.monotonic()-started,'mpmath_version':mp.__version__}
 for item in result['inputs'].values():assert bind(item['path'])==item
 with (B/'predictions.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'saved':bind(B/'predictions.json'),'rows':len(rows),'pade_models':len(pades),'failures':len(failures),'wall_seconds':result['wall_seconds'],'reference_reads':0}))
if __name__=='__main__':main()
