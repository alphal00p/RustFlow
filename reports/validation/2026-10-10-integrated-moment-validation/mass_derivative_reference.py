"""One-sided physical mass diagnostic of the fixed original numerator.

No producer values are read. The moving radius and fixed-radius companion
separate the positive Fermi surface from the bulk mass derivative.
"""
import argparse,hashlib,json,time
from pathlib import Path
import mpmath as mp
import finite_eta_reference as reference

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dec(v):return mp.nstr(v,mp.mp.dps)
def mixed_unraised(D,eta,mu,mass,n,moving):
 alpha=(D-2)/2;nu=D/2-2
 radius=mp.sqrt(mu*mu-mass) if moving else mu
 radial1=reference.radial_rule(n,D-4,radius)
 radial2=reference.radial_rule(n,D-4,mu)
 angular=reference.angular_rule(n,alpha);parameter=reference.parameter_rule(n)
 area=2*mp.pi**((D-1)/2)/mp.gamma((D-1)/2)
 A=area/(2*mp.pi)**(D-1);gaussian=(4*mp.pi)**(-D/2)*mp.gamma(2-D/2)
 rows=[]
 for r1,w1 in radial1:
  E1=mp.sqrt(r1*r1+mass)
  for r2,w2 in radial2:
   for t,wt in angular:
    z=1-2*t
    h=-mass+2*(E1*r2-r1*r2*z)
    assert eta+h>0
    N=(h-mass)/4-E1*E1/2+E1*r2/2
    K,_=reference.kernel(h,eta,nu,gaussian,parameter)
    rows.append(w1*w2*wt*(r1*r1/E1)*r2*K*N)
 return A*A/4*mp.fsum(rows)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--receipt',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 receipt=json.loads(a.receipt.read_text());profile=receipt['profile']
 assert receipt['prediction_artifacts_frozen_before_reference'] is True
 assert receipt['reference_evaluator_sha256']==sha(__file__)
 assert receipt['support_evaluator_sha256']==sha(reference.__file__)
 assert profile['dimension']=='6.5' and profile['eta']=='2' and profile['mu']=='1'
 assert profile['digits']==50 and profile['nodes'] in [24,32]
 assert profile['physical_mass_squared'] in ['0','1/256','1/1024','1/4096']
 assert not a.output.exists()
 mp.mp.dps=profile['digits'];started=time.monotonic()
 D,eta,mu=map(mp.mpf,['6.5','2','1'])
 pieces=profile['physical_mass_squared'].split('/');mass=mp.mpf(pieces[0])/mp.mpf(pieces[1] if len(pieces)==2 else 1)
 assert 0<=mass<min(mu*mu,eta)/16
 fixed=mixed_unraised(D,eta,mu,mass,profile['nodes'],False)
 moving=fixed if mass==0 else mixed_unraised(D,eta,mu,mass,profile['nodes'],True)
 result={'schema':1,'status':'physical-mass diagnostic quadrature completed','profile':profile,
  'fixed_radius_unraised_original_numerator':dec(fixed),'moving_radius_unraised_original_numerator':dec(moving),
  'moving_radius':dec(mp.sqrt(1-mass)),'elapsed_seconds':time.monotonic()-started,
  'receipt_sha256':sha(a.receipt),'evaluator_sha256':sha(__file__),'support_evaluator_sha256':sha(reference.__file__),
  'scope':'Validation only. Original g1_2+u1*u2 and mu fixed; physical mass varied before shell substitution. No producer coefficients read.'}
 tmp=a.output.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(a.output)
if __name__=='__main__':main()
