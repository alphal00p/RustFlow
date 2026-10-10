#!/usr/bin/env python3
"""Exact source-expression specialization proving a fixed-ansatz obstruction."""
from pathlib import Path
from fractions import Fraction as Q
import json,ast,hashlib,re
R=Path(__file__).resolve().parent
source=R.parent/'2026-10-10-integrated-moment-expansion/scalar-source-kernel.json'
S=json.loads(source.read_text());a=S['audit'][0]
# Legitimate interior simplex/compact point: alpha0=alpha1=1/2, r0=r1=1/2,c=0.
# On C1: g00=g11=0,g01=r0*r1*(1-c)=1/4. This is a polynomial nonidentity witness,
# not a selected-invariant reduction, integration, or physical reference value.
vals={'integrated_moment_gauss::{}::alpha0':Q(1,2),'integrated_moment_gauss::{}::alpha1':Q(1,2),
      'rustflow_occupied::{}::q_0':Q(0),'rustflow_occupied::{}::q_3':Q(0),'rustflow_occupied::{}::q_1':Q(1,4),
      'integrated_moment::{}::eta':Q(-1,8)}
def value(raw):
 expr=raw.replace('^','**')
 for name,q in vals.items():expr=expr.replace(name,f'({q.numerator}/{q.denominator})')
 def go(n):
  if isinstance(n,ast.Constant) and isinstance(n.value,int):return Q(n.value)
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -go(n.operand)
  if isinstance(n,ast.BinOp):
   l,r=go(n.left),go(n.right)
   if isinstance(n.op,ast.Add):return l+r
   if isinstance(n.op,ast.Sub):return l-r
   if isinstance(n.op,ast.Mult):return l*r
   if isinstance(n.op,ast.Div):return l/r
   if isinstance(n.op,ast.Pow) and r.denominator==1:return l**r.numerator
  raise ValueError(ast.dump(n))
 return go(ast.parse(expr,mode='eval').body)
U,V,F=[value(a[k]) for k in ['U','V','F']];eta=vals['integrated_moment::{}::eta']; lam=Q(5,4)
assert U==1 and V==Q(1,8) and F==0 and eta+V/U**2==0
pure=value(S['physical_factors'][2]);assert pure==Q(-1,2) and eta-pure!=0
lead=-Q(32,1053)*eta**2*(eta+1)*(eta+4)
fall=Q(1)
for k in range(4):fall*=lam-k
assert lead!=0 and fall!=0
result={'status':'fixed_direct_eta_ansatz_obstructed_not_no_identity','identity_certified':False,
 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'candidate_sha256':hashlib.sha256((R/'target-0/recurrence-candidate.json').read_bytes()).hexdigest(),
 'kernel':'T = weight(r,c,x) * (eta+V/U^2)^(5/4)/(eta-physical_factor[2]); projective U=1',
 'allowed_certificate_denominator_power_per_irreducible_factor':2,
 'L_over_T_generic_pole_on_K':4,'div_R_plus_R_dlogT_max_pole_on_K':3,
 'witness':{'alpha':['1/2','1/2'],'r':['1/2','1/2'],'cos_angle':'0','eta':str(eta),'U':str(U),'V':str(V),'F':str(F),'other_denominator':str(eta-pure),'leading_deta4_coefficient':str(lead),'falling_kernel_power':str(fall),'nonzero_K_minus4_residue':str(lead*fall)},
 'proof':'At a generic smooth K=0 point the highest eta derivative contributes a4*falling(5/4,4)/K^4. Lower derivatives have pole<=3. Its numerator is not identically zero on K=0, witnessed by direct exact evaluation of the actual U/F/V source polynomials. R_i with irreducible K-pole<=2 gives derivative-plus-log-derivative pole<=3. Other distinct allowed factors are generically regular on K=0. Thus this fixed direct-eta rational certificate ansatz cannot certify this candidate.',
 'scope_limits':['Does not rule out higher-pole certificates, recurrence telescopers, other representations, boundary-coupled systems or an actual identity.','No telescoping linear system was run because this necessary condition fails.','No enlarged ansatz, new fit, endpoint answer or raised-target holdout was used.','Forcing must be treated as boundary terms; no independently prescribed interior forcing kernel with a compensating K^-4 pole was allowed.']}
(R/'certificate-pole-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
