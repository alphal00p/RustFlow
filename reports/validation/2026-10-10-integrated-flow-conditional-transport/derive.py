#!/usr/bin/env python3
"""Exact recurrence-to-forced-ODE conversion; no fitting or reference values."""
import json, hashlib, math
from fractions import Fraction as Q
from pathlib import Path
R=Path(__file__).resolve().parent
P=R.parent/'2026-10-10-integrated-flow-physical-fit'
D=R.parent/'2026-10-10-integrated-moment-expansion/fit-data'
candidate=json.loads((P/'target-0/recurrence-candidate.json').read_text())['candidate']
train=json.loads((D/'target-0-training.json').read_text())
held=json.loads((D/'target-0-holdout.json').read_text()) # This scalar holdout was already frozen-and-validated.
assert candidate['kind']=='recurrence' and candidate['status']=='finite_training_candidate_not_identity'
r=candidate['order']; d=candidate['degree']; beta=Q(train['channels'][0]['beta'])
p=[[Q(x) for x in candidate['coefficients'][j*(d+1):(j+1)*(d+1)]] for j in range(r+1)]
c=[Q(x) for x in train['channels'][0]['coefficients']+held['channels'][0]['coefficients']]
def peval(poly,n): return sum(a*n**k for k,a in enumerate(poly))
# L_beta=sum_j t^(r-j) p_j(theta-beta-j); coefficient storage [theta_power][t_power].
A=[[Q(0) for _ in range(r+1)] for _ in range(d+1)]
for j in range(r+1):
 for k in range(d+1):
  for ell in range(k+1): A[ell][r-j]+=p[j][k]*math.comb(k,ell)*(-beta-j)**(k-ell)
B=[Q(0) for _ in range(r)]
for j in range(r+1):
 for k in range(j): B[r-j+k]+=peval(p[j],Q(k-j))*c[k]
# Verify every available coefficient of L_beta[t^beta F] equals t^beta B.
for n in range(len(c)):
 lhs=sum(A[k][s]*(Q(n-s)+beta)**k*c[n-s] for k in range(d+1) for s in range(r+1) if n>=s)
 rhs=B[n] if n<r else Q(0)
 assert lhs==rhs,(n,lhs,rhs)
assert A[4]==[Q(-32,1053),Q(-160,1053),Q(-128,1053)]
# This checks only conversion consistency with known data, not the unknown infinite tail.
def poly(v): return '+'.join(f'({a})*(1/eta)^{i}' for i,a in enumerate(v) if a) or '0'
# y=(g,theta_t g,...theta_t^(d-1)g,h=t^beta); d/deta=-theta_t/eta.
m=[['0']*(d+1) for _ in range(d+1)]
for k in range(d-1): m[k][k+1]='-1/eta'
for k in range(d): m[d-1][k]=f'({poly(A[k])})/(eta*({poly(A[d])}))'
m[d-1][d]=f'-({poly(B)})/(eta*({poly(A[d])}))'
m[d][d]=f'-({beta})/eta'
record={'schema':'candidate-forced-companion.v1','identity_certified':False,'target_index':0,'q':1,'beta':str(beta),
 'recurrence_order':r,'differential_order':d,'p_n_ascending':[[str(x) for x in row] for row in p],
 'L_beta_theta_rows_t_ascending':[[str(x) for x in row] for row in A], 'startup_B_t_ascending':[str(x) for x in B],
 'equation':'sum_k A_k(t) theta_t^k g(t) = t^beta B(t), g=t^beta F, theta_t=t*d/dt',
 'conversion_checked_coefficients':len(c),'matrix_eta':m,
 'highest_theta_coefficient':'-32/1053*(1+t)*(1+4*t)',
 'leading_eta_derivative_coefficient':'-32/1053*eta^2*(eta+1)*(eta+4)',
 'ordinary_singularities':['0','-1','-4'],
 'normalization':'Native hard master and compact angular normalization remain factored; output is g only',
 'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [P/'target-0/recurrence-candidate.json',D/'target-0-training.json',D/'target-0-holdout.json']}}
(R/'companion.json').write_text(json.dumps(record,indent=2)+'\n')
(R/'scalar-coefficients.json').write_text(json.dumps({'beta':str(beta),'coefficients':[str(x) for x in c]},indent=2)+'\n')
print(json.dumps({'startup':[str(x) for x in B],'matrix_eta':m,'checked':len(c)},indent=2))
