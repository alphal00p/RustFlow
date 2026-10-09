# Independent massive sunset reference at D=12/5

Validation-only derivation; no oracle answers, AMF reductions, or AMF boundary
values enter the integrands below. Target definitions are the native example:
q1 and q2 carry the same chemical potential mu=1, their masses are m1=m2=1/2,
and the neutral q1-q2 line has mass M=1. Original shifted Euclidean targets are
S=1/(rho1 rho2 rho3) and R=(g12+u1*u2)/(rho1^2 rho2 rho3).
D=12/5 and d=D-1=7/5. Both targets are UV convergent; all masses are positive.
The auxiliary numerator target with powers111 is used only as a differentiation
identity; its divergent vacuum part is never numerically integrated.

Let a=m1², b=m2², c=M², E_i=sqrt(r_i²+m_i²),
A_d=2*pi^(d/2)/(Gamma(d/2)*(2*pi)^d), and rFi=sqrt(mu²-m_i²).
All displayed values use unscaled Euclidean measures. MSbar and requested scale
normalization may be applied only after forming the complete sum.

## Vacuum terms, direct Schwinger Gaussian integration

On the simplex x1+x2+x3=1, write U=x1*x2+x1*x3+x2*x3 and
F=a*x1+b*x2+c*x3. Then

V_S = Gamma(3-D)/(4*pi)^D * int_simplex U^(-D/2)*F^(D-3) dx,
V_R = (D+1)/2 * Gamma(3-D)/(4*pi)^D
      * int_simplex x1*x3*U^(-D/2-1)*F^(D-3) dx.

The numerator Gaussian covariance is
<g12+u1*u2>=(D+1)*alpha3/(2*(alpha1*alpha2+alpha1*alpha3+alpha2*alpha3)).
The alpha1 factor from the raised propagator gives V_R directly, avoiding the
UV-divergent unraised numerator vacuum integral.

For each of the six permutations(j,k,l), set
x_j=1/h,x_k=t/h,x_l=t*w/h,h=1+t+t*w, with t=s^5.
The scalar kernel on (s,w) in [0,1]^2 is
5*s^3*(1+w+s^5*w)^(-D/2)*(m_j²+s^5*m_k²+s^5*w*m_l²)^(D-3).
For V_R multiply it by x1*x3/U, expressed through these ratios; this ratio is
1/(1+w+t*w), w/(1+w+t*w), or t*w/(1+w+t*w), according to the permutation.
All kernels are smooth at the unit-square boundary.

## One-cut terms, Feynman-parameter integration only in the reference

Define K=Gamma(2-D/2)/(4*pi)^(D/2), alpha=D/2-2,
F1(x)=x*b+(1-x)*c-x*(1-x)*a,
F2(x)=x*a+(1-x)*c-x*(1-x)*b.

B_i = K*int_0^1 F_i(x)^alpha dx,
A_i = K*int_0^1 (1-x)*F_i(x)^alpha dx,
A1_a = K*(2-D/2)*int_0^1 x*(1-x)^2*F1(x)^(alpha-1) dx,
minus_A2_a = K*(2-D/2)*int_0^1 x*(1-x)*F2(x)^(alpha-1) dx.

W_i=A_d/2*int_0^rFi r^(d-1)/E_i dr,
T_i=A_d/2*int_0^rFi r^(d-1)*(m_i²+E_i²)/E_i dr,
T1_a=A_d/4*int_0^rF1 r^(d-1)*(3*E1²-a)/E1³ dr
     - A_d*rF1^(d-2)*(a+mu²)/(4*mu).

The original one-cut numerator N111 contributions are +A1*T1 and +A2*T2.
Applying -d/da at fixed original momentum numerator, before equal-mass
specialization, yields

C1_S=-W1*B1, C2_S=-W2*B2,
C1_R=-A1_a*T1-A1*T1_a,
C2_R=minus_A2_a*T2.

The final term in T1_a is the explicit upper Fermi surface. It must be reported
separately and remain nonzero at the chosen point. The scalar bubble here is a
validation-only parameter integral, not an evaluator terminal or production
analytic bubble formula.

## Fully occupied compact integrals

With z the spatial angle cosine, put
h=c-a-b+2*E1*E2-2*r1*r2*z,
n=-2*E1*E2+r1*r2*z,
h_a=-1+E2/E1, n_a=-E2/E1.
The bound h>=c-(m1-m2)^2>0 is uniform on support.
Use normalized angular average
< f >_d=Gamma(d/2)/(sqrt(pi)*Gamma((d-1)/2))
        * int_-1^1 (1-z²)^((d-3)/2)*f(z) dz.

C12_S=A_d²*int dr1 dr2 r1^(d-1)*r2^(d-1)/(4*E1*E2) * <1/h>_d.

The raised fixed-original-numerator contribution is the sum of bulk and surface:

C12_R_bulk=A_d²*int dr1 dr2 r1^(d-1)*r2^(d-1)/(4*E1*E2)
          * < n/(2*E1²*h) + E2/(E1*h) + n*h_a/h² >_d,
C12_R_surface=A_d²*rF1^(d-2)/(4*mu)
          * int_0^rF2 dr2 r2^(d-1)/(2*E2)
          * < n/h evaluated at E1=mu,r1=rF1 >_d.

The surface sign follows from -dW1/da, including delta(mu-E1)/(4*E1²).
Differentiate the complete on-shell expression derived from the original
momentum numerator, retaining E1(a), n(a), h(a), and moving support; do not freeze
mass-dependent shell substitutions.
The derivative formulas provide an independent check of native raised C_n/H_n
normalization without using the native distribution integration owner.

For radial integrals use r=rF*s^5, so r^(d-1)dr=5*rF^d*s^6 ds.
For the angular integral split z=+(1-s^5) and z=-(1-s^5). In either branch,
(1-z²)^(-4/5)*abs(dz)=5*(2-s^5)^(-4/5) ds.
These maps remove all algebraic endpoint weights at the fixed D. Quadrature
nodes/weights can be shared across scalar, raised, bulk and surface integrands.

## Numerical owner and verification plan

RustFlow src/numeric.rs provides fixed-working-precision MPFR/Symbolica
operations, gamma, powers and comparisons. The pinned Numerica integration API
lib/numerica/src/numerical_integration.rs provides adaptive sampling grids, not
a deterministic multiprecision quadrature suitable for ten-digit references.
No deterministic quadrature owner was found in RustFlow or the pinned HEPKit
crates. A small validation-only Gauss-Legendre implementation can therefore use
RustFlow Precision arithmetic, with Legendre nodes found by Newton recurrence,
and independently validate weight sum, polynomial moments and n-vs-2n results.
It belongs in tests/tools, never in the production finite-density evaluation.

Run independent tensor orders 32,48,64 (increase if needed), then repeat at
40/60 working decimal digits. Record every cut and each explicit surface,
node count, working precision, wall time, memory and absolute/relative changes.
At each refinement form S=V_S+C1_S+C2_S+C12_S and
R=V_R+C1_R+C2_R+C12_R_bulk+C12_R_surface before comparison.
Require at least1e-12 absolute and1e-10 relative stability, tightened if
cancellations make component precision insufficient. Preserve the individual
contributions to catch an omitted nonzero vacuum or Fermi surface.
A further independent Schwinger vacuum parametrization or numerical mass
derivative check may be used to cross-check these validation kernels; it is not
part of the AMF feature. The reference must remain unlabelled as certified until
its actual refinement results pass.
### Extension toward the Laurent point

The validation-only `reference_at_dimension` integrator also represents the
massive integral near D=4−2 epsilon. In each of six Cheng–Wu sectors take
`x_j=1, x_k=t, x_l=t*u`, with t,u in [0,1]. The scalar vacuum integrand is
`t^(-1+epsilon) B(t,u)`, where

```text
B(t,u) = (1+u+t*u)^(-D/2)
         * (m_j²+t*(m_k²+u*m_l²))^(D-3).
```

Subtract B(0,u) under the t integral and add its exact integral divided by
epsilon. With `J(D)=(2^(1-D/2)-1)/(1-D/2)`, the scalar counterterm is
`2*sum_j (m_j²)^(D-3)*J(D)/epsilon`. The raised numerator's Gaussian weight
is `x_0*x_2/U`. Its two sector limits sum to one when j is not 1 and vanish
when j is 1, giving `sum_(j!=1) (m_j²)^(D-3)*J(D)/epsilon`, followed by the
same overall `(D+1)/2` tensor factor used in the convergent reference.
Both include `Gamma(3-D)/(4*pi)^D`. This gives an independent double-pole
bound; single-cut bubbles have only `Gamma(epsilon)`, and the massive
double-cut terms are analytic near epsilon zero.

The subtracted remainder is integrable for Re epsilon > −1. The displayed
angular representation requires D>2; actual gamma poles are excluded from
sampling. Fifth-power endpoint maps leave powers `s^(4+5*epsilon)` in the
vacuum remainder and `s^(4-5*epsilon)` in the angle integral. Their noninteger
endpoint behavior requires fresh quadrature-order checks near epsilon zero.
The reference generator independently varies quadrature order, the exact
epsilon grid, and working precision. Exact rational Lagrange weights extract
the first three Taylor coefficients of epsilon² times each contribution.
It uses neither the production Laurent fitter nor oracle coefficients.
The saved [Laurent reference report](../reports/validation/2026-10-09-finite-density-native-assembly/independent-laurent-reference/README.md)
passes all 108 coefficient refinement checks and six separately derived
analytic ultraviolet-residue checks. Its largest observed relative change is
3.481e−17 under quadrature refinement. This is empirical reference accuracy,
not a rigorous interval bound or a finite-density AMF prediction.

### Comparing complete native predictions

The ignored native tests `complete_massive_sunset_assembles_vacuum_and_all_occupied_sectors`
and `complete_massive_sunset_laurent_refinement` save predictions before reading
any reference. Only after a successful run, compare the saved artifacts with:

```sh
python3 tools/finite_density/compare_complete_massive_reference.py sample \
  --predictions PATH_TO_SAMPLE_RUN --output PATH_TO_SAMPLE_COMPARISON.json
python3 tools/finite_density/compare_complete_massive_reference.py laurent \
  --predictions PATH_TO_LAURENT_RUN --output PATH_TO_LAURENT_COMPARISON.json
```

The sample test must use `RUSTFLOW_DENSITY_FLOW_EPSILON=4/5` and save all four
precision/order/start configurations `(18,60,8)`, `(28,60,8)`, `(28,80,8)` and
`(28,80,12)`. Laurent mode requires those four configurations at epsilon-grid
denominator 1000 **plus** `(28,80,12)` at denominator 2000. The latter file is
`prediction-28-80-12-grid-2000.json`; the first four retain their original
`prediction-DIGITS-ORDER-START.json` names. Each Laurent record must declare
its `epsilon_grid_denominator` explicitly. The fifth comparison changes only
the grid denominator.

All profiles must have identical `guard_digits` and `search_frontier_sectors`
metadata within each comparison run, so these settings cannot silently change
along with the independently varied precision, order, start scale or grid.
The exact input definition must match the independent reference. Sample
comparisons include every vacuum/cut contribution,
the explicit nonzero raised surface terms in the reference, and the assembled
sum. Laurent comparisons cover the assembled coefficients at powers −2, −1, 0.
All compared inputs are hashed; this utility cannot generate predictions or
fill missing sectors. Availability of the command is not a comparison result.
