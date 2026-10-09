# Raised prism: a concrete local subtraction representation

Validation-only construction. The complete raised finite-density prism
reference remains ungenerated. A bounded independent upstream run now supplies
one profile of five virtual coefficient functions at epsilon=1/5, independently
checked against Gamma functions; it does not supply the missing complete
Laurent amplitude. See the partial-reference report linked below. This does not establish the original common thermal endpoint
prescription. It supplies implementable UV and isolated-collinear subtraction
operators for the three-cut kernel in
[the remaining-reference derivation](finite-density-missing-references.md), derived from the original fixed
numerator. The two-cut coefficient functions and overlapping endpoint sectors
must still be included before a complete independent reference exists.

## Independent neutral indices and the virtual UV piece

Assign the three neutral triangle lines (physical prism slots 2, 6, 8) powers
`beta_i=1+lambda_i`, i=0,1,2, and introduce the same line regulators in every
cut contribution before splitting the amplitude. Write

```
A = beta_0+beta_1+beta_2, Lambda=A-3, nu=A-D/2,
w_beta(x) = product x_i^(beta_i-1) / product Gamma(beta_i),
F = x0*x1*h + x0*x2*b + x1*x2*c,
h=h12, b=h13, c=h23.
```

An optional fixed reference scale multiplies the regulated amplitude by
`kappa^(2 Lambda)`; it is an analytic-index scale, not an AMF deformation.
For the supplemental numerator, which is at most linear in the remaining
virtual momentum, the exact triangle parameter numerator is

```
N(a)=(b-a)^2/4 + (h-a)/2 * [x1*(h+b-c)/2+x2*b].
T_beta = kappa^(2 Lambda)/(4*pi)^(D/2)
         * integral_simplex w_beta(x) N(a) Gamma(nu) F^(-nu).
```

Here a is the independent squared mass of occupied q1, and h,b depend on a.
This follows from the same Gaussian mean as at unit powers; no covariance
term occurs for a numerator linear in the virtual loop.

The angular raised pole requires Re D>6 before subtraction. Increasing A
separates that requirement from the triangle UV bound Re D<2 Re A. For
example, the isolated generic collinear/UV patch has an open overlap at
6<Re D<6+2 Re Lambda with Re Lambda>0. This is a local overlap statement;
it is not a claim that all soft and overlapping parameter regions converge
there. A complete initial domain must be checked for every resolved sector.

For any fixed positive scale S (with units of mass squared), use the exact
Schwinger identity

```
Gamma(nu) F^(-nu)
 = Gamma(nu) S^(-nu)
   + integral_0^infinity ds s^(nu-1) [exp(-s F)-exp(-s S)].
```

The difference is integrable at s=0 for Re nu>-1, and at infinity for F,S>0.
Equivalently, evaluate the difference as
`Gamma(nu)*(F^(-nu)-S^(-nu))`, using a cancellation-safe exponential-minus-one
formula near nu=0. Its nu=0 value is `log(S/F)`.
The exact UV counterterm can integrate the parameter simplex analytically:

```
N_UV(a)=(b-a)^2/4
 +(h-a)/2 * [(beta_1/A)*(h+b-c)/2+(beta_2/A)*b],
T_UV = kappa^(2 Lambda)/(4*pi)^(D/2)
       * Gamma(nu)/Gamma(A) * S^(-nu) * N_UV(a).
```

The remainder is the same simplex integral with the subtracted Schwinger
kernel. The identity is exact where the original parameter integral converges
and then defines its meromorphic continuation. The pole at nu=0 is isolated
without introducing a fitted constant or changing the original target. At
unit indices, `N_UV(0)=b^2/4+h*(h+3*b-c)/12`.

Apply the complete fixed-original-numerator mass derivative to both pieces.
This includes d_a N_UV, the outside factors 1/(h*b*c), the shell Jacobian,
and moving support. S is fixed and is not differentiated. Changing S must
cancel between the two pieces; this provides a useful numerical cross-check.

## Exact isolated q1-parallel-q2 subtraction

Keep r1,r2,r3 positive, q3 separated from the collinear pair, and the triangle
parameters away from their singular boundary for this local derivation. Let

```
t=(1-z12)/2, L=4*r1*r2, h=L*t,
z=z13, b=2*r1*r3*(1-z),
c(t,w)=c0+c1*t+sigma(w)*sqrt(t*(1-t)),
c0=2*r2*r3*(1-z), c1=4*r2*r3*z,
sigma(w)=-4*r2*r3*sqrt(1-z^2)*w,
d_h = partial_a h|0 = -1+r2/r1,
d_b = partial_a b|0 = -1+r3/r1,
F_a=x0*x1*d_h+x0*x2*d_b.
```

At a=0 define, including any common parameter prefactor C,

```
N=b^2/4+h*x1*(h+b-c)/4+h*x2*b/2,
N_a=b*(d_b-1)/2
 +(d_h-1)*[x1*(h+b-c)/4+x2*b/2]
 +h*[x1*(d_h+d_b)/4+x2*d_b/2],
G(h,c)=C*N/(b*c*F^nu),
H(h,c)=C/(b*c*F^nu)
       * [N_a-N*(d_b/b+nu*F_a/F+1/(2*r1^2))].
```

The massless raised bulk, before angular measure, is exactly

```
B(t,w) = -d_h*G(h,c)/h^2 + H(h,c)/h,
f(t,w)=t^2*B(t,w) = -d_h/L^2*G(L*t,c(t,w))+t/L*H(L*t,c(t,w)).
```

The `-K/(2*r1^2)` shell-Jacobian term is present in H. It has not been
dropped because a shell relation was used. The upper Fermi-surface term is
still separate and has the previously derived negative sign.

The azimuth weight is even in w. Pair the evaluations at w and -w before
subtracting: `f_even(t,w)=(f(t,w)+f(t,-w))/2`. This removes odd powers of
sqrt(t); simply applying an integer Taylor subtraction to a single w sample
would leave an unremoved half-power singularity.

The first two paired coefficients are explicit:

```
f0 = -d_h/L^2 * G(0,c0),
f1(w) = -d_h/L^2 * [L*G_h(0,c0)+c1*G_c(0,c0)
                              +sigma(w)^2/2*G_cc(0,c0)]
         +H(0,c0)/L.
```

The normalized azimuth average uses
`<w^2>=1/(D-2)`, so the analytic counterterm replaces sigma(w)^2 by
`16*r2^2*r3^2*(1-z^2)/(D-2)`. No azimuth quadrature is needed for f0 or the
averaged f1. For implementation, at h=0,

```
N0=b^2/4, F0=x0*x2*b+x1*x2*c0,
N_h0=x1*(b-c0)/4+x2*b/2,
G_h0=C/(b*c0*F0^nu)*(N_h0-nu*N0*x0*x1/F0),
G_c0=-G0*(1/c0+nu*x1*x2/F0),
G_cc0=G0*[2/c0^2+2*nu*x1*x2/(c0*F0)
                       +nu*(nu+1)*(x1*x2)^2/F0^2].
```

Set alpha=(D-2)/2. The exact normalized polar measure is
`t^(alpha-1)*(1-t)^(alpha-1)/B(alpha,alpha) dt`. The continued integral is

```
1/B(alpha,alpha) * integral_0^1 dt t^(alpha-3)*(1-t)^(alpha-1)
                   * [<f_even(t,w)>-f0-t*<f1(w)>]
 + f0*B(alpha-2,alpha)/B(alpha,alpha)
 + <f1>*B(alpha-1,alpha)/B(alpha,alpha).
```

The bracket is O(t^2), so its isolated lower endpoint is integrable for
Re D>2. The exact normalized counterterms simplify before numerical use:

```
R0=(2*alpha-1)*(2*alpha-2)/[(alpha-1)*(alpha-2)]
   =2*(2*alpha-1)/(alpha-2),
R1=(2*alpha-1)/(alpha-1).
```

At D=4-2*epsilon these are
`R0=-2*(1-2*epsilon)/(1+epsilon)` and `R1=-1/epsilon+2`.
The apparent alpha=1 ratio in R0 is removable and must be cancelled exactly.
The leading inverse-square pole has a D=6 pole; its next coefficient has the
D=4 collinear pole. Counterterm coefficients retain their full D and analytic
index dependence before the Laurent series is taken.

For a local cutoff t<delta, replace complete beta functions by incomplete
ones and apply the same cutoff to the subtractions. This avoids subtracting
through an unrelated t=1 singular patch. The identity with the full interval
is useful only after overlapping endpoints have been separately resolved.
The upper Fermi-surface kernel has one fewer inverse h power: subtract its
constant term with the same beta-function procedure, retaining its own sign.

## Resolving the overlapping angular zeros

Use stereographic polar radii u=tan(theta12/2), v=tan(theta13/2), and
y=(1-w)/2. Then all angular denominators become positive rational functions:

```
h12=4*r1*r2*u^2/(1+u^2),
h13=4*r1*r3*v^2/(1+v^2),
h23=4*r2*r3*Q/[(1+u^2)*(1+v^2)],
Q=(u-v)^2+4*u*v*y.
```

Split u,v into finite and reciprocal charts so each integration variable lies
in [0,1]. In a finite chart with u>=v, put v=u*(1-z), giving
`Q=u^2*[z^2+4*(1-z)*y]`. The reverse ordering is analogous. Near the pair
collision z=y=0, set y=s^2 and split z>=s versus s>=z:

```
s=z*k: Q=u^2*z^2*[1+4*(1-z)*k^2],
z=s*k: Q=u^2*s^2*[k^2+4*(1-s*k)].
```

Both bracketed factors are bounded below by one for s,z,k in [0,1]. The
z=1 boundary (v=0) is a different collinear face and must first be isolated
by an interval split. Reciprocal charts have the same Q polynomial after
clearing positive denominators. Mixed finite/reciprocal charts have a pair
collision only at their matching boundary; a translation resolves that face.

These substitutions turn the simultaneous collinear zeros into explicit
monomials. Combine them with the largest-radius sectors and positive-simplex
parameter sectors from the existing derivation. The remaining positive F
polynomial has to be sector-resolved too; local angular subtraction alone does
not settle its vanishing parameter faces or soft radial overlaps.

On a final sector, collect each integrand as
`product z_j^(a_j+b_j*epsilon+sum c_jk*lambda_k) * R(z;epsilon,lambda)`
with R smooth and nonzero denominators on the closed cube. Apply the exact
one-variable Taylor/plus-distribution formula to every exponent whose real
part is <=-1; integrate each subtracted monomial analytically. A forest/product
of these operators accounts for overlaps without double counting. All UV,
collinear, radial and Fermi-surface terms must carry one consistent regulator.

## Numerical implementation and remaining gates

### All three two-cut sectors reduce to five single-scale virtual coefficients

This supplies the missing numerator/routing maps, but not their evaluated
coefficients. Use the native Minkowski family

```
D0=K^2, D1=L^2, D2=(K-p)^2, D3=(L-p)^2,
D4=(K-q)^2, D5=(K-L)^2,
p^2=-1, q^2=0, p.q=-1/2.
```

The optional seventh factor `(L-q)^2` is a numerator completion with index
zero in every physical target. For pair [1,5], take K=P1,L=P3,q=P2,p=P4.
For pair [1,7], take K=P4,L=P3,q=P2,p=P1. The original supplemental numerator
then gives respectively

```
NA0=[(D0-D4)^2-(D0+D1-D5)]/4,
NB0=[1+(D1-1-D3)*(D0-D4)]/4.
```

Let s=q_E^2, keep h=p_E^2=1 and (q-p)_E^2=0. At s=0 the required momentum
variation is `dq/ds=p-q`; in native Minkowski products `d_s q^2=-1` and
`d_s(p.q)=-1/2`. It follows exactly that `d_s D4=D2-D4` and

```
NA_prime=(D0-D4)*(D4-D2-1)/2-(D0+D1-D5)/4,
NB_prime=1/2+(D1-1-D3)*(D4-D2-1)/4.
```

For A or B the native virtual C0 integral has numerator N0 over
`D0*D1*D2*D3*D4*D5`. Its off-shell derivative C1 has numerator
`D4*N_prime+N0*(D4-D2)` over
`D0*D1*D2*D3*D4^2*D5`. The latter expression already contains the sign from
differentiating D4. Both C0 and C1 convert to the unscaled Euclidean values
with `(4*pi)^(-D)`; do not apply another seven-denominator sign to C1.
These are the two jets in `I_E(h,s)=h^(D-4)*(C0+C1*s/h+...)`.
This jet notation alone does not exclude nonanalytic branches: derive it in
the common analytic-index domain before taking the original massless limit.

The pair [5,7] has the raised slot 1 uncut. Set
`q1=P2-P4`, `q3=P2-P1`, `p=q1-q3`, `q=-q3`,
`K=P1=P2-q3`, `L=P1-P3`. It has the same six-denominator family, with D4
the original charged slot 1 and hence squared. Its numerator is

```
NC=[(D0+D4)^2+(D0-D1+D5)*(D4+D2)]/4.
```

Here the Euclidean coefficient is **minus** `(4*pi)^(-D)` times the native
seven-power integral, because this is an actual raised Euclidean denominator,
not the derivative combination above. The outside neutral factor is 1/h for
every pair. Thus the complete two-cut part needs C0A,C1A,C0B,C1B,CC, with a
single explicit native-family normalization.

The compact integrations are exact meromorphic Beta moments. Define

```
alpha=(D-2)/2,
K_b = B(alpha+b,alpha)/B(alpha,alpha),
b=D-5, n=2*D-7,
A_d = area(S^(d-1))/(2*pi)^d, d=D-1.
```

For either pair containing the raised cut, including its positive upper
surface, the contribution is

```
(A_d^2/4) * 4^(b-1) * mu^(2*n-2)
 * [2*C0*K_b*(n-1)/(n*(n-2))
    +(b*C0+C1)*K_(b-1)/(n-1)^2
    -b*C0*K_(b-1)/(n*(n-2))].
```

This follows directly from
`W1*W2*[C0*h^b/(2*r1^2)-b*(-1+r2/r1)*C0*h^(b-1)+C1*h^(b-1)]`
plus `+A_d*mu^(d-3)/4 * integral W2*C0*h^b|r1=mu`.
For pair [5,7], the virtual raised line instead gives

```
(A_d^2/4) * 4^(D-6) * mu^(4*D-16)
  * K_(D-6)/(2*D-8)^2 * CC.
```

These formulas must be expanded jointly: the factors near n=1 and the angular
Beta poles can require higher epsilon orders of the virtual coefficients than
the requested final finite part. A fixed-D numerical sample of the virtual
coefficients is useful evidence, but it cannot supply those Laurent orders.

The concrete code path can be a validation-only integrand generator:

1. Generate N, N_a, F and F_a exactly from the fixed numerator and routing;
   attach common neutral-line indices to every cut amplitude.
2. Split the triangle into the explicit UV term and finite Schwinger remainder.
   Apply the polar subtraction above to each, including the upper surface.
3. Resolve remaining positive-polynomial angular/radial/simplex overlaps,
   retain the exact monomial exponents, and emit finite coefficient integrands.
4. Assemble all three two-cut sectors, this three-cut sector and any regulated
   vacuum/one-cut continuation terms before removing common auxiliary indices.
5. Only then expand at D=4-2*epsilon and integrate finite coefficients, with
   precision, node order, angular/parameter sector partition and subtraction
   scale S varied independently. Retain real/imaginary branch conventions.

This gives explicit local counterterms rather than a request to integrate an
unregulated singular kernel. It does not yet provide the required epsilon expansions of all five single-scale
two-cut tensor/jet coefficient integrals, a complete sector list, or a proved common
thermal mass-derivative continuation, or an empirical precision estimate.
Those are substantive remaining tasks. No supplied coefficient, AMF prediction,
or target-dependent fit enters this construction, and no generated reference
should be claimed until the complete refined calculation succeeds.

## Actual bounded independent virtual check

The [partial upstream report](../reports/validation/2026-10-09-finite-density-native-assembly/upstream-prism-virtual-attempt-2/README.md)
records one completed 18-digit-request profile and a timed-out independent
26-digit-request profile. At epsilon=1/5, exact combination reduction cancels
all but three standard massless virtual masters. Independently evaluating
their Gamma formulas agrees with all five saved virtual combinations within
1e-12. This is neither a full prism amplitude nor a Laurent coefficient
reference. The missing common continuation, analytic-index cancellation and
overlapping endpoint subtractions remain necessary.
