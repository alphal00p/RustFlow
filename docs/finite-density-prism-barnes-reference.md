# Independent prism three-cut Barnes construction

This is a validation-only derivation, written before any full prism reference
values. It reads the original input and exact numerator, not supplied oracle
coefficients or native predictions. The angular identity below has independent
direct quadrature checks. The complete raised-target Laurent reference through
order zero now passes the independent refinements recorded below. Its
uncertainty is an empirical estimate, not a rigorous enclosure. This reference
does not establish production prism admission or native four-loop acceptance.

## Original raised target and measure

Use the exact triangular-prism routing in
`examples/finite_density/triangular_prism.json`. Its required supplemental target
has powers `[1,2,1,1,1,1,1,1,1]` and the original Euclidean numerator
`g1_2^2+g1_3*g2_4`. Set the three future occupied momenta to
`q1=P2`, `q2=P2-P4`, `q3=P2-P1`, with positive energies `r1,r2,r3`.
Let `h=h12`, `b=h13`, `c=h23`, where
`hij=4 ri rj tij`, `tij=(1-ni.nj)/2`. Each direction is normalized to unit
angular measure. Write

```
p = D-2, alpha=p/2, nu=3-D/2,
Ad = area(S^(D-2))/(2*pi)^(D-1),
Wi = Ad/2 * ri^(D-3) dri,
F = x0*x1*h+x0*x2*b+x1*x2*c,  x0+x1+x2=1.
```

The remaining virtual triangle has mean `x1*P4+x2*P1`. Before differentiating
its occupied mass `a=m1^2`, the numerator and invariant derivatives are

```
N(a) = (b-a)^2/4 + (h-a)/2*[x1*(h+b-c)/2+x2*b],
h_a = -1+r2/r1,  b_a = -1+r3/r1,
N0 = b^2/4+x1*h^2/4+x1*h*b/4-x1*h*c/4+x2*h*b/2,
N_a = b*(b_a-1)/2
    +(h_a-1)*[x1*(h+b-c)/4+x2*b/2]
    +h*[x1*(h_a+b_a)/4+x2*b_a/2].
```

The unraised three-cut term is `-integral W1 W2 W3 K`, with
`K=(4*pi)^(-D/2)*Gamma(nu)*N/(h*b*c*F^nu)` integrated over the simplex.
The original raised target is its negative mass derivative. The bulk kernel
is therefore `K_a-K/(2*r1^2)`, and the moving upper boundary is
`-Ad/4*mu^(D-4)*integral W2 W3 K(r1=mu)`. The upper boundary must be retained.
The original numerator is fixed when taking this complete mass derivative.

The exact definition-only script in `prism-triple-cut-reference/derive_monomials.py`
expands this into 42 monomials: 18 bulk terms proportional to `Gamma(nu)`,
19 bulk terms proportional to `Gamma(nu+1)`, and five upper-surface terms.
The common normalization is `Ad^3/[8*(4*pi)^(D/2)]`. Surface coefficients already
include `-1/2`; their `r1^(D-3) dr1` measure is replaced by `mu^(D-4)`.
Every term has total mass dimension `4D-16`. These are original-integrand
identities, not numerical references.

## Three independent directions

For normalized directions on `S^p`, define
`A(a,b,c)=average[t12^a*t13^b*t23^c]`. Then

```
A(a,b,c) = Gamma(p)^2 Gamma(alpha+a) Gamma(alpha+b) Gamma(alpha+c)
             Gamma(p+a+b+c)
           /[Gamma(alpha)^3 Gamma(p+a+b) Gamma(p+a+c) Gamma(p+b+c)].
```

One derivation diagonalizes each zonal kernel `t^a`. Its degree-l eigenvalue is
`Ka*(-a)_l/(p+a)_l`, where
`Ka=Gamma(alpha+a)*Gamma(p)/(Gamma(p+a)*Gamma(alpha))`.
The harmonic multiplicity is
`(2*l+p-1)/(p-1)*(p-1)_l/l!`. The trace of three kernels is consequently the
Rogers–Dougall very-well-poised 5F4 with upper parameters
`p-1,(p+1)/2,-a,-b,-c` and lower parameters
`(p-1)/2,p+a,p+b,p+c`. Its Gamma evaluation is
[NIST DLMF 16.4.9](https://dlmf.nist.gov/16.4#E9).
The ordinary angular convergence conditions are
`Re(alpha+a)>0`, `Re(alpha+b)>0`, `Re(alpha+c)>0`,
and `Re(p+a+b+c)>0`; the last is the simultaneous three-direction collision
condition and also the Dougall convergence condition. First derive the formula
in an absolutely convergent domain, then continue meromorphically.

The normalization checks are `A(0,0,0)=1`, `A(a,b,0)=Ka*Kb`, and
`A(1,1,1)=(1-1/(p+1)^2)/8` (equal to `1/9` at p=2).
Independent direct integration uses `x=n1.n2`, `y=n1.n3` and the relative
transverse cosine z. The third cosine is
`x*y+sqrt(1-x^2)*sqrt(1-y^2)*z`. The normalized weights are
`(1-x^2)^((p-2)/2)`, the same for y, and
`(1-z^2)^((p-3)/2)`. Orders 24,40,64,96 in saved binary64 diagnostics converge
to the Gamma expression for positive and negative noninteger powers; the
largest final relative difference across the three selected cases is 2.63e-13.
This checks the angular identity, not the full Laurent integral or a rigorous
error bound. Higher-precision Barnes generation remains separate.

## Two-dimensional Barnes formula

For an individual saved monomial let its transfer powers be A,B,C, energy
powers k1,k2,k3, simplex powers X0,X1,X2, and Gamma shift l=0 or1.
Independent positive virtual indices beta_i and transfer regulator delta
provide an initial common convergence domain. Set
`m=sum(beta_i)-D/2+l` and

```
a=A+delta+z,
b=B+delta+w,
c=C+delta-m-z-w.
```

The virtual triangle identity is

```
Gamma(m)*F^(-m) = integral dz dw/(2*pi*i)^2
  Gamma(-z) Gamma(-w) Gamma(m+z+w)
  (x0*x1*h)^z (x0*x2*b)^w (x1*x2*c)^(-m-z-w).
```

This follows by two applications of the inverse Mellin transform of
`(1+t)^(-m)`; see the contour-separation discussion in
[NIST DLMF 5.19(ii)](https://dlmf.nist.gov/5.19#ii).
The simplex integration multiplies the integrand by

```
Gamma(beta0+X0+z+w)
Gamma(beta1+X1-m-w)
Gamma(beta2+X2-m-z)
/[Gamma(sum(beta_i)+sum(X_i)-2*m)*product Gamma(beta_i)].
```

The angular integration supplies `A(a,b,c)`, and the compact radial integrations
at mu=1 supply

```
4^(a+b+c)/[(p+a+b+k1)*(p+a+c+k2)*(p+b+c+k3)].
```

For upper-surface terms omit the first radial denominator. All these formulas
retain the common factor and signed monomial coefficients above. At the
physical indices, set beta_i=1 and delta=0 by analytic continuation. The
identity is proved first where the simplex, angular, radial, and Barnes
integrals all converge. A draft continuation starts at D=7, beta_i=2, delta=3,
where straight separated contours exist for each original monomial.

Straight-contour substitution at D=4-2*epsilon alone is not that continuation.
Every crossed pole must contribute its oriented residue. The exact rational crossing diagnostic encounters repeated pole families.
The symbolic analytic regulator and its finite-part removal below resolve
those collisions. At the initial draft checkpoint, the high-D angular check
alone did not validate this continuation step. The subsequent pole/residue and
uniform-tail checks, complete cut assembly, and numerical refinements now
establish the independent reference recorded below. This validation algorithm
remains separate from production, and no native AMF comparison is claimed.

## Regulator removal and continuation checks

The continuation diagnostic now retains an independent symbolic regulator rho:
`beta_i=1+rho/(7,11,13)_i`, transfer shifts
`delta_ij=rho/(17,19,23)_ij`, and radial shifts
`kappa_i=rho/(29,31,37)_i` at the endpoint. During its initial path these
are interpolated from D=7, beta_i=2, delta_ij=3 and kappa_i=0.
The added radial factors are applied to the complete differentiated kernel;
they are auxiliary analytic indices and are removed before a reference is
reported. They are not differentiated as part of the original mass jet.

For the sample epsilon=1/101, exact rational inspection shows that every
nonconstant Gamma argument stays away from its pole set as rho goes from its
small positive seed to zero. The smallest real contour distance at zero is
larger than 0.286. The surviving two- and one-dimensional integrals are regular
there. All remaining regulator poles are simple and belong to isolated
zero-dimensional residues. Their finite parts use

```
Gamma(-n+c*rho) = (-1)^n/(n!*c*rho)
                 * [1+c*(H_n-EulerGamma)*rho+O(rho^2)],
Gamma(x+c*rho) = Gamma(x)*[1+c*psi(x)*rho+O(rho^2)].
```

Here n is a nonnegative integer and x is not a Gamma pole. These follow from
Gamma recurrence, its [pole residues](https://dlmf.nist.gov/5.2), and the
[integer digamma values](https://dlmf.nist.gov/5.4#E14). For a product with net
simple pole, its finite part is the leading residue times the sum of the
signed logarithmic derivatives, including the rho derivative of `4^(a+b+c)`.
Reciprocal Gamma zeros participate before counting the net pole. No tiny-rho
numerical extrapolation is required. The standalone native-arithmetic draft
checks cancellation of the total 1/rho residue in every original monomial.
It uses only Precision and Symbolica Gamma/digamma, not AMF or feature values.

A two-dimensional Barnes test whose value is a product of two applications of
[Barnes's first lemma](https://dlmf.nist.gov/5.13#E3) checks the crossed-pole
orientation and continuation. A separate test continues the elementary
three-term Mellin identity through negative noninteger m. Direct high-D
original-kernel integration supplies an additional check of the actual prism
numerator, shell Jacobian and moving upper surface.

The diagnostics also exposed and retain a failed contour optimization.
When moving a w contour across a mixed z+w pole, eliminating z leaves the
wrong surviving contour. The corrected procedure moves coordinates
sequentially, eliminates the coordinate being moved, and separately moves
any newly generated one-dimensional contour away from its poles. Old
high-precision values computed on the defective representation are marked
failed diagnostics and are not reference values.

At D=19/2 and rho=1/7, the corrected Barnes result and 1,048,576-point direct
original-kernel diagnostic differ by 4.62e-5 relative, within the direct
integration's observed refinement changes. After exact rho removal, the
35-digit native-arithmetic sample likewise agrees with the independent
unregulated high-D direct check to about 4.8e-5. These are useful checks of the
continuation; they do not replace high-precision refinement of the final
Laurent coefficients. The subsequent complete assembly and Laurent refinements
are recorded below; native AMF comparisons remain pending.

## Uniform regulator removal and complete finite-epsilon samples

The reference covers the supplemental target
`[1,2,1,1,1,1,1,1,1]` with original numerator
`g1_2^2+g1_3*g2_4`, at mu=1. It does not yet cover the other prism target
`g2_3^3`. Its complete sector list is the vacuum, singletons `[1]`, `[5]`,
`[7]`, pairs `[1,5]`, `[1,7]`, `[5,7]`, and triple `[1,5,7]` (zero-based
physical slots). The vacuum is dimensionally scaleless. Each singleton has
only the one on-shell massless external momentum and is zero in the common
analytic-index/dimensional continuation, including the original finite mass
jet; the [singleton germ proof](finite-density-massless-germ.md) specifies the
ordered thermal limits. For this actual routing no uncut row is proportional
to the selected occupied row, so there is no external-only q-squared factor.
The virtual quotient has rank three. Cut `[1]` has P=8 and mass-jet budget
J=1; cuts `[5]` and `[7]` have P=9 and J=0. Their common germ margin is
`3*Re(D)/2-9>0`; origin bounds are D>4 for the raised cut and D>2 for
the other cuts. A regular neighborhood of D=13/2 satisfies all three.
These zeros are not separate real principal values.

The three pair contributions use the five independently native-verified
virtual identities and the compact formulas in
[the original reference construction](finite-density-prism-reference-construction.md).
The pair `[5,7]` contains the raised propagator inside the virtual integral;
the other pairs contain its complete shell mass jet and upper contact.
The fixed original two-virtual-loop targets now also have an independently
reviewed UV-continuation bridge, preserved under
`prism-support-wide-proof/two-virtual-thermal-bridge/`. Its exact audit checks
all 192 active-support queries and the primary/transverse UV charts. The proof
retains independent same-sign line-mass Feynman regulators at fixed real T,
the Fermi pole-strip condition, regulator removal before the thermal limit,
one-sided physical mass jets at positive eta, a generic high-D UV Taylor
continuation, the separate joint-eta Mellin remainder, and compact origin and
angular bounds. This resolves the identified original-target reference
prescription obligation. The production bridge and sealed-permit integration
remain unmerged; this is not new numerical admission.

The triple contribution uses the 42 original-integrand mass-jet monomials
above. No derivative of a massless answer is used. All terms first use the
same raw Euclidean measure. Conversion to the stated reference convention
at mu=1 multiplies the sum by
`(4*pi)^8*(exp(EulerGamma)/pi)^(4*epsilon)`.

For a surviving two-dimensional Barnes leaf the imaginary Gamma arguments
are linear forms in y, v, and y+v. Let their signed numerator-minus-denominator
absolute-slope totals be A, B, S. The exact coercivity constant is

```
gap = min(A+S, B+S, (A+B)/2),
A*abs(y)+B*abs(v)+S*abs(y+v) >= gap*(abs(y)+abs(v)).
```

Stirling's formula therefore bounds the Gamma ratio by a polynomial times
`exp(-pi*gap*(abs(y)+abs(v))/2)`. At epsilon=1/101 every two-dimensional leaf
has gap 2 and every one-dimensional leaf has gap 4. Affine real Gamma
arguments vary over a compact set for small rho, so the polynomial power
bound is uniform. Together with the exact pole-separation check, this
justifies taking rho to zero under every surviving contour integral.
The generation controller repeats both checks at **each** epsilon sample;
the earlier numerical contour-distance bound is not asserted for other
samples.

The independently generated complete finite-epsilon value at epsilon=1/101
changes by `7.1063e-17` relative between 45-digit/step-1/20/cutoff-12 and
55-digit/step-1/30/cutoff-14 Barnes calculations. Both check all 42 regulator
residue cancellations. Exact inputs, arithmetic source, resource records,
hashes, and the observed difference are saved under
`reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference/finite-epsilon/`.
This is empirical finite-epsilon refinement, not a rigorous error enclosure,
not a Laurent reference acceptance, and not a comparison with native AMF.
The complete Laurent generation is a separate saved grid calculation.

The physical analytic-index point is also independent of the auxiliary rho
**direction**. Choose sufficiently high nonresonant D with strict virtual,
radial, and angular bounds. These persist in an open complex neighborhood of
all physical indices; the UV-meromorphic triangle is holomorphic there after
avoiding its Gamma pole hyperplanes. This proves an ordinary multivariate
holomorphic limit before continuing in D. As a separate computational check,
replacing the three regulator denominator triples by `(41,43,47)`,
`(53,59,61)`, `(67,71,73)` changes the complete epsilon=1/101 value by only
`2.985e-34` relative at fixed 55-digit quadrature. The exact alternative
continuation, pole/tail checks, source and resource records are saved under
`regulator-direction-check/`, with comparison in
`regulator-direction-comparison.json`.

At epsilon=1/10000 the original complete sample also refines from
70 digits/step-1/30/cutoff-14 to 80 digits/step-1/40/cutoff-18 by
`1.209e-21` absolute (`2.785e-26` relative). This checks the substantial
cancellation among cut terms close to four dimensions. Its saved
`near-zero-pilot-comparison.json` is a finite-epsilon check only; a
one-sample controller output is not used as a Laurent fit.


## Complete Laurent reference and refinement

The saved [reference](../reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference/reference.json)
contains all eight cut subsets and orders `[-4,0]` for the supplemental raised
polynomial target. The [saved-output validation](../reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference/laurent-reference-validation.json)
passes all 20 raw coefficient refinement comparisons across five complete
12-sample profiles. Epsilon denominators 10000, 500000, 1000000, arithmetic
precision 70/80 digits, and Barnes step/cutoff 1/30,14 versus 1/40,18 are varied
independently. Ten-versus-twelve-node interpolation and imaginary-zero checks
also pass. All60 samples retain exact pole-separation and exponential-tail
checks and all 2520 original-monomial regulator-residue cancellations.

For the comparison normalization stated above, the nonzero fitted
coefficients are

| Order | Coefficient | Empirical absolute uncertainty |
| --- | ---: | ---: |
| -1 | 4.33925042047248508077 | 2.88e-21 |
| 0 | 17.03999268944775288649 | 5.03e-16 |

The fitted coefficients at orders -4,-3,-2 are retained in full precision in
the artifact; they were not set to zero. The uncertainty is the largest
nearby/fine-grid, precision, quadrature, or ten/twelve-node change. It is an
empirical refinement estimate, not a rigorous bound. The initial coarser-grid
finite-coefficient change 1.43e-12 remains a convergence diagnostic; an
additional nearby-grid run was completed to establish the final estimate
without loosening the requested 1e-12 absolute reference criterion.

The raw Euclidean finite coefficient is approximately
`4.32431526349508047544e-8`, with empirical uncertainty `1.66e-24`.
The complete definition, exact continuation inputs, source/executable hashes,
per-sample resources and failed-contour diagnostics are preserved. The other
prism target `g2_3^3` still uses its separate supplied reference. No native AMF
prediction, supplied oracle coefficient, Kira or Wolfram invocation enters
this generator or its required comparison gate.


The final evidence audit also replays both exact regulator and tail analyses
against all 60 archived numerical contours. All 120 regenerated reports equal
the saved reports. The comparator binds each evidence epsilon, source hash,
contour archive hash and actual pre-launch numerical input hash. This audit
changes no numerical sample or fitted coefficient.
