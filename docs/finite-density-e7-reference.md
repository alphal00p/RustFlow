# Independent reference construction for the seven-edge four-loop family

This is a validation-only derivation for the explicit
`chain_of_three_parallel_pairs` input. It is not a production evaluator, boundary
condition or source of reduction rules. No supplied numerical answer enters the
derivation. The raised compact kernel, neutral tensor/routing contractions,
dimensions and reference-normalization factor have received an independent
algebraic audit. Independent native reference generation and its precision checks
passed; no native AMF prediction or oracle comparison has been made. Do not mark
either target numerically accepted from reference generation alone.

Write D=4-2 epsilon, d=D-1, with common positive chemical potential mu. The exact
edge routing is

```text
P1, P2, P3, P4, P1-P4, P2-P4, P3-P4,
```

and only P1 and P1-P4 carry chemical charge. All physical masses are zero. The
first target is I37, with powers `[1,1,1,2,1,1,1]` and numerator `(P1.P2)^2`.
The supplemental target has powers `[2,1,1,1,1,1,1]` and original numerator
`(P1.P2)^2 + (P1.P3)*(P2.P4)`. Measures here are raw Euclidean loop measures
d^D P/(2 pi)^D, before MSbar or reference normalization.

## Neutral virtual integrations

The P2 and P3 loops contain no occupation functions. Let

```text
G(D) = (4*pi)^(-D/2)
       * Gamma(2-D/2)*Gamma(D/2-1)^2/Gamma(D-2).
```

Feynman parameter integration and the rotational Gaussian covariance give

```text
integral dk 1/[k²(k-p)²] = G(D)*(p²)^(D/2-2),
integral dk k_mu/[k²(k-p)²] = p_mu/2 * G(D)*(p²)^(D/2-2),
integral dk k_mu*k_nu/[k²(k-p)²]
 = G(D)/(4*(D-1)) * (p²)^(D/2-2)
   * [D*p_mu*p_nu-p²*delta_mu_nu].
```

For example the two tensor coefficients before simplification are

```text
A = (4*pi)^(-D/2)/2 * Gamma(1-D/2)*Gamma(D/2)^2/Gamma(D),
B = (4*pi)^(-D/2) * Gamma(2-D/2)
    * Gamma(D/2-1)*Gamma(D/2+1)/Gamma(D).
```

Their exact identities are B=-D*A and A+B=G/4. A cancelled denominator leaves
a translation-invariant scaleless neutral integral. Such terms vanish under the
same analytic regularization before any finite-density spatial integration.

There is no common unregulated real-D convergence band for every intermediate
tensor bubble and the final compact integral. Introduce independent analytic
indices for the neutral propagators to justify integration and tensor identities
in a convergence domain, then continue the complete meromorphic identities.
Deleting UV-divergent terms as ordinary convergent integrals would not justify
these formulas. The original graph's common energy-contour prescription must be
retained when continuing the virtual factors to occupied momenta.

## Two occupied momenta

Set q1=P1 and q2=P1-P4. In the doubly occupied contribution both are future
massless vectors. Its neutral Euclidean quadratic is
h=2*r1*r2*(1-z), where z is the cosine between their spatial momenta. Define

```text
A_d = 2*pi^(d/2) / [Gamma(d/2)*(2*pi)^d],
K_d(c) = angular_average [(1-z)^(-c)]
       = 2^(d-2-c)*Gamma(d/2)*Gamma((d-1)/2-c)
         / [sqrt(pi)*Gamma(d-1-c)],
n = d-1-c.
```

The occupied measure of each simple line is
`A_d*r^(d-2)*dr/2`; the two Euclidean cut signs multiply to plus. The zero-cut
massless vacuum is scaleless. A single-cut virtual amplitude has one massless
external invariant and is also scaleless. With those statements interpreted as
meromorphic regulated identities, the complete simple charged-pair integral is

```text
J(1,1,c) = A_d² * 2^(-c-2) * mu^(2*n) * K_d(c)/n².
```

The compact beta integrals themselves require Re(n)>0 and
Re((d-1)/2-c)>0. Values outside this domain mean analytic continuation of the
complete regulated expression, not divergent direct spatial integration.

## Raised original charged line and moving support

Keep the first charged squared mass a independent while differentiating. At
fixed spatial radii,

```text
E1 = sqrt(r1²+a),
h(a) = -a + 2*E1*r2 - 2*r1*r2*z,
h_a(0) = -1+r2/r1.
```

Apply `-d/da` to the complete simple-line expression, including its shell energy
and theta support, before taking a to zero. The doubly occupied bulk kernel is

```text
W1*W2 * [h^(-c)/(2*r1²) + c*h^(-c-1)*(-1+r2/r1)].
```

The upper-surface term is positive and contains
`delta(mu-r1)*h^(-c)/(4*r1²)` times the second occupied measure. Radial integration
gives, with the common prefactor
`P=A_d²*2^(-c-3)*mu^(2*n-2)`,

```text
bulk/P = K_d(c)/[n*(n-2)]
       + c*K_d(c+1)*{1/[n*(n-2)]-1/(n-1)²},
surface/P = K_d(c)/n,

J(2,1,c) = P/[n*(n-2)]
           * [(n-1)*K_d(c) + c*K_d(c+1)/(n-1)²].
```

For c=3-D the radial derivative requires D>7/2; the angular derivative requires
D>10/3. The band 7/2<D<4 therefore defines these compact terms directly.
In this band, small-a single-cut derivatives carry the positive power
`a^(D/2-2-c)`, and the vacuum derivative carries `a^(D-3-c)`. They vanish in the
massless limit. The possible collinear boundary region of h(a) has derivative
scaling `a^((d-1)/2-c-1)`, also a positive power in this band. Lower radial
boundary terms vanish by n>2. These conditions are essential to taking the
massless limit after the complete mass derivative. A direct substitution into
a massive Laurent series at epsilon=0 does not establish this reference.

## Original four-loop targets

Using `P1.P4=(P1²+P4²-(P1-P4)²)/2` inside the exact neutral tensor result, the
scaleless contacts cancel and the proposed expressions are

```text
I37 = G(D)²*D/[16*(D-1)] * J(1,1,4-D),

raised_supplement
 = G(D)²/[16*(D-1)]
   * [2*(2*D-3)*J(1,1,4-D) + (3*D-2)*J(2,1,3-D)].
```

Both have mass dimension 4*D-12. The powers of mu in the two supplemental terms
agree. The original numerator is retained during the raised-line construction;
shell substitutions are differentiated with their mass dependence, not frozen.

Gamma recurrence further gives the useful cross-check

```text
J(2,1,3-D)/J(1,1,4-D) = (D-3)*(6*D-17)/(2*D-5),
raised_supplement/I37
 = (18*D³-109*D²+191*D-72)/[D*(2*D-5)].
```

For four loops the conventional MSbar integral divided by the unexpanded factor
`(4*pi)^(-8)*(Lambda_bar/2)^(8*epsilon)` equals the raw expression above times
`(4*pi)^8*(exp(gamma_E)/pi)^(4*epsilon)`. This follows from the original four
loop measures and is independent of the number of occupied cuts.

At mu=1, duplication of the Gamma function puts the normalized I37 into

```text
I37_normalized = r(epsilon)*H(epsilon)/epsilon²,
r(e) = (2-e)/[2*(3-2*e)*(1-2*e)^5*(1-4*e)],
H(e) = exp(4*gamma_E*e)*Gamma(1+e)^2*Gamma(1-e)^5*Gamma(1-3*e)
       / [Gamma(1-2*e)^3*Gamma(1-4*e)].
```

Here `r=1/3+85*e/18+1066*e²/27+O(e³)` and
`log(H)=-pi²*e²+O(e³)`; all Euler-constant linear terms cancel. Multiplying the
supplemental ratio, whose expansion is
`25/3-137*e/9+119*e²/54+O(e³)`, gives the independently derived analytic references

| Target | epsilon^-2 | epsilon^-1 | epsilon^0 |
| --- | --- | --- | --- |
| I37 | 1/3 | 85/18 | 1066/27-pi²/3 |
| Raised occupied supplemental target | 25/9 | 617/18 | 20887/81-25*pi²/9 |

Orders below -2 vanish in these meromorphic expressions. Higher orders are not
claimed by this table. Another agent has independently checked the duplication,
rational expansions and both coefficient triples without using oracle answers.
`tests/finite_density_e7_reference.rs` evaluates both unreduced beta expressions
and the simplified Gamma expression using native arithmetic at separate
precisions, and checks these finite coefficients with independently changed
exact epsilon values. Explicit execution passed both tests, with no tests left
ignored, in 0.03672 seconds and peak resident memory 12332 KiB (prebuilt process,
compilation excluded). The largest 50-to-80-digit relative change was 1.132e-54;
the largest beta/tensor-versus-Gamma-expression difference was 2.292e-84. The
finite-coefficient checks at epsilon=1e-20 and 1e-30 passed, with a largest relative
extrapolant difference of 5.096e-18. Exact analytic expressions and all measured
differences are retained in the
[reference report](../reports/validation/2026-10-09-finite-density-native-assembly/independent-e7-reference/reference.json).

The remaining numerical acceptance gates require native AMF predictions and
their independent precision, boundary and epsilon-grid refinements. Reference
values remain outside the production pipeline and may be compared only after AMF
predictions are saved. This derivation does not supply references for the other
two mandatory graph families.
