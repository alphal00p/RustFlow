# Correlated boost: exact polynomial source assessment

This is a mathematical and existing-source audit. No production source, ordering, admission or zero theorem is changed; native performance/closure controls are separate.

Use native Minkowski vectors P_i, the fixed medium vector U with U²=1, energies E_i=U.P_i and Gram coordinates g_ij=P_i.P_j. For any fixed rational routing vector r, let Q=sum_i r_i P_i and E_Q=U.Q. The correlated vector field is

```
v_i = E_i Q - (Q.P_i) U, for every loop i.
```

It applies the same antisymmetric infinitesimal Lorentz map to every loop, although the generator Q itself depends on the loop variables. Its exact coordinate action and total Cartesian divergence are

```
delta g_ij = 0,
delta E_i = E_i E_Q - Q.P_i = X_i,
sum_i div_(P_i) v_i = (D-1) E_Q.
```

The last identity includes differentiation of Q: for each i the two E_Q terms cancel, and the surviving divergence is `(D-1) r_i E_i`. Freezing Q while computing the divergence would incorrectly give zero. For a future occupied coordinate a, Q=P_a. A tilted cut with rational routing is covered directly by the general Q formula, or by the existing occupied-coordinate transformation. A fixed nonzero external momentum that is not transformed along with the loops is outside this statement; the finite-density vacuum graph input has no such external vector.

Every actual Lorentz-scalar physical factor, including arbitrary routed quadratic denominators, cut shells, mass constants and their fixed scalar eta shifts, is invariant pointwise. This holds before shell integration and for all cut orders. Scalar Gram completion factors are invariant as well. Only medium-energy completions, an explicit medium-dependent numerator, and occupation factors contribute derivatives. No inverse energy is used and there is no new massless-origin prescription.

## Distribution signs and finite labels

For compact loop b, put `X_b=E_b E_Q-Q.P_b`. Upper and lower factors are `h_b=mu_b-E_b` and `l_b=E_b`, so their actions are `-X_b` and `+X_b`. With `H_0=theta`, `H_s=C_s` for s>=1:

| factor index | upper contribution | lower contribution |
|---|---|---|
| 0 | -X_b H_1(h_b) | +X_b H_1(l_b) |
| s>=1 | +s X_b H_(s+1)(h_b) | -s X_b H_(s+1)(l_b) |

These are the existing `WeightedMeasure::ibp` rules. Negative occupation indices remain invalid. Lower terms may be dropped only through the already bound massive empty-support or massless dimensional-origin proof. No theta index is itself a zero.

For an original polynomial numerator N(g,E), add `delta N=sum_i X_i dN/dE_i`. With scalar N=1 and one occupied loop a, the bulk identity reads

```
0 = (D-1) M[E_a] - S[E_a²-g_aa] + L[E_a²-g_aa].
```

M, S and L share the unchanged virtual factors and all other indices; S/L have the corresponding upper/lower occupation raised from zero to one. On a massless simple cut, g_aa*C_1(g_aa)=0, yielding `(D-1)M[E_a]-S[E_a²]+L[E_a²]=0`. Once the existing lower-origin proof applies, the lower term is zero.

Raised cuts require the original g_aa term. For `f_a=g_aa-m_a²`, distribution multiplication gives

```
g_aa C_n(f_a)=m_a² C_n(f_a)+C_(n-1)(f_a),
```

where C_0 is the required-cut zero. Therefore the upper term contains a lower-cut correction. This source never raises a physical cut or virtual denominator; the polynomial coordinate-to-factor map may lower physical indices. Replacing g_aa by m_a² before multiplication would lose a real contribution for n>=2.

As an independent one-shell check at mass zero and N=E^r, lower contacts removed by their stated prescription, the identity is

```
0 = (D-1+r) M_n(E^(r+1)) - r M_(n-1)(E^(r-1))
    - S_n(E^(r+2)) + S_(n-1)(E^r).
```

The exact radial binomial formula verifies this identity, including both lowered-cut terms. For n=1 it gives `S_1(E^(r+2))=(D-1+r)M_1(E^(r+1))`. The r=0 upper surface is nonzero. These references test source signs; they are not virtual analytic reduction formulas.

For an already positive upper index s, write M_(n,s)(r) for the same raw compact integral with E^r and H_s(mu-E). The sign reverses relative to differentiating theta:

```
0 = (d+r) M_(n,s)(r+1) - r M_(n-1,s)(r-1)
    + s M_(n,s+1)(r+2) - s M_(n-1,s+1)(r),
d=D-1, s>=1.
```

After dividing by A_d/2, its exact reference is

```
M_(n,s)(r) = (-1)^(n+s-2) binom(d/2-1,n-1)
             falling(d+r-2n,s-1)/(s-1)!
             * mu^(d+r-2n-s+1).
```

All four terms have the same chemical-potential exponent. The independent check includes mu=1 and3/2 and retains the n-1 terms; C_0 alone is zero. Lower-contact terms use the existing joint high-D continuation, not a pointwise apex assignment at the sampled dimension.

## Relation to the existing source corpus

`lorentz_ibps_partitioned` emits separate fields `v_i=P_j` and `v_i=U`, each with one differentiated loop. It does not emit the correlated vector field or its full sum. `compact_tangent_ibps` emits the individual field `E_i P_j-g_ij U` only when the differentiated loop i is compact. In a singleton sector with virtual loops it thus lacks the virtual components needed to annihilate every Gram coordinate. If every loop is compact, the boost is a finite sum of the existing radial tangent fields, but is still not emitted as one combined source.

The new source is nevertheless in the polynomially shifted original Lorentz span:

```
sum_i [ IBP_i(Q ; E_i W) - IBP_i(U ; (Q.P_i) W) ].
```

Here `IBP_i(v ; F)=integral dP div_(P_i)(v F)`. Expand E_i and Q.P_i in the exact **deformed** non-occupation factor basis and constants. For each monomial, shift the whole original integrand label, substitute those shifted labels into coefficient polynomials/nonzero conditions, and pull back its domain. Eta/mass constants in the affine inverse map must remain. This is a finite exact source presentation, not a new elimination procedure.

Simply multiplying an already differentiated row by E_i or Q.P_i inside an integral is not the same operation. The product rule for the basic Q=P_a case produces the net correction `-E_a W` when written as `sum_i [E_i div(P_a W)-g_ai div(U W)]`; shifting the whole integrand includes this correction automatically. An exact original-row-span check should therefore use shifted integrands, not uncorrected coefficient multiplication.

Use the existing `WeightedMeasure::ibp` with zero actions on all Gram coordinates, energy actions X_i and divergence `(D-1)E_Q`; it will partition only active occupation axes and preserve storage tails. Intersect the original polynomial completion domain and retain the original zero-domain corpus. New source identity must bind the presentation and exact routing; discover and replay all rules fresh for each control.

## Potential benefit and limit

The generated scalar row can contain only a compact energy moment and an upper surface, after proved required-cut/lower-contact reductions. It avoids intermediate raised virtual factors and cut orders which appear while native elimination reconstructs the correlation from separate rows. This is a concrete source-presentation improvement to test on the same targets and native budgets.

There is no closure claim: native physical-degree ordering may still prefer eliminating energy moments in favor of surfaces, and higher occupation indices require the existing multiplication/source combinations. The row is not allowed to override ordering or conditions. Compare the original full corpus against that corpus plus the exact correlated rows, with identical domains, point fixtures and budgets. A lower frontier alone does not prove final derivative closure.
