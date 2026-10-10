# Polynomial singleton C1 normalized-boost identity

This is a validation-only source proposal. It changes neither physical
finite-density admission nor the ordinary-completion domain. It is a new
guarded source derivation; it must not be described as replay of an original
polynomial source at a forbidden reciprocal-energy seed.

## Narrow domain

* There is exactly one occupied loop, with future momentum q, positive chemical
  potential mu and zero physical mass. Its required-cut index is exactly one.
* Every energy-completion index is exactly zero. Other completion factors are
  Lorentz scalars with the original polynomial index bounds. Physical virtual
  factors and Gram numerators are Lorentz scalars; their integer indices may
  range over the existing admitted domain.
* The lower occupation theta(E) has index zero. The upper occupation theta(mu−E)
  is split into index s=0 and s>=1 branches. Positive lower indices are already
  zero under the existing bound joint dimensional-origin prescription and do
  not need a new row.
* Retain the exact original family, masses, deformation, source conventions,
  and the sealed finite-label origin prescription. The proof is pointwise in
  finite integer labels, with a sufficiently high-dimensional neighborhood per
  finite row. It does not assert a common dimension for all unbounded indices.

## Derivation without reciprocal integrals

Write q=E(1,n), |n|=1 on the future massless shell. At fixed n, give every
virtual loop momentum the same infinitesimal Lorentz boost in direction n.
On the occupied shell this boost acts as E*d/dE, leaves n unchanged, and is
therefore well-defined away from E=0. All pairwise Gram products, scalar masses,
and scalar deformations remain invariant. Virtual Lorentz-invariant measures
have zero boost divergence. Polynomial Gram numerators remain invariant too.

The occupied on-shell measure is proportional to E^(D−3) dE dOmega. Its radial
Euler divergence is D−2. Thus the complete virtual divergence plus occupied
radial integration by parts acts only on the upper and lower occupations.
No energy numerator is differentiated under the stated fixed-zero guards.

One may make the origin step explicit with a smooth cutoff chi(E/delta), equal
to zero near zero and one above 2*delta. The extra origin term has the form
E*chi'(E/delta)/delta and support E approximately delta. Existing finite-label
high-D origin estimates make its integral vanish as delta tends to zero, at
fixed positive eta in the stated regulated prescription. UV meromorphic
continuation of virtual integrals precedes this high-D estimate as in the
bound singleton origin owner. The usual real-Fermi finite-T tail gives no
infinity boundary; its T-to-zero limit yields the compact upper endpoint.
This uses a shell tangent field and does not define E^-p times an apex delta.

As an algebraic check only, the full-space correlated boost divided by E has
divergence D−2+q²/E², and its energy action is delta E_i=E_i−(q.P_i)/E.
On C1(q²), the q² terms vanish away from the regulated origin. The shell proof
above justifies the resulting polynomial identity without admitting the
intermediate reciprocal terms as integrals.

## Exact source rows

Let I_s denote the same full integral with upper occupation index s and all
other indices unchanged. The convention is H0(x)=theta(x) and
Hs(x)=(-1)^(s−1) delta^(s−1)(x)/(s−1)! for s>=1.

For the upper bulk branch s=0:

    (D−2) I_0 − mu I_1 = 0.

The raw upper derivative is −E H1(mu−E). The existing multiplication identity
E H1(mu−E)=mu H1(mu−E) gives the displayed polynomial row. The lower derivative
is E H1(E), zero after the same regulated origin limit (and also the simple
distribution multiplication identity away from the apex).

For the positive upper branch s>=1:

    (D−2−s) I_s + s*mu I_(s+1) = 0.

Here the upper derivative is +s E H_(s+1)(mu−E), and
E H_(s+1)(mu−E)=mu H_(s+1)(mu−E)−H_s(mu−E). No new occupation convention or
theta-at-zero prescription is used.

The rows have no cut shift, no virtual-factor shift, and no energy-completion
shift; only the upper index can increase by one. The original role validity,
guarded source instantiation, exact replay, condition extraction and common
native integral order remain authoritative. Solving these rows can introduce
nonzero conditions such as mu or D−2−s, which must be retained normally.

For search presentation, an exact base shift s=H−1 gives rows that explicitly
lower the current upper index:

    H=1:   (D−2) I_(H−1) − mu I_H = 0,
    H>=2:  (D−1−H) I_(H−1) + (H−1)*mu I_H = 0.

Their upper-index guards are fixed one and at least two, respectively. Every
other guard stays unchanged. This recentering shifts the complete coefficient
and domain; it is not an assertion that the bulk identity holds at every H.

## Excluded cases

At raised cut index n>1, q² Cn=C(n−1), so the reciprocal contact terms do not
vanish. At nonzero energy-completion indices, differentiated numerators can
leave (q.P_i)/E terms. At multiple occupied loops, other occupations likewise
produce mixed energy ratios. None of these cases is covered by this source.
Positive cut mass also leaves mass²/E² terms. Extending any such case needs
its own polynomial cancellation or the separate reciprocal continuation.

The source is valid at finite eta and is not an endpoint-zero shortcut.
Whether adding it improves a bounded native search remains an experiment.
