# Sufficient massless occupied-channel continuation

This derives a candidate common-contour certificate for the uncut amplitude at
strictly positive auxiliary mass. It does not by itself admit massless shell
products, justify the physical eta=0 limit, close a weighted system or establish
numerical acceptance. Production still requires positive occupied shell masses.
The [companion massless distribution construction](finite-density-massless-distributions.md)
specifies an origin prescription; its integration into the complete regulated
weighted measure and an endpoint proof are still required before changing that
restriction.

## Incidence determines all relevant channels

Let G be the admitted vacuum graph, with directed incidence B and a connected
independent occupied cut set C. Remove those edges to obtain the connected
amputated graph H=(V,E\C). Each removed edge carries a physical momentum Q_c and
an oriented future occupied momentum q_c=s_c Q_c, where s_c is the sign of its
nonzero chemical shift. Thus Q_c=s_c q_c. The original graph's exact incidence
and routing identity B R=0 is required; routing rank alone supplies no channel
evidence.

For a spanning two-forest F of H, let S be either component's vertex set. Its
external momentum, up to an irrelevant global sign, is

```text
P_F = sum_c kappa_c(S) q_c,
kappa_c(S) = s_c * [1(head(c) in S)-1(tail(c) in S)].
```

These coefficients are exact integers in {-1,0,1}. They come from cut incidence,
not a guessed momentum picture. Every two-forest with this vertex partition has
the same channel. A partition occurs if and only if both H[S] and H[V\S] are
connected, counting a single isolated vertex as a connected component.

A concrete finite certificate algorithm is therefore:

1. Obtain the native admitted cut sets and their future orientations.
2. Enumerate nonempty proper vertex subsets S, fixing vertex zero in S to remove
   duplicate complements, within an explicit partition budget.
3. Retain partitions whose two induced uncut subgraphs are connected.
4. Compute kappa exactly from the supplied directed cut edges.
5. Admit this sufficient certificate only if every channel is zero, has one
   nonzero coefficient, or has exactly two opposite nonzero coefficients.

For future massless q_i, the admitted channel has respectively square 0, 0, or

```text
(q_i-q_j)_M² = -2*q_i.q_j <= 0.
```

Unequal positive chemical potentials change radial support but do not change this
inequality. A channel with two equal signs or more occupied momenta receives an
explicit unproved result. This does not prove it is pinched; another contour
argument would be required.

## Positive auxiliary mass fixes one common continuation

All uncut physical quadratic factors receive the same positive Euclidean squared
mass eta, with any nonnegative physical uncut squared masses retained. In native
Minkowski convention these factors are `p_e²-m_e²-eta+i0`. With positive Schwinger
parameters alpha_e, the first and second Symanzik polynomials of H are

```text
U = sum_T product_(e not in T) alpha_e,
F = U * sum_e alpha_e*(m_e²+eta)
    - sum_(spanning two-forests F2) P_F2,M²
        * product_(e not in F2) alpha_e.
```

The first sum runs over spanning trees. For a connected graph U>0 in the
parameter interior. For projective parameters sum_e alpha_e=1, the certified
channels give F>=eta*U>0. Without projective normalization the bound is
`F>=eta*U*sum_e alpha_e`. A tree uncut graph has U=1 and is included.

To define the common energy continuation, keep all occupied spatial momenta and
physical shell energies fixed and interpolate **all external energies of the
uncut amplitude together** by q_i0 -> lambda*q_i0, 0<=lambda<=1. The occupation
and shell measures remain fixed. For any admitted channel,

```text
P_F(lambda)_M² = lambda²*P_F0²-|P_F,spatial|²
               <= P_F(1)_M² <= 0.
```

At lambda=0 the uncut amplitude has real Euclidean external energies and admits
ordinary Wick rotation. The positive complete F polynomial then defines its
analytic continuation to lambda=1 without an interior physical threshold. This
continues the full amputated amplitude, including its numerator; it assigns no
independent principal values to its denominators. For complex eta with Re eta>0,
first make this continuation at real eta>0, then analytically continue the
already Wick-rotated parameter representation. Its real part stays positive.
This does not assert a pole-free Minkowski-to-Euclidean homotopy for arbitrary
complex eta at every intermediate contour.

Ultraviolet parameter boundaries require dimensional/analytic regulation before
this argument is used as an integral identity. Polynomial numerators and finite
raised powers change powers and parameter weights, but do not create a new F
polynomial. Zero-index or pinched sectors must inherit the same regulated
parameter representation or receive a separate terminal proof.

## Independent physical-mass derivatives at fixed eta

The common continuation also has a useful one-sided mass neighborhood before the
massless limit. Give the occupied lines independent positive squared masses a_i
and let a_max be their maximum. Future massive vectors obey

```text
(q_i-q_j)_M² <= (sqrt(a_i)-sqrt(a_j))² <= a_max,
q_i² = a_i <= a_max.
```

For a connected graph the unsigned two-forest weights satisfy

```text
sum_F2 product_(e not in F2) alpha_e <= U*sum_e alpha_e.
```

One proof chooses a connecting edge e for each forest and forms a spanning tree
T=F2+e. The forest monomial is alpha_e times that tree monomial. Distinct forests
mapped to one tree are obtained by removing distinct tree edges, so their sum is
bounded by the tree monomial times sum_e alpha_e. Summing over trees gives the
inequality.

It follows that `F >= (eta-a_max)*U*sum_e alpha_e` for real eta>a_max, throughout
the same lambda continuation. Thus every fixed eta>0 has an independent positive
cut-mass neighborhood tending to zero, without changing the virtual branch.
The original momentum numerator and chemical potentials are fixed while taking
physical-mass derivatives. Mass-dependent shell substitutions and moving support
must be differentiated. This neighborhood alone does not prove that the massless
limit commutes with eta->0: it shrinks as eta approaches zero.

## Exact graph evidence

`tools/finite_density/massless_channel_certificates.py` is a read-only mathematical
checker with respect to production code. It reads the supplied graph definitions
and already established cut-set certificates, computes exact incidence channels,
and saves all vertex partitions, coefficients and failures in
[`massless-channels.json`](../reports/validation/2026-10-09-finite-density-native-assembly/massless-channels.json).

| Family | Nonempty occupied cut sets satisfying this sufficient certificate |
| --- | --- |
| Seven-edge chain | 3/3 |
| Eight-edge five-vertex graph | 12/29 |
| Nine-edge triangular prism | 7/7 |

For example, the eight-edge graph's cut `[0,6]` has connected two-forest partition
`{0,2} | {1,3,4}`, with channel coefficients `[-1,-1]` in the future cut basis.
The channel is `-(q_0+q_6)` and can be timelike. The certificate therefore does
not cover the full amplitude of this mandatory family. No result for that family
may omit its other 17 contributing cut sets.

The checker also independently enumerates spanning two-forests to verify every
retained vertex partition, rechecks exact incidence/routing and charge
conservation, and passes all 316 individual edge-orientation reversals over the
39 cut sets.

The native `finite_density::massless_contour::massless_channel_evidence` query
uses the admitted input and occupied-cut owners, native undirected graph
connectivity and an explicit partition budget. It requires precisely all uncut
physical quadratics to be declared shifted by the same unit eta, zero occupied
masses and nonnegative uncut masses. Its typed evidence retains cut orientation,
chemical magnitudes, mass assignments, routing/determinant, every channel and
explicit unsupported witnesses. The deformation descriptor passed to this query
must describe the actual flow when it is eventually integrated into evaluation.
A failed channel is reported as unproved, not as a demonstrated pinch.

Both implementations supply graph/channel evidence only. The native query sets
`endpoint_admitted=false` and `numerical_evaluation_admitted=false`; the numerical
assembly's positive-shell restriction is unchanged. Shell-distribution
admission, common UV/IR regulation, closure and numerical accuracy remain
separate obligations.

All four native certificate unit tests passed in the shared 31-test finite-density
unit run, covering E7/prism sectors, the E8 witness, bounded/partial-deformation
rejection, orientation invariance and rejection of a positive occupied mass by
this specifically massless query. The actual output is retained in
[`normal-unit-resources.log`](../reports/validation/2026-10-09-finite-density-native-assembly/normal-unit-resources.log),
with source identity in `normal-source-hashes.json`. These checks validate the
algebraic query; they do not change either numerical-admission flag.

## Remaining order-of-limits obligations

For fixed eta>0 and positive chemical potentials, a massless shell and its finite
raised/occupation derivatives can be defined first at sufficiently large Re D,
with separate analytic indices controlling uncut UV boundaries, then continued
meromorphically. The companion compact-moment/lower-endpoint proof must specify
this distribution extension and exclude unsupported coincident endpoints at
mu=0. The [explicit zero-jet construction](finite-density-massless-distributions.md)
gives the sufficient finite-index domain Re D>2*n+ell for shell index n and
positive lower-contact index ell. Its jointly continued zero is an origin
distribution prescription. Real empty-support arguments used for positive mass
do not establish that massless contact identity.

The certificate above is intentionally restricted to eta>0. At eta=0, U or F can
vanish at soft or collinear boundaries. Before numerical physical endpoint
projection is admitted, each family needs a compatible convergence/analytic
regulator construction and a complete-amplitude endpoint argument. Terms
scaleless only at eta=0, dimension-dependent powers, logarithms and cancellations
between cut sectors cannot be decided by setting eta to zero in separate
integrands or in a massive Laurent series. The intended order is to establish
the common regulated amplitudes, define the massless fixed-eta distributions,
continue the complete cut sum to its physical endpoint with symbolic dimension,
and only then extract Laurent coefficients. Every exchange of these operations
requires its own convergence or meromorphic identity. Weighted source closure,
integrated boundaries and stable native numerical predictions remain mandatory.

In particular, a massless on-shell single-cut virtual amplitude may have no
ordinary simultaneous UV/IR convergence strip and may develop nonanalytic eta
powers. Positivity for eta>0 is insufficient grounds to discard such branches
in a shared endpoint projector. Each omitted branch needs either convergence
inequalities for the original target and its required raised derivatives, or a
scaleless-region identity established in dimensional/analytic regularization
after the common amplitude continuation. The independent E7 reference's compact
convergence band does not by itself establish that criterion for every generated
master or separately regulated cut sector.

The raised triangular-prism three-cut sector gives a concrete obstruction to
using a single unregulated real-D convergence strip. Keep q3 generic and take
q1 parallel to q2 with unequal positive radii r1,r2, away from every Fermi
surface. For Feynman parameters inside the triangle simplex, h12 is proportional
to the squared relative angle, F stays positive, and the supplemental numerator
at h12=0 is h13²/4, which is nonzero. Differentiating the independent q1 mass
gives h12_a=-1+r2/r1, so the raised bulk contains a nonzero h12^(-2) term.
The angular measure is proportional to theta^(D-3)dtheta: local absolute
convergence therefore requires Re D>6. Even the scalar virtual triangle requires
Re D<6 for ultraviolet convergence; an absolute bound on its linear numerator
is stricter. These conditions cannot hold together. This is an obstruction to
a per-sector dominated-limit argument, not a proof that the dimensionally
regulated complete amplitude fails to exist. Independent analytic propagator
indices or a proved complete-sector subtraction/cancellation must precede the
endpoint and mass-derivative limit. The exact kernels and missing-reference
status are recorded in [the reference diagnosis](finite-density-missing-references.md).

## Future single-cut endpoint proof obligation

A possible generic zero criterion concerns exactly one future massless cut
momentum q, at least one remaining virtual loop, zero masses on every uncut
physical line, and a medium vector appearing only in polynomial numerators.
After virtual integration there is only one external momentum: Lorentz tensor
coefficients can depend on q², and contractions with the medium multiply them
by polynomial powers of u.q. The compact chemical-potential scale therefore
cannot by itself supply a scale inside these virtual tensor coefficients.
This observation is **not** a current admission or a proof that the thermal
contribution vanishes. Without a virtual loop, the one-cut polynomial terminal
is a nonzero counterexample to an unqualified single-cut zero claim.

One prospective construction starts with the off-shell virtual tensor
amplitude in a convergent domain, with independent analytic indices on virtual
lines. Its tensor coefficients form a meromorphic homogeneous family in q².
Continue that family to an open parameter domain in which every coefficient
vanishes with more derivatives at q²=0 than the largest finite cut-derivative
order. Only there restrict it to the massless shell and its derivatives. The
resulting zero family may then be jointly continued in dimension and analytic
indices. Directly multiplying delta derivatives by a singular boundary value
such as (-q²-i0)^lambda, or merely setting q²=0 in a scale-counting argument,
does not implement this construction. The high-dimensional origin prescription
is also required when the compact integration reaches q=0.

An alternative starts from the common positive-eta virtual contour. If all
virtual masses are the same auxiliary eta and q²=0, covariance and scaling
suggest that each tensor coefficient, including any fixed finite number of
q² derivatives, has an eta power hD/2 minus an integer, where h>0 is the number
of virtual loops. Polynomial compact moments carry the remaining medium scale.
Independent ultraviolet index continuation could then provide an open domain
with all these eta powers positive, followed by a vanishing eta->0 limit and
joint meromorphic continuation. A proof must establish this tensor scaling for
the complete retained numerator and all shell/upper/lower distributions; it
cannot silently replace a raised original numerator by its shell value.

Either construction still needs equality with the original regulated thermal
product, its common contour and its prescribed order of line regulators,
thermal limit, cut sum and eta limit. A vanishing meromorphic family chosen by
one order of restriction is not sufficient evidence for that equality.
Overlapping pole contacts and exceptional homogeneity values must be treated
before specializing dimension. No single-cut physical-zero owner or massless
flow endpoint admission follows from this note.
