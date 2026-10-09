# Finite-density measure, deformation and limits

This document defines the intended mathematical contract and distinguishes it
from capabilities demonstrated by the implementation. The acceptance plan is
[`finite-density-plan.md`](finite-density-plan.md), with exact inputs and open
reference requirements in
[`finite-density-acceptance.json`](finite-density-acceptance.json). Defining a
measure or compiling an input does not claim its AMF system or boundaries close.

## Physical inputs and original contour

Write D = 4 - 2 epsilon, d = D - 1, and
P_i = (k_i0 + i nu_i, **k**_i), with real integration coordinates k_i. Each
physical edge has routing Q_e = sum_i R_ei P_i, independent squared mass m_e²,
and inverse propagator rho_e = Q_e² + m_e². Independent species potentials enter
nu_i through integer conserved charge cycles. Incidence, edge orientation and
charges remain part of the measure when an index is zero. Scalar numerators use
g_ij = P_i.P_j and u_i = P_i0, with Euclidean bilinear products. An unshifted input
must substitute all offsets explicitly:

```text
g_ij = k_i.k_j + i*nu_i*k_j0 + i*nu_j*k_i0 - nu_i*nu_j
u_i_coordinate = k_i0 + i*nu_i
```

The second line uses `u_i_coordinate` for the medium coordinate, not the chemical
parameter nu_i. This avoids conflating an inverse-propagator numerator with an
unshifted polynomial.

The original integral is the T -> 0+ limit of the regulated thermal contour,
with masses kept independent. The complete contour, its orientation and regulator
limits are retained through cutting and raising. For simple lines the
[cutting-rules paper, sections 2 and 4 and Appendix C](https://arxiv.org/html/1609.04339v2)
provides connected cut complements, a zero-cut vacuum term and regulated energy
integration. Its numerator and independent-chemical-potential extensions require
preserving energy substitutions and charge assignments. Its mass differentiation
argument does not license setting masses to zero before treating infrared limits.

The contour-level principal value is not a collection of independently assigned
principal values for denominators in an on-shell amplitude. Continue the complete
amplitude from real external energies using the common inherited contour before
taking physical on-shell limits. A pinched or unspecified continuation must be
reported explicitly. Merely choosing a value for theta(0) is insufficient.

## Occupation, shell and surface distributions

For one positively oriented charged line let E = sqrt(q² + m²), mu >= 0 and
W(q;m²,mu) = theta(mu-E)/(2E). A simple occupied cut supplies -W times the
continued vacuum amplitude, integrated with d^d q/(2 pi)^d. The same measure is

```text
W(q;m²,mu) = integral_0^infinity de theta(mu-e) delta(e²-q²-m²).
```

This identity fixes shell normalization. In the independent-energy form, a mass
derivative differentiates the shell distribution. In the integrated form it
differentiates E, the occupation and the amplitude's on-shell energy. These are
two representations of one derivative; adding both double-counts the surface.

Direct differentiation at fixed q and mu gives

```text
dW/d(m²) = -theta(mu-E)/(4E³) - delta(mu-E)/(4E²).
```

For a smooth test amplitude A evaluated at q0=iE, a raised original Euclidean
line follows from -d/d(m²) of the complete simple-line result. Its occupied piece
therefore contains d(W A)/d(m²), including the two negative terms above and
W i/(2E) times the derivative of A with respect to its continued energy. Other
explicit mass dependence and input-basis coefficients are differentiated too.
Only after all derivatives are taken may different line masses be equated.

For a smooth spatial vector v and test function F, differentiating the integrated
weighted measure gives the exact distributional source

```text
0 = integral d^d q {
    W [div(v)*F + v.grad(F) - (v.q)/E²*F]
    - delta(mu-E)*(v.q)/(2E²)*F
}.
```

Here the product has compact support and is distributionally differentiated as
a whole. Radial reformulations must retain lower-energy endpoints. This formula
is a source identity, not a proof that a chosen finite reduction basis closes.
[Integrating by parts at finite density](https://arxiv.org/html/2304.05427v2)
explains why differentiated occupation functions produce additional terms and why
the zero-temperature limit cannot simply discard them. Guarded RustRed replay
verifies algebraic consequences of supplied sources; it does not prove those
sources physically valid.

Occupation index zero denotes theta support. A required-shell index at zero has
a different distributional meaning; its vanishing rule must not be applied to an
occupation. Caches must distinguish these semantics and retain support,
normalization, orientation, physical masses and deformation identity.

For mu < m, W vanishes. For fixed m > 0 near mu=m, its integral against a regular
nonzero test function scales as (mu²-m²)^(d/2) times a smooth factor. One physical
mass derivative lowers this threshold exponent by one; repeated derivatives may
be singular or distributional. This elementary estimate follows by q=sqrt(mu²-m²)
r and requires Re(d)>0 before continuation. It does not define a generic
threshold limit for singular test amplitudes. At m=0 the lower endpoint changes
and must be analyzed with symbolic dimension, before Laurent expansion.

## Fixed-shell auxiliary deformation

First express the physical input numerator exactly in completed inverse
propagator and medium coordinates. Negative powers refer to those same slots.
Along auxiliary flow hold the physical-point expansion coefficients and indices
fixed. This convention differs from physical-mass differentiation, which holds
the original momentum polynomial fixed and differentiates mass-dependent basis
coefficients.

For each admissible occupied cut set C and a nonempty placement S of physical
uncut quadratic factors, the proposed deformed term is the same inherited
contour integral with

```text
rho_e(t) = Q_e² + m_e² + t,  e in S,
rho_e(t) = Q_e² + m_e²,      e not in S.
```

All physical occupied shell masses, positive-energy orientations and occupation
factors stay fixed. The native Minkowski convention is
D_e(eta) = q_e² - m_e² - eta; therefore t=eta under Wick rotation, consistently
with `IntegralFamily::deform` in `src/family.rs`. The path of eta is part of the
continued-amplitude prescription; a freely chosen detour around an individual
cut denominator does not establish the required branch.

At a large positive t, fully occupied compact terms whose every uncut denominator
is shifted admit a simple candidate nonsingular domain. If the continued energy
of an uncut line obeys |Q_e0/i| <= B_e on support, then
Re(rho_e(t)) >= m_e² + Re(t) - B_e². Choosing Re(t) larger than every B_e²-m_e²
keeps all these denominators away from zero. For mixed terms the uncut virtual
integrations are not compact; their Landau/contour admission needs additional
analysis. Partial placements can leave soft poles and require separate checks.
This bound alone is neither a universal continuation theorem nor a boundary
value.

For the specified massive two-loop example there is a stronger sufficient
admission. Write the two charged masses as m1, m2 and the neutral mass as M.
The physical values m1=m2=1/2 and M=1 lie in the open positive-mass domain
|m1-m2| < M. Both charged masses also remain below mu=1 in a neighborhood of
this point, so independent derivatives stay away from their support thresholds.

In either one-cut amplitude, the remaining bubble has Feynman-parameter
polynomial

```text
F(x) = x*m_uncut² + (1-x)*M² - x*(1-x)*m_cut²
     = [x*m_uncut-(1-x)*M]²
       + x*(1-x)*[(m_uncut+M)²-m_cut²] > 0,  0 <= x <= 1.
```

The mass inequality gives the strict interior bound and positive masses give
positive endpoints. Deforming either admitted uncut squared mass by t >= 0 adds
x*t or (1-x)*t, preserving positivity. In the fully occupied amplitude the two
future shell momenta satisfy q1.q2 >= m1*m2 in Minkowski signature, hence

```text
rho_neutral = M² - (q1-q2)_M² >= M²-(m1-m2)² > 0.
```

Its fixed-shell deformation adds t and also stays positive. The zero-cut vacuum
has positive Euclidean denominators. These bounds establish a common real
positive-t continuation domain for this selected small graph, including an open
neighborhood for independent physical-mass differentiation. Define it first in
a convergence domain and continue in dimension as needed. The thermal boundary
prescription is inherited; no independent denominator principal values are
introduced. This proves neither generic mixed four-loop admission nor weighted
closure, integrated boundaries or numerical assembly for the example.

The same graph's uncut unit-numerator vacuum contribution is strictly positive
in a real convergence domain such as D=2: all three masses are positive, both
subintegrals and the overall integral converge, and the Euclidean integrand is
positive on a set of nonzero measure. Its meromorphic continuation therefore
cannot be discarded as a scaleless zero.

The zero-cut contribution uses ordinary vacuum AMF. Terms with no admissible
uncut physical factor need a justified compact-measure terminal. No dummy
propagator is introduced. Deforming a shell mass at fixed mu would instead
collapse its occupation support for large positive t; a zero large-mass boundary
would then lose the threshold information. That is not the deformation defined
here.

## Boundaries, target projection and order of limits

An occupied momentum must annihilate every hard-region subspace. Branch-based
region enumeration must nevertheless retain correlated hard virtual directions,
all admitted placements and exact determinants. `factor_region` and HEPKit's
tensor projector must retain the projected hard/soft expressions and coordinate
maps before assigning ordinary and weighted owners. A nonpolynomial soft factor
requires recursive weighted evaluation with a strictly decreasing complexity;
a formal Taylor expansion or replacement by ball moments is not an integrated
boundary value.

Use `RecursiveBoundary` and ordinary AMF/FT for hard factors, with bubble subloop
shortcuts disabled and the tadpole-only recursive terminal policy. Preserve
symbolic epsilon exponents and explicit logarithmic region coefficients in
matching. An unspecified log coefficient is unknown, not zero. Insufficient rank
or expansion depth must remain a failure.

Reconstruct targets by multiplying their full rational weights into the power/log
series before `frobenius::project_limit`. Componentwise physical limits can lose
cancellations. The evaluation path needs a finite-density example with cancelling
singular weights, an uncancelled divergence and insufficient-order failure. Tests
of the ordinary projector alone do not satisfy that requirement.

A valid calculation must establish a common convergence domain, continue the
complete regulated amplitude, remove the auxiliary parameter, take the specified
physical-mass limits with any necessary regions/distributions, and only then
extract the requested Laurent coefficients. Exchanging these operations requires
proof for the admitted family. The ordinary dimensional endpoint rule does not
by itself provide such a proof for occupation-weighted integrals. No general
four-loop admission theorem or successful evaluation is claimed here.

## Whole-amplitude normalization check

For an ordinary vacuum contribution with physical powers a_e, define the native
Minkowski integral with one d^D k/(i pi^(D/2)) measure per loop. Begin with the
standard Wick map k0=i P0, then reverse all native Minkowski momenta. We choose
q_M=(-i*P_E0,-P_Espatial), so the Euclidean crossed pole P_E0=+iE maps to future
q_M0=E. The global Minkowski loop reversal preserves its measure and +i0
prescription. In this future-shell convention the numerator coordinates map as

```text
g_E -> -g_M,  u_E -> +i*u_M,  D_E -> -D_M.
```

If I_M denotes the native integral with this transformed numerator, the Euclidean
measure in the oracle gives

```text
I_E = (-1)^(sum a_e) * (4*pi)^(-L*D/2)
      * (exp(gamma)*Lambda_bar²/(4*pi))^(L*eps) * I_M.
```

Dividing by the unexpanded reference normalization
(4*pi)^(-2L)*(Lambda_bar/2)^(2L eps) leaves
(-1)^(sum a_e)*(4*exp(gamma))^(L eps)*I_M. This is derived for the ordinary
vacuum term, including its numerator mapping. It is a testable normalization
identity, not an assumed blanket conversion for occupied cuts. Each occupied
factor must retain the explicit -theta(mu-E)/(2E) measure above and the remaining
virtual-loop measure. Existing final-state phase-space formulas have a different
normalization and a total-channel constraint; they cannot be substituted for an
independent product of occupied measures.

The full massive assembly is required to validate this adapter numerically. The
massless oracle alone cannot detect an omitted vacuum term.

## Current native input path and supported scope

The native Rust API is `finite_density::DensityInput::prepare`, returning
`PreparedDensityInput`. It validates supplied scalar incidence with Linnet,
checks exact routing and conserved nonbranching species cycles, completes affine
inverse-propagator/medium coordinates with Symbolica, and retains independent
line masses for differentiation. `cut_decomposition(budget)` enumerates connected
complements within an explicit subset budget. `preparation_report(budget)` emits
those definitions and states that numerical evaluation is unavailable.

A runnable example is:

```sh
nix develop --command cargo run --locked --release --bin rustflow -- \
  finite-density-prepare examples/finite_density/massive_two_loop_sunset.json
```

The JSON fields are `name`, `loops`, `vertices`, `edges`, `loop_charges`,
`chemical_potentials`, `targets`, `numerator_convention`, `laurent_orders` and
`digits`. Each physical edge supplies directed `vertices`, exact rational
`routing`, `mass_squared` and species `charges`. Targets contain physical-edge
`powers` and a native Symbolica `numerator`. The scalar coordinate names are
`g1_1`, `g1_2`, ... and medium coordinates `u1`, `u2`, ... . Completion slots are
created algebraically; they are not extra graph edges. Numerator convention is
`shifted_euclidean` or `unshifted_euclidean`. Unknown JSON fields are rejected.
Names are user labels and never select a topology-specific evaluator.

The current JSON charge representation uses integer `loop_charges` and unit
signed edge species charges; an edge can carry at most one species. Exact
rational routings are admitted only when compatible with these integer loop
charges and conserved nonbranching cycles. A rational change of loop basis that
requires fractional loop charges is currently unsupported, so this interface
does not yet satisfy unrestricted arbitrary-routing support.

The Python host exposes `prepare_finite_density(input_json, cut_budget=65536)`
with the same preparation report. A successful preparation does not establish
contour admission or evaluation accuracy. The documented limitation is explicit:
there is no complete native finite-density AMF evaluator in this implementation
stage. Weighted derivative closure, integrated recursive soft boundaries,
complete-amplitude continuation and the numerical acceptance gates remain open.
Standalone polynomial compact-shell seeds do not evaluate retained virtual
poles or assemble a multiloop amplitude.

`finite_density::compact::CompactShell` provides separate seed calculations for
polynomial one-shell moments with d^d q/(2 pi)^d measure. `moment` integrates
energy and radial monomials in its stated convergence domain. `raised_moment`
uses physical-mass differentiation of a simple occupied line, including upper
surface distributions. Its real energy insertion is E^r, corresponding to an
original Euclidean numerator (q0/i)^r; an original Euclidean medium insertion
q0^r requires the factor i^r after q0=iE. This seed convention is not a complete
amplitude normalization adapter. Below support it returns zero. At a positive-mass
threshold it returns zero only when the derived threshold exponent proves a
continuous zero limit; singular/noncontinuous cases receive an explicit
exponent-based diagnostic. A massless infrared-divergent seed requires a separate
dimensional continuation and is rejected. These scope limits do not validate
arbitrary singular test amplitudes at threshold. Working-precision stability must
be checked independently of the bounded compact-series stopping criterion.

| Mandatory result | Current scope and missing evidence |
| --- | --- |
| Complete nonfactorized small graph | Input and cut combinatorics available; no assembled 0/1/2-cut numerical amplitude |
| Generic weighted reduction | Guarded sources/replay infrastructure available; derivative closure for the requested vacuum families unresolved |
| Integrated occupied boundaries | Projected factors and matching extensions available; generic retained soft poles and recursive weighted values unresolved |
| Regulated common-contour continuation | Sufficient positive-real-t admission derived for the massive sunset and its independent-mass neighborhood; general admission and pinch/branch resolution unimplemented |
| Full massive normalization | Ordinary vacuum mapping derived; full massive assembly and independent numerical check missing |
| Physical endpoint reconstruction | Existing ordinary rational projection retained; finite-density evaluation-path cancellation/divergence/depth integration tests missing |
| Three four-loop families | Exact distinct graph certificates and target definitions available; predictions and required supplemental raised-line references absent |
| Ten stable Laurent digits | Precision criteria recorded; no finite-density multiloop stability report or numerical oracle comparison |
