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
cancellations. The finite-density integration test now prepares the actual
massive one-cut sunset from generated guarded sources, changes two components
exactly to `y_j=I_k`, `y_k=I_k-eta*I_j`, and reconstructs the original targets
with `(y_j-y_k)/eta`. Its infinity boundary is the original integrated occupied
region data, matched after the exact polynomial `z=1/eta` transformation of
every known power/log coefficient. Shared transport reproduces the original
targets at start scales 8 and 12. The same source-derived solution rejects an
uncancelled target pole and a target whose pole exceeds the known endpoint
series depth. These integration checks pass in
`src/finite_density/flow_endpoint_tests.rs`; the historical run is recorded in
`reports/validation/2026-10-09-finite-density-native-assembly/positive-energy-unit-resources.json`.

A valid calculation must establish a common convergence domain, continue the
complete regulated amplitude, remove the auxiliary parameter, take the specified
physical-mass limits with any necessary regions/distributions, and only then
extract the requested Laurent coefficients. Exchanging these operations requires
proof for the admitted family. The ordinary dimensional endpoint rule does not
by itself provide such a proof for occupation-weighted integrals. No general
four-loop admission theorem or successful evaluation is claimed here.

For the implemented **strict positive-mass** sunset and heavy-edge domains,
the endpoint admission is stronger than a dimensional-limit convention. The
strict bounds in the contour certificates persist in a complex neighborhood
of eta=0. Compact shells and their Fermi surfaces stay fixed; their masses are
positive. After introducing independent analytic indices to control virtual
ultraviolet integrals, the uniformly nonvanishing virtual denominators permit
Taylor expansion in eta under each compact/distribution pairing and virtual
integral. Meromorphic continuation in those indices and D preserves this
local analytic dependence. Finite derivatives of the physical line masses
and polynomial medium numerators preserve it as well. Thus the physical
solution in these domains has an ordinary Taylor endpoint. Noninteger
epsilon-dependent homogeneous solutions of its differential equation have
zero physical coefficients; they do not supply missing physical regions.
The shared projector selects the constant only **after** multiplying and
summing the exact target weights, including cancellations of apparent poles
introduced by reduction or a change of basis. This argument uses the strict
mass/contour margins and does not extend the admission to massless shells,
threshold pinches, or unproved changes in the order of limits.

## Whole-amplitude normalization check

For a contribution with physical powers a_e, define a native virtual loop with
d^D k/(i pi^(D/2)) measure. Begin with the
standard Wick map k0=i P0, then reverse all native Minkowski momenta. We choose
q_M=(-i*P_E0,-P_Espatial), so the Euclidean crossed pole P_E0=+iE maps to future
q_M0=E. The global Minkowski loop reversal preserves its measure and +i0
prescription. In this future-shell convention the numerator coordinates map as

```text
g_E -> -g_M,  u_E -> +i*u_M,  D_E -> -D_M.
```

Define each native occupied loop with the independent measure

```text
d^D q/pi^(D/2) * theta(q0) * theta(mu-q0) * C_n(q²-m²),
C_n(x) = (-1)^(n-1)/(n-1)! * delta^(n-1)(x).
```

This measure has no total timelike channel or final-state momentum-conservation
delta. A simple shell integrates to theta(mu-E)/(2E). The original Euclidean
occupied sign is supplied by (-1)^n on that physical index; higher indices follow
by independent physical-mass differentiation. Applying an additional per-cut
minus sign would double-count it. The same quadratic phase applies to virtual
denominators and negative scalar numerator indices. Medium slots instead retain
their exact powers of i from the stated Wick map.

For L original loops and k occupied loops, direct comparison of d-dimensional
spatial measures gives

```text
native_measure_to_euclidean = (2*pi)^k * (4*pi)^(-L*D/2).
```

The extra 2*pi per occupied loop comes from replacing an energy contour measure
by the independent shell distribution. If I_M includes the transformed numerator
and native C_n distributions, the Euclidean measure in the oracle gives

```text
I_E = (-1)^(sum a_e) * (2*pi)^k * (4*pi)^(-L*D/2)
      * (exp(gamma)*Lambda_bar²/(4*pi))^(L*eps) * I_M.
```

Dividing by the unexpanded reference normalization
(4*pi)^(-2L)*(Lambda_bar/2)^(2L eps) leaves
(-1)^(sum a_e)*(2*pi)^k*(4*exp(gamma))^(L eps)*I_M. Exact routing Jacobians
are separate from these measure and index factors. Existing final-state
phase-space formulas have a different normalization and a total-channel
constraint; they cannot be substituted for independent occupied measures.

`finite_density::normalization` implements the measure conversion, quadratic
index phase, polynomial Wick map and unexpanded MSbar scale conversion as
separate exact Symbolica functions. `OccupiedCutFamily` already includes target
Wick/index phases in its coefficients; callers must apply only the measure and
routing factors to those targets. `IntegratedOccupiedBoundary` returns the native
mixed measure above. `CompactShell::raised_moment` instead returns an original
Euclidean occupied seed, including its physical-line phase and d^d q/(2*pi)^d
measure. It must not receive the native conversion again.

The normalization regressions exercise massive one-loop vacuum plus occupied
assembly, raised powers, below/at/above thresholds, and actual mixed and fully
occupied sunset large-mass coefficients. At D=3, m1=m2=1/2 and mu=1, the leading
Euclidean mixed and fully occupied sunset coefficients are respectively
-1/(64*pi²) and +1/(64*pi²), obtained independently from occupied and uniform
Gaussian seeds. These are boundary coefficients, not the sunset at its physical
endpoint. Complete massive assembly is a separate validation gate; the selected
sunset now passes it at D=12/5, including its nonzero vacuum term. A massless
oracle cannot detect an omitted massive vacuum term. Current command results are
recorded in the separate
[`native-assembly` validation report](../reports/validation/2026-10-09-finite-density-native-assembly/README.md).

## Current native input path and supported scope

The native Rust API is `finite_density::DensityInput::prepare`, returning
`PreparedDensityInput`. It validates supplied scalar incidence with Linnet,
checks exact routing and conserved nonbranching species cycles, completes affine
inverse-propagator/medium coordinates with Symbolica, and retains independent
line masses for differentiation. `cut_decomposition(budget)` enumerates connected
complements within an explicit subset budget. `preparation_report(budget)` emits
schema-version-2 definitions, states that numerical evaluation was not performed,
and identifies the separate numerical interfaces and their documented scope.

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
contour admission or evaluation accuracy. The native
`finite_density::assembly::PreparedDensityFlow` owner attempts the complete cut
sum at a common regulator sample and reconstructs Laurent coefficients only
after summation. Shared evaluation entry points are `rustflow finite-density`,
`rustflow finite-density-sample INPUT.json EPSILON` and Python
`evaluate_finite_density(input_json, epsilon=None)`, with detailed controls in
the [CLI and Python protocol](cli.md#native-finite-density-input-and-evaluation).
They retain raw Euclidean normalization and propagate any failed required sector.
The runtime weighted-AMF dispatch preserves exact storage arities 7 and 9,
then uses capacities 12, 16, 20, 24 or 32 for other physical arities that fit.
Physical factors are never padded: only native index arrays acquire verified
fixed-zero tails, which are removed before boundary and transport interfaces.
The const-generic Rust API can compile larger capacities, subject to physical
admission and closure. See the [capacity invariants and validation](finite-density-native-capacity.md).
The default closure requests provisional frontier sectors, allows 12 rounds,
and uses native search depth 3 with 8192 domains. Guard refinement is bounded
to 3 additional passes, 256 added domains and interval width 2. Before that dispatch, a strict disjoint-support proof
can return zero, and `PreparedOccupiedTerminal` handles surviving polynomial
noncut factors using native compact moments without an auxiliary parameter or a
compiled guarded index arity. It retains unrestricted polynomial virtual factors
only when their dimensional integral is scaleless. It does not replace retained
virtual poles by moments.

For assigned strictly positive rational shell masses, weighted AMF uses the
open massive-sunset or routing-dependent heavy-edge certificate, with the common
+i0 prescription and all uncut physical quadratic factors shifted. A separate
sealed massless permit covers zero physical masses, positive occupied chemical
magnitudes and polynomial completions. Singleton cuts admit a uniformly
UV-continued external-invariant germ at arbitrary admitted virtual rank, with
every uncut momentum depending on virtual loops. Two-cut sectors admit pure
spacelike transfer factors and independent rank-one virtual quadratic blocks
whose external direction is q1-q2. Exact routing, all-uncut deformation, origin
convention and every generated finite basis label are bound to the
[endpoint proof](finite-density-massless-endpoint.md) and
[uniform germ/block derivation](finite-density-massless-germ.md). These classes
are selected by physical structure, not graph names or targets. All E7 cut
sectors pass this structural admission, but their bounded native weighted
closure remains unresolved and no E7 amplitude has been evaluated. General
coupled virtual blocks, E8 pinched continuations and partial deformation
placements remain unsupported. Compact terminals retain their
separate distributional admission. The algebraic parser can still represent
families outside numerical admission.

The optional Rust `PreparedOccupiedFlow::prepare_with_source_options` path can
admit integer inverse powers of completion factors proved exactly equal to
`c*E_i`, with nonzero rational c and an occupied shell of assigned mass m_i>0.
Its source and boundary domains retain the certificate `E_i>=m_i>0`; virtual
energies and sums of energies receive no such admission. Signed energy moments
use the original numerator before independent mass differentiation. This
option also enables explicit shell-normal source presentations derived from
the smooth vector `u/(2E_i)`. The default remains the validated Lorentz
presentation with polynomial completions. The option and exact certificates
participate in program identity and evaluation provenance. It is an extension
of the intermediate integral domain, not permission to replace a virtual pole
by a compact moment. The successful selected two-cut system below uses the
default polynomial completion domain.

The complete massless two-loop sunset now passes all four native profiles at
D=13/2 and D=15/4, and five Laurent profiles through order zero. All occupied
sectors use native weighted AMF, with closed bases 2, 2 and 3. The 110 independent
reference comparisons and 84 refinements retain the raised original numerator
and its nonzero upper surface; the vacuum and both singleton endpoints are
checked as zeros. Native D=7 retains an explicit additional indicial resonance
before producing a prediction. In the current implementation, native vacuum
zero certificates skip vacuum flow and both singleton cuts succeed; the
remaining resonance occurs in the double-cut infinity recurrence. Independent
convergent D=7 quadrature validates the reference and support derivative. See the
[massless numerical report](../reports/validation/2026-10-09-finite-density-native-assembly/massless-native-validation.json)
for counts and the [provenance record](../reports/validation/2026-10-09-finite-density-native-assembly/massless-native-provenance.json)
for the original-source executable-digest limitation and final-order check.
The current structural extension was revalidated at D=13/2 and through Laurent
order zero: 70 independent references and 54 refinements pass. A separate
massive D=12/5 profile passes ten independent references and reproduces all ten
historical component/total decimal strings exactly. Its current
[numerical provenance](../reports/validation/2026-10-09-finite-density-native-assembly/massless-blocks-numerical-provenance.json)
captures each executable digest before launch; the
[validation aggregate](../reports/validation/2026-10-09-finite-density-native-assembly/massless-blocks-numerical-validation.json)
keeps these regression counts separate from four-loop acceptance.

The complete generic numerical evaluator has not passed its required acceptance
gates. The selected massive sunset now passes complete numerical assembly at
epsilon=4/5 (D=12/5), including its vacuum, both single cuts and the two-cut term
for the scalar and raised original-numerator targets. All 40 independent
reference comparisons and 30 precision/order/start-scale refinement comparisons
pass, with 40 guard digits held fixed. Native frontier-sector discovery closes
bases of sizes 6, 6 and 11. These are closed generating sets, with no minimality
or full index-domain coverage claim.

The earlier 64-element N9 basis also evaluates successfully when its guard
precision is increased from 20 to 40 at otherwise unchanged numerical settings.
Its saved predictions agree with the reduced bases at the matching profile;
the earlier endpoint-divergence diagnostic was not evidence that the original
massive integral diverges. Historical failures and their exact resource records
remain in the validation report. Generic four-loop closure, integrated recursive
soft boundaries and the remaining numerical acceptance gates stay open.
Independent references now also cover both complete massive targets through
Laurent order zero, with separate vacuum, cut and Fermi-surface coefficients;
the [subtracted reference report](../reports/validation/2026-10-09-finite-density-native-assembly/independent-laurent-reference/README.md)
passes quadrature, epsilon-grid and precision refinements. Complete native AMF
Laurent coefficients through order zero now pass all 30 independent-reference
comparisons and 24 precision, order, start-scale and epsilon-grid refinement
comparisons across five profiles. The largest relative reference discrepancy is
1.430e-18 and the largest refinement difference is 9.250e-45; these are empirical
agreements, subject to the independent reference accuracy estimates.
The public sample and Laurent commands have also evaluated the complete massive
one-loop tadpole above and below support. Its scalar and raised sample values
and four Laurent coefficients pass independent analytic comparisons; see the
[interface validation report](../reports/validation/2026-10-09-finite-density-native-assembly/native-interfaces/README.md).
This one-loop success does not complete the required nonfactorized multiloop case.
Public assembly also passes an actual two-species calculation with independent
chemical potentials: both singleton cuts use native source-closed weighted AMF,
and the double cut uses the compact polynomial terminal. Scalar and original
Euclidean `u1*u2` contributions are checked separately in all four sectors and
in their sums. The graph factorizes into two one-loop cycles, so this validates
species-dependent support and Wick normalization without completing the required
nonfactorized sunset. The same release batch passes multihard occupied-boundary
and three-compact-loop coefficient tests, plus public thermal-threshold assembly;
see the [boundary and assembly gates](../reports/validation/2026-10-09-finite-density-native-assembly/new-boundary-gates.json).
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
threshold it returns zero when the derived threshold exponent proves a continuous
zero limit. The finite-jump case uses the regulated thermal limit derived below;
more singular cases receive an explicit exponent-based diagnostic. A massless infrared-divergent seed requires a separate
dimensional continuation and is rejected. These scope limits do not validate
arbitrary singular test amplitudes at threshold. Working-precision stability must
be checked independently of the bounded compact-series stopping criterion.

For a seed with energy power `r`, spatial radial power `j`, and raised line power
`k+1`, write `a=m²>0`, `gap=mu²-a`, `alpha=d/2+j`, and
`A_d=2*pi^(d/2)/(Gamma(d/2)*(2*pi)^d)`. A simple occupied seed near threshold is
`-A_d/4 * a^((r-1)/2) * gap^alpha/alpha` on the occupied side. When
`alpha=k>=1`, regulate the occupation before taking the `k` independent mass
derivatives: with `tau=2*m*T` its leading radial integral is

```text
-A_d/4 * a^((r-1)/2) * tau^k
    * integral_0^infinity x^(k-1) n_F(x-gap/tau) dx,
n_F(x)=1/(exp(x)+1).
```

The `k` gap derivatives give `Gamma(k)*n_F(-gap/tau)`. Dividing by `k!` and
taking `T->0+` at fixed `mu=m` therefore gives
`-A_d/(8*k) * a^((r-1)/2)`. Derivatives of the smooth prefactors vanish with
positive powers of `T`. This specifies a thermal finite jump, without assigning
an arbitrary value to `theta(0)` or claiming continuity.

Independent finite-temperature radial identities test this case. At `d=2`, the
scalar raised seed is `-n_F((m-mu)/T)/(8*pi*m)`. With the original spatial `q²`
numerator, the power-three seed is `-n_F((m-mu)/T)/(16*pi*m)`. These follow by
differentiating the fixed-spatial-momentum thermal kernels before integration;
the unit test numerically integrates their smooth derivative kernels at three
temperatures and three chemical offsets. At `mu=m=1/2`, their limits are
`-1/(8*pi)` and `-1/(16*pi)`. A power-three scalar seed has gap exponent `-1`
at `d=2` and retains an explicit singular-threshold diagnostic.

| Mandatory result | Current scope and missing evidence |
| --- | --- |
| Complete nonfactorized small graph | Full D=12/5 vacuum-plus-all-cuts assembly passes for both targets: 40 independent sample reference checks and 30 refinements; Laurent [-2,0] passes 30 reference checks and 24 refinements |
| Generic weighted reduction | Exact requested-target/derivative closure and source replay pass for all massive-sunset cut sectors with bounded guard refinement; generic four-loop closure remains unresolved |
| Integrated occupied boundaries | Polynomial compact and recursively integrated hard coefficients verified; generic retained soft poles and recursive weighted values unresolved |
| Regulated common-contour continuation | Massive sunset and routing-dependent heavy-edge sufficient domains derived; see [contour certificates](finite-density-contours.md); thresholds, massless pinches and general admission remain separate |
| Full massive normalization | Complete massive-sunset scalar and raised-numerator values pass independent raw-Euclidean comparisons, including the nonzero vacuum and raised surface contributions |
| Physical endpoint reconstruction | Complete massive and massless two-loop values/Laurent coefficients pass; singleton germs and independent rank-one blocks structurally admit all E7 cuts, with native closure/numerics unresolved; D7 double-cut infinity resonance remains explicit |
| Three four-loop families | Exact distinct graph certificates and E7 endpoint admission available; independent E7 raised-line reference generated; all native predictions and the other two complete supplemental references absent |
| Ten stable Laurent digits | Complete massive-sunset Laurent comparison passes five profiles, independently changing digits, order, start scale and epsilon grid; no four-loop AMF-to-oracle comparison |

An independent validation-only reference for both complete massive sunset targets
has been generated at D=12/5 using direct Schwinger/Feynman parameters and compact
quadrature. All 12 components/totals passed node and working-precision refinement;
the generator reads no AMF values or oracle answers. A separate validation-only
script subsequently compared the complete saved native predictions for both
targets across four configurations. All 40 vacuum/cut/total comparisons and
30 independent refinement checks pass the 1e-12 relative criterion; the largest
relative reference discrepancy is 9.449e-36. Agreement is empirical, subject to
the independent reference's refinement estimates rather than a rigorous interval
error bound. The separate complete Laurent comparison also passes all 30
reference checks and 24 refinements; no four-loop numerical gate is completed
by this small-graph result. See the
[reference derivation](finite-density-reference.md) and the
[native-assembly report](../reports/validation/2026-10-09-finite-density-native-assembly/README.md)
for exact measurements, source provenance and failure records.

The independent four-loop E7 reference was subsequently compared with the
supplied I37 reference: all three Laurent coefficient expressions agree exactly
after resolving the supplied constant. This is a
[reference-versus-reference check](../reports/validation/2026-10-09-finite-density-native-assembly/independent-e7-reference/README.md),
with zero AMF predictions read and no native four-loop acceptance claim.
