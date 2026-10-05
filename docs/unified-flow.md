# One numerical engine for auxiliary and physical flow

RustFlow aims to evaluate a reduced integral family independently of its loop
count, external multiplicity and mass pattern. This is an architectural target,
not a claim of current upstream parity. Automatic boundary construction and
singular algebraic endpoints still restrict the supported problem classes.
An exactly singular requested point may require an asymptotic expansion instead
of a finite numerical value. Reducer capabilities and resource budgets remain
separate from numerical-solver capabilities.

## Reduction to numerical evaluation

The native HEPKit graph/model layer owns graph routing, momentum bases and tensor
operations. RustRed or another `ReductionBackend` supplies exact integral
relations, ordered masters and nonzero conditions. RustFlow differentiates the
masters, requests the additional reductions needed for closure, constructs the
connection, obtains boundaries, and evaluates the requested target combinations.

A reduction of the original targets alone does not necessarily include those
derivative integrals. Auxiliary-mass deformation also changes the family, and
boundary recursion can introduce further families. A supplied-table collection
must identify each exact family, denominator order, parameter namespace and
dimension. Missing coverage must produce an incomplete-reduction error. It must
never reuse a table just because another family has the same number of slots.
The backend callback remains available to extend reductions on demand.

`ScopedTableBackend::from_table` binds one supplied reduction to a family;
`insert` adds further family entries. Admission retains uncanceled denominator
conditions and validates the table through the existing native graph planner.
The legacy `TableBackend` leaves family dispatch to its caller. The scoped
adapter adds no RustRed arity restriction to externally supplied tables.

Physical-family preparation also accepts `skip_reduction` and `refine_basis`.
The former retains requested integrals as differential coordinates and reduces
their derivatives; dependent coordinates remain valid when the supplied
boundary satisfies their relations. Independent partial derivatives enter
closure separately, even when their sum cancels. Both closure modes retain the
nonzero conditions from every reduction round.

Basis refinement uses the same candidate row swaps as AMF. Each exact basis
change is applied to every coordinate's connection, including its own derivative
of the transformation. `basis_transformations()` exposes the ordered changes;
`basis_refinement()` reports whether the candidate reduction coefficients have
factorized dimension denominators. This does not assert an epsilon-canonical
form for every connection matrix. Symbolica factors the exact denominators over
its native Gaussian field, so the imaginary unit is a coefficient, not an extra
physical variable. Master order follows the resulting basis and need not be
lexicographic; target projection uses that actual order.
Automatic seeding and projection require epsilon-independent family kinematics
and propagators for their conservative pole bound. Exact differential preparation
can accept epsilon-dependent input, but its pole range needs separate analysis;
the seed API rejects that case before reducing or mutating the cache.

An optional [diagonal epsilon rescaling](epsilon-shearing.md) now regularizes
rational and registered-square-root physical matrix poles when one integer-power
transformation works for all coordinates. Exact root-quotient normalization
precedes epsilon valuation and expansion, preserving original norm domains and
root germs. The supplied-system and prepared-family interfaces retain explicit
scaled basis identities, original domain conditions and coefficient errors.
They request extra source orders before restoring targets; this is separate
from dimension-factorizing refinement. A genuine equal-mass bubble basis with
a `1/epsilon` coupling has independently certified native IBPs and passes the
AMF-seed, physical-cache and target-projection route against a fresh direct AMF
evaluation. General nondiagonal epsilon regularization remains incomplete.

`PreparedPhysicalFamily` already constructs all requested invariant and mass
derivatives in a common closed basis. Its `seed_cache` method calls the ordinary
`solve_integrals` workflow on that basis, retaining the resulting coefficient
uncertainties in `RustFlowCache`. Thus AMF can already provide a direct result or
a boundary for subsequent physical transport within the supported classes.

The [unequal-mass sunset regression](../tests/multimass_physical.rs) exercises
this entire public route for two connected loops. It varies two mass squares
from `(1,2)` to `(3/2,5/2)`, checks scalar and raised-power targets through
epsilon zero against a fresh direct AMF solve, and confirms an exact cache hit.
The target is 20 digits, with absolute scaling near zero. This validates one
nontrivial family, not unrestricted automatic boundary construction.

## Shared solver, different preparation

For coordinates including auxiliary masses and physical parameters, a path
reduces the connection to

```text
dI/dt = [ (dη/dt) Aη + Σa (dxa/dt) Aa ] I.
```

The desired common interface describes this pulled-back system, an ordered
basis, normalization, boundary data, branch prescriptions and accuracy target.
AMF supplies an infinity boundary and physical endpoint projection; physical
transport supplies a kinematic path and usually a cached boundary; FT supplies
a propagator-combination parameter and a final integration functional. These
preparation and endpoint operations need distinct policies, while continuation,
series arithmetic, error propagation and reusable results share numerical owners.
General simultaneous auxiliary/physical paths are a design direction, not a
completed automatic API.

The existing `ode::SeriesSystem` and `transport_series_observed` already serve
ordinary AMF/FT systems, direct epsilon hierarchies and registered algebraic-root
charts. Algebraic roots carry accepted branch state and use Symbolica series
operations. Ordinary and epsilon-channel finite-polynomial Taylor recurrences
now use one implementation without a dense augmented matrix. An independent
dense rational-series reference checks the shared recurrence with noncommuting
epsilon-zero, epsilon-one and epsilon-two matrices at complex expansion centers.

AMF evaluates exact nonzero rational epsilon samples and reconstructs Laurent
coefficients after applying target weights. It checks additional samples and
independent precision/order profiles. Direct physical transport can instead
evolve the Laurent coefficients, requiring an epsilon-regular matrix either
directly or after the supported common diagonal rescaling.
These are coefficient representations of the same differential-equation method.
A general fallback for epsilon-singular physical bases remains outstanding.

## Single-mass vacuum recursion

The recursive boundary provider can remove one loop from a single-mass vacuum
integral by homogeneity and an analytic radial integration. After an exact
routing makes the massive momentum the first loop momentum, the remaining
massless problem has one fewer loop and an external momentum with square −1.
The same boundary provider evaluates that child and memoizes it, including reuse
between different powers of the massive denominator. The radial gamma factors
retain the routing determinant and denominator normalization.

This route follows AMFlow's `SingleMass` ending scheme and has no topology or
loop-count list. Symbolica owns exact matrix and polynomial operations; native
HEPKit kinematics and family completion own the child scalar products and any
missing zero-power numerator slots. Its current admission requires nonnegative
powers, one active exact rational-complex mass square off the nonpositive real
axis, and real rational rank-one quadratic denominators with positive
normalization. The principal mass power analytically continues the radial
formula from positive mass square; the lower-loop massless child is unchanged
and reused across different masses and massive powers. Values directly on the
mass cut need a separate limiting prescription and retain the existing fallback.
Other unsupported cases also retain the existing
fallback. Custom terminal providers and certified scaleless sectors take
precedence.

The radial boundary formula and the automatic AMF contour have separate branch
contracts. For recognized single-mass vacuum targets, including one-line
tadpoles with positive quadratic normalization, the default imaginary-axis
contour is rejected when `Re(M²) <= 0` and the nonzero imaginary mass has the
incompatible sign: positive for `PlusI0`, negative for `MinusI0`. This check runs
before reduction and symbolic-cache lookup, and is retained and rechecked if
evaluation options change. Compatible conjugate contours pass principal-power
and precision-refinement checks. Native certified-zero sectors bypass this
irrelevant sheet constraint. This is a deliberately scoped guard, not a
general complex-mass homotopy planner. A pole-free contour alone need not end on
the requested sheet; direct values on the negative real mass cut also remain a
separate limiting-prescription problem.

The [finite-sample comparison](../reports/validation/2026-10-05-single-mass-vacuum-finite-samples.json)
checks a connected three-loop six-line vacuum and its raised massive power
against the separate FT route at two epsilon values and two precision/order
settings. Those checks validate the new recursive formula; a full Laurent fit
through the public route now also [passes through epsilon zero](../reports/validation/2026-10-05-single-mass-vacuum-laurent.json)
for both targets at 20 requested digits. It uses 33 samples checked against a
disjoint 29-sample grid, with higher working precision and series order. A
separate public FT solve agrees on all 14 coefficients; the largest difference
is 7.20e−54. This is refinement and cross-method evidence, not an interval proof.
The full Laurent comparison is between the two native routes. The recorded AMF and FT times are
120.115 s and 52.336 s with different precision/order profiles.

A separate [live original AMFlow check](../reports/validation/2026-10-05-single-mass-vacuum-upstream.json)
evaluates both targets automatically at exact epsilon `1/10`. Its 40-digit outputs
pass a 35-digit comparison with the native finite sample, including the raised
power checked through exact homogeneity. The original run took 175.447 s with
Kira and Mathematica; that finite-sample time is not directly comparable to the
native full Laurent timings above. The optional `single_mass` Laurent regression
runs without Mathematica and preserves the independent native fit comparison.

## Current generality gaps

- Region enumeration, normalization and recursion have explicit search budgets.
  Boundary routing currently expects rank-one quadratic momentum propagators;
  linear and cut automatic boundary support covers only part of upstream's scope.
- Rational Frobenius matching supports logarithms and resonances, but requires
  supported indicial roots. General algebraic Puiseux/log endpoint matching is
  incomplete; the new algebraic origin initializer covers compatible analytic
  sectors.
- Standalone continuation accepts complex paths. The reusable algebraic physical
  cache accepts exact rational-complex coordinates on regular affine segments,
  certifying root-sheet transitions with native exact polynomial arithmetic.
  General algebraic coordinate constants, automatic causal-path inference,
  degenerate Gram charts and algebraic singular endpoints require further work.
- General Möbius/Padé transport controls remain unfinished.

Benchmarks test these interfaces; family-specific fixtures do not define the
solver's mathematical scope. New coverage must exercise the same public route
from reductions to boundaries, transport, target projection and uncertainty.

## Relevant recent work

[AMFlow 2.0, arXiv:2607.08477](https://arxiv.org/html/2607.08477v1)
motivates the existing FT recursion, basis refinement, symbolic restart cache,
skipped initial reduction and compiled solver work. Its FT mode only guarantees
Euclidean evaluations because a general physical contour prescription is lacking.
That limitation must not be concealed by a shared API. Broad validation of all
these features remains necessary.

[Analytic regression, arXiv:2507.17815](https://arxiv.org/html/2507.17815v1)
uses high-precision integral samples to reconstruct rational combinations of a
chosen function basis. Multipoint lattice reduction trades sampling count against
precision, with conditioning limiting the gain. A progressively populated cache
could support that downstream use, with fresh points reserved for validation.
It does not supply a general replacement for AMF boundary recursion.

[SOFIA](https://github.com/StrangeQuark007/SOFIA) is a singularity-analysis
reference. Analytic reconstruction and SOFIA integration are not active parity
gates or new runtime dependencies. Symbolica already provides integer-relation
machinery; any future extension must audit existing implementations before adding
new algorithms. Generic numerical evaluation remains the priority.
