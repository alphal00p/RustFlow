# Prescribed physical routes in the boundary cache

This API extends the existing `RustFlow<S>` and mutable `RustFlowCache` to planner-generated threshold detours. It does not introduce another cache or numerical solver. The existing `evaluate_to` API retains its real, regular, straight-path restrictions.

Bind a `PhysicalContinuation` with `flow.with_prescribed_continuation(...)`. It contains exact ordered polynomial prescriptions in the physical variables, the explicit default side for otherwise unprescribed real singularities, and a continuation-domain label. This typed record is separate from the original branch-domain label and from the dense/canonical system owner. Its representation discriminator and exact contents join the cache identity. A flow bound this way requires `evaluate_prescribed_to` and a mandatory `HomotopyAdmission` callback; ordinary `evaluate_to` rejects that routed identity.

The request is a value in the named physical continuation domain. The admission callback must establish that the particular candidate-to-target route belongs to that integral's global branch convention. Root signs do not determine logarithmic monodromy: a loop for `Y'=epsilon/s Y` can change the answer while leaving coordinates and all root signs unchanged. Arbitrary additional windings are outside this entrypoint. Use low-level explicit-contour transport and a distinct domain when those are intended. Geometry and root-sheet checks do not certify a caller's global physics convention.

The API currently accepts regular exact real rational physical endpoints and an affine kinematic chart. The existing `PrescribedContour` planner supplies disjoint triangular detours around its exact singularity polynomials. Original source/reduction conditions and algebraic norm domains remain obstacles even if matrix coefficients cancel. The planner retains its exclusions for multiple prescribed roots, conflicting local signs, unsupported non-real polynomial coefficients and unresolved root geometry. Endpoints on singularities are errors.

The actual rounded contour vertices become exact binary-rational Atoms and remain fixed across independent precision/order runs. Accepted root Taylor charts come from the existing shared ODE engine; rejected trials do not emit metadata. Endpoint roots are checked against the explicitly requested germ. A saved checkpoint is eligible only when both its solution coefficients and discrete root sheet agree with an independently refined run. Each retained point stores the exact kinematic image of its rounded parameter and its own transported germ. Nonreal contour-detour endpoints are omitted from the physical bank.

Cached input uncertainty still propagates through the existing weighted disk/Gronwall estimator along every actual complex segment. Its floor comes from source accuracy evidence, not working precision. Extra numerical agreement cannot manufacture new source accuracy. The final bank update remains transactional; cancellation, a wrong destination sheet, failed homotopy admission or insufficient accuracy leaves the bank unchanged. Accuracy fallback uses the same bounded candidate loop as ordinary transport.

The binary format is schema 5. It stores the optional typed `PrescribedAffineV1` convention beside the existing identity and rejects earlier schemas explicitly. Loading reconstructs and validates the exact declarations before checking the identity digest. No old bank is silently reinterpreted. Source/dependency fingerprints, payload integrity checks and atomic replacement remain enforced.

Focused regressions verify analytic logarithmic and square-root threshold crossings at a 20-digit requested mixed-error target, transported checkpoint sheets, progressive reuse, independent refinement, binary restart, wrong-sheet rejection, homotopy refusal, cancellation, and original-domain obstacles. Existing regular-path and accuracy-fallback tests continue to use the shared implementation. The full thirteen-component five-point cache acceptance is an explicit ignored test; it is not part of the default suite, and its scientific result is recorded separately. No performance claim is made here.

The accuracy convention is `error <= 10^(-digits) * max(1, |coefficient|)`. It is an absolute target for coefficients below unit magnitude and a relative target above it; requested digits are not significant digits for arbitrarily small coefficients.

Root geometry reuses Symbolica isolation certificates. When a broad certified-real disk from one polynomial obscures a separate crossing, the planner refines that obstructing certificate a bounded number of times, then repeats the same exact clearance and rounded-triangle checks. Complex obstacles remain conservative. Unresolved separation returns an accuracy error. Native isolation/refinement calls are synchronous; cancellation is checked around them.

For a regular endpoint `e`, endpoint classification uses an exact exclusion disk.
Symbolica shifts the polynomial to `P(e+z) = c0 + sum(c_k z^k)`; the radius
`delta = min(1, |c0| / (2 sum_{k>0} |c_k|))` guarantees
`|P(e+z)| >= |c0|/2` throughout the disk; a nonzero constant uses radius one.
A real root enclosure that overlaps an
endpoint is refined to one quarter of this bound before repeating the strict
inside/outside check. Thus a nearby regular endpoint does not fail merely because
a fixed number of refinements was exhausted. Exact endpoint zeros still fail.
The [endpoint regression report](../reports/validation/2026-10-05-contour-endpoints.json)
records the original linear-threshold failure and its focused checks.

A separate [Symbolica-only root-isolation reproducer](../repros/symbolica-cross-factor-isolation/README.md)
preserves a bounded slowdown on a degree-eight polynomial from the nonplanar
Higgs+jet path. The planner now asks Symbolica to factor each exact polynomial
before isolating its roots. Every factor's real and complex root enclosure stays
an obstacle; the original full polynomial still determines its prescribed side.
This avoids unnecessary native separation between unrelated factors without
weakening the subsequent geometry checks. The new regression certifies all
positive roots and full-polynomial signs of the actual degree-eight example.
The side certificate evaluates the full derivative modulo the root's native
defining polynomial, shifted exactly to its enclosure center. Directed interval
rounding then certifies its sign, including when nearby complex roots make the
unshifted polynomial badly conditioned. Native cached roots can have a different
variable label; rebinding that label preserves their exact coefficients.
The [factorization report](../reports/validation/2026-10-05-contour-factorization.json)
records the successful complete NP contour geometry. Numerical physical
continuation remains a separate validation step.

For entire differential systems with no pole list, the shared solver initially proposes the full remaining segment, then uses up to twice the last accepted step length. This avoids repeatedly rejecting an oversized proposal at each expansion point. Taylor-tail, midpoint, endpoint and branch checks remain unchanged. Regressions cover the zero connection, exponential and polynomial systems along changing complex directions, and a sparse Taylor series with a misleading zero tail. Step-budget exhaustion remains a limit error, so cache accuracy fallback does not retry it.

For a system whose exact epsilon-zero coefficient matrix vanishes, the flattened
hierarchy is strictly triangular. With N retained epsilon coefficients, every
product of N connection matrices vanishes, including when physical blocks do
not commute. The inherited-error estimator reuses the same disk norm bounds,
accumulates their integral L over the entire path with fixed source weights,
and bounds growth by the finite sum `sum(L^j/j!, j=0..N-1)`. It does not multiply
a truncated scalar polynomial separately for every numerical step. Any exact
nonzero epsilon-zero term retains `exp(L)`. This sharper estimate does not
change supplied boundary evidence, numerical comparisons or accuracy caps.

The full thirteen-component planar one-loop PH1-to-PH6 cache acceptance passed
for all 65 coefficients through epsilon^4 at a 20-digit absolute reference target.
The same supplied 128-digit source evidence was transported at 50 then 70 working
digits with orders 64 and 96. The largest discrepancy from the original DiffExp
output is 6.17e-23; the discrepancy from the ancillary reference is 4.26e-61.
The scientific claim remains 20 digits. The bank holds 427 physical entries,
including the source and accepted checkpoints/destination; different provenance
keys need not mean distinct coordinates.

The report retains the earlier conservative geometry and exponential-error
rejections as well as the final pass, all 65 values and propagated errors, exact
contour vertices, source/binary/log hashes, and precision metadata. See the
[validation report](../reports/validation/2026-10-04-prescribed-cache.json).
This is supplied-boundary transport, not automatic boundary construction or an
amplitude calculation. The scientific test remains explicitly ignored by the
default suite and is run separately; central release verification is distinct
from the recorded standalone scientific trial.

The [JSON CLI](cli.md#prescribed-threshold-continuation) exposes the same typed
continuation. Its required admission declaration asserts that all planner routes
belong to the stated domain; applications needing route-specific decisions use
the native callback. The [integrated validation report](../reports/validation/2026-10-04-prescribed-integration.json) records release checks and the separately executed scientific tests.
