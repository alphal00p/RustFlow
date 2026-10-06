# RustFlow: AMFlow and DiffExp parity goal

This is the active development goal, not a claim of completed parity. The
original AMFlow paper's four-target 20-digit acceptance passed; see
[coverage](coverage.md). The expanded goal also includes linear propagators,
cut phase-space integrals, all upstream benchmark families, and measured
performance at matched achieved accuracy. RustRed is checked for upstream updates
at validation milestones and newer versions are integrated between builds while
preserving local work. The upstream arity fix is consumed through RustRed’s shared dispatcher; no
independent arity registry is maintained here. Broad coverage, including
the complete Higgs+jet application, takes priority over general performance
tuning; performance blockers may be addressed to make a coverage test practical.

## Current implementation snapshot

As of 2026-10-06, the descriptive Python API and complete native gg→hg notebook
acceptance have passed on the publication wheel. Fresh numerical boundaries,
4,360 transported coefficients, eight W/Z form factors, three coherent
observables, independent 40-digit regeneration and interruption/restart/cache
reuse are covered. The 48/61 parent equations and basis maps are certified
supplied symbolic inputs; independent full parent reduction with RustRed is not
demonstrated. The external EW-square reference supports 19 relative digits.
See [the notebook report](python-notebook-status.md) for exact runtime provenance,
publication status and accuracy limits.

The shared Rust engine also supports recursive Euclidean FT for tested mixed
linear families and local registered-root Frobenius endpoints; their combined
gate passes 74 targeted release tests on top of the preceding 541-test baseline.
Checked rational/epsilon Padé transport and a mixed-cut adapter have since
passed a [combined full release run](../reports/validation/2026-10-06-pade-cut-integration/summary.json):
574 unit/integration tests and one doctest, formatting and strict all-target
Python/stub-generation Clippy, with 13 scientific tests explicitly opt-in.
Mixed-cut boundaries cover complete positive-energy
final states with uniform causal virtual directions, massive two-body or
massless N-body leaves and certified nonzero real-only denominators. Broader
cut recursion, remaining DiffExp controls, constrained endpoint sectors and
broad performance parity remain unfinished. Supplied finite singular endpoints
now integrate with regular transport and a separate terminal cache. Native cut
diagram projection is exposed in Rust, Python and the CLI; bounded adaptive
integer residual arithmetic is opt-in. These three additions pass a
[combined targeted gate](../reports/validation/2026-10-06-endpoint-cut-arithmetic-integration/summary.json)
of 328 top-level release tests across 29 targets, six child-process checks,
formatting and strict Python/stub-generation Clippy. They have not been rebuilt
into the frozen notebook publication wheel. The
[comparison table](benchmark-comparison.md) records measured successes and
remaining coverage, without treating successful benchmark families as proof of
arbitrary kinematic, mass or loop-count support.

The [Python prescribed-contour interface](python-prescribed-transport.md) now
binds native exact polynomial declarations through `ContinuationPrescription`.
It reuses the same continuation, branch and endpoint-cache owners; nonconstant
paths require explicit admission, while exact compatible hits retain source
policy without planning a new route. Its isolated milestone passes 58 release
tests, formatting and strict Python/stub-generation Clippy, including both
logarithm sides, root sheets, cancellation and endpoint restart.
The subsequent main integration passes 62 release tests across nine targets
with the same static checks and unchanged source hashes.

Three further validated milestones extend this coverage. Exact
[correlated rational endpoint constraints](supplied-singular-endpoints.md)
retain an affine space of integration constants and check every reconstructed
matching component. Schema 7 stores their asserted relations and provenance;
constrained registered-root spaces and dimensional-sector selection remain open.
[Partial cut mass placements](cuts.md) reuse the existing region, tensor,
recursive vacuum and phase-space owners, including a nonfactorized three-loop
comparison with original AMFlow. The
[general registered-root Python constructor](python-algebraic-transport.md)
handles dense coupled connections as well as canonical ones. These newer
interfaces are separate from the frozen community notebook runtime.

The [native series-construction comparison](../reports/performance/2026-10-06-native-series-construction/README.md)
uses balanced additions of native Symbolica series. All numerical outputs and
uncertainty data agree exactly with the preceding implementation. Matched
checked first-destination medians improve from 6.160 to 5.464 seconds for 48
masters and from 15.963 to 14.770 seconds for 61 masters. The full 4,360-coefficient
transport acceptance also passes. These are native revision comparisons;
full-amplitude timing against the original plugin and broad performance parity
remain open. Immutable exact connection preparation is the next optimization.

The [unified-flow design and generality audit](unified-flow.md) makes reduced-family
evaluation the central interface. AMF direct evaluation, AMF-generated physical
boundaries and DiffExp-style transport share the numerical engine. Supplied
reductions must cover derivative closure and each deformed/recursive family;
upstream parity is not established by process-specific benchmark success.

## Pinned references

| Reference | Revision | Role |
|---|---|---|
| AMFlow 2.0 | `26005517a288086c4cb4d1b26d829691bc088485` | Automatic auxiliary-mass boundaries, AMF/FT recursion, oracle |
| DiffExp 1.1 | `784c8229bf92369a03f011a48e161522c8c54bbd` | Physical-variable transport, examples, oracle |
| Symbolica community | `c3408e4ba1d3bdd4ea55678fad50e27009be13d4` | Current shared algebra and arbitrary-precision dependency; older benchmark reports retain their measured revisions |

DiffExp is checked out at `/common/dev/diffexp`. Its source is GPL-3.0-or-later;
the Rust algorithms are implemented independently on this project's existing
MIT core, with attribution to the papers and numerical provenance. No DiffExp
implementation is copied or translated into MIT source. The separate
`diffexp_refactored` repository is not the reference selected by the user.

The tool is named RustFlow. Its native reusable boundary bank is
`RustFlowCache`, passed by mutable reference and persisted in binary form.
The native library remains in this repository along with a HEPKit-compatible
DOT-input CLI. Optional PyO3 bindings also live in this crate and register in the
community host under `symbolica.community.hep.integration`, using descriptive
Python class names. See the [Python API and notebook status](python-notebook-status.md).
Prefer existing Symbolica,
Linnet and Spenso/HEPKit capabilities over independent implementations of graph
or numerator-tensor operations. The reference checkout under
`/common/dev/fastsecdec/DO_NOT_PUSH_FOR_REFERENCE_ONLY` is read-only.

The application at
[mg5_higgs_ew_plugin](https://bitbucket.org/aschweitzer/mg5_higgs_ew_plugin/src/master/)
is an additional design and validation reference: repeated mixed QCD–EW
`g g -> h g` evaluations make reuse of intermediate transport points essential.

## Reusable physical boundaries are central

Every accepted local series is useful beyond its final endpoint. Retain its
coordinate chart, traversed interval, branch, epsilon coefficients, and precision
metadata. Physical points sampled along that interval can become starting
boundaries for later evaluations. Complex detour points must remain distinguishable
from physical phase-space points.

A reusable boundary is identified by the complete differential system, ordered
basis, normalization, epsilon range, prescription, and branch domain. Matching
coordinates alone does not establish compatibility. Working precision alone
does not establish accuracy. Unverified intermediate values must not enter a
search requiring verified digits, and cached boundary uncertainty must remain
visible during subsequent transport.

Search should minimize estimated continuation cost among compatible entries.
Scaled kinematic distance is a fallback; pole distances, path feasibility,
required order, and available precision can change which boundary is cheapest.
The design provides a cost policy so applications can include physical cuts and
their own phase-space metrics. A physically disallowed or incompatible path is
excluded, not merely assigned a large distance.

Retained points and series should survive process restart through versioned,
authenticated, atomic persistence. Preserve exact coordinates when supplied;
label coordinates obtained by numerical transport as rounded numerical points.
Do not pass coordinates or values through binary64 during serialization or
search. Cache restarts and repeated MC-like queries need their own numerical
and timing tests, separately from cold boundary construction.

## Implementation gates

1. **Exact physical equations and paths.** Supplied partial-derivative matrices,
   canonical dlog construction, exact chain-rule pullback, and invariant
   derivatives including the external Gram variation. Validate masses versus
   squared masses and preserve nonzero Gram/basis conditions.
2. **Direct epsilon transport.** Integrate the triangular coefficient hierarchy
   using sparse finite-polynomial recurrences; retain common Laurent leading
   powers. Reuse matrix preparation across epsilon orders rather than building
   a dense augmented system. Validate coupled noncommuting matrices.
3. **Reusable transport cache.** Store and search compatible physical boundaries,
   retain series, verify intermediate checkpoints, restart from persistent
   entries, and demonstrate reduced work for repeated destinations.
4. **General boundary matching.** Partial and asymptotic power/log constraints,
   infinity boundaries, compatible boundary-line changes, and underdetermined
   affine solution spaces. Reuse AMF boundaries when no supplied point exists.
5. **Algebraic coefficients and continuation.** Square roots of polynomial
   radicands, per-polynomial iδ prescriptions, Puiseux series, and branch state.
   Include radicand zeros in path planning. Original DiffExp does not support
   higher algebraic roots or elliptic-function prefactors.
6. **Transport controls.** Dynamic/predivision segmentation, convergence-radius
   controls, Möbius coordinates, Padé with explicit fallback/error checks,
   reusable piecewise evaluation, and degenerate paths. Regular rational/epsilon
   [Padé trials](pade-transport.md) now share checked accepted/saved output,
   exact-source defects, fallback and genuinely different refinement orders;
   registered-root candidates remain Taylor.
7. **AMFlow extensions.** Linear propagator regions/terminals and reverse-unitarity
   cut semantics; broaden normalization, boundary recursion and basis refinement
   to the supported upstream surface. Track limitations separately from the
   native backend’s search and representation limits.
8. **Benchmark and performance gate.** Every reference records source, settings,
   input precision, normalization, path, prescription, output precision and
   runtime provenance. Numerical parity precedes timing claims. Compare cold
   preparation, transport, cache restart and repeated-query throughput separately,
   at matched work or matched achieved accuracy.

## DiffExp benchmark ladder

| Family | Dimension / task | Required comparison |
|---|---|---|
| Multiple polylogarithms | Small words and a 20-letter word | Analytic branch values, orders 50/75/100, saved points |
| Equal-mass three-loop banana | 4 masters, 3 coupled | Partial boundary at infinity, -1 then 32, epsilon 0–4 |
| Unequal-mass three-loop banana | 15 masters, 11 coupled | Five variables; independent mass/momentum paths; published 55-digit coefficients |
| Nonplanar massless five-point | 108 masters | `1812.11160v2` ancillary matrices and at least 50-digit references |
| Planar one-mass five-point | 13, 75, 74, 86 masters | 1loop, zmz, mzz, zzz; point-to-point checks including 128-digit zzz profile |
| Higgs + jet application | Mixed two-loop QCD–EW amplitude | Intermediate-point reuse and repeated physical queries |

The bundled banana matrices contain only epsilon orders zero and one; higher
order files are zero placeholders. Preserve the distinction between notebook
settings and paper precision profiles. The planar notebook's displayed URL and
actual download differ in versioning; pin and verify ancillary data before use.
Published timings are historical observations, not present-machine comparisons.

The first live DiffExp smoke oracle uses `G(1,0,1;4)` with an asymptotic zero
boundary and a contour below `t=1`. Orders 50 and 75 agree to 24 digits; the
refined result agrees with the notebook to 30 digits. This is not yet the full
DiffExp benchmark gate or a Rust performance comparison.

The equal-mass banana gate now has a production regression through ε⁴ at both
`-1` and `32`, with analytic partial infinity data and independent precision/order
refinement. See [its provenance and measured scope](banana-equal.md). The growing
cache, native HEPKit DOT input, AMF-generated physical seeds, and two-body cut
terminal also have offline release regressions. These milestones do not complete
the algebraic-kernel, full five-point, full amplitude, or performance gates.

The unequal-mass banana now has a complete 15-master, 75-coefficient production
regression with independent mass/momentum routes and an original DiffExp oracle;
see [its measured scope](banana-unequal.md). Standalone [algebraic transport](algebraic.md)
tracks registered square-root sheets on regular contours, including winding and
rejected-step tests. Root-sum denominators use native formal quotient inversion
and preserve source-domain conditions. [Algebraic cache reuse](algebraic-cache.md)
now supports regular rational-complex affine paths with explicit root germs and
exact principal-cut crossing checks. The
[prescribed-contour interface](prescribed-contours.md) constructs detours from
caller-supplied polynomial prescriptions, including admitted cached threshold
crossings; automatic inference from an integral and general algebraic Frobenius
endpoints remain incomplete.
The [33,000-boundary cache benchmark](cache-selection.md) measures lazy path
validation separately from integral evaluation. Exact target projection now
checks positive epsilon orders and carries cached master uncertainty into the
requested target coefficients.

## Runtime arity and independent limits

Both RustFlow reduction adapters now follow RustRed’s compiled registry. The
factorized path delegates to `rustred::dispatch_arity!`, retaining its existing
native solver and every search option. The sparse path calls the checked runtime
bridge. This consumes RustRed’s upstream `33fd03ec` fix and removes the former
independent 12-slot dispatch mismatch.

The registry defaults to 1–16 slots. `rustred::compiled_runtime_arities()` reports
the actual build. Set `RUSTRED_RUNTIME_ARITIES` before building to select another
finite registry; it is not a runtime option. Cache compatibility includes the
build script and actual compiled registry. Unsupported entries return a typed
capability error. Native sector enumeration is still exponential in slot count;
indices that cannot fit its host-word enumeration return a typed limit error.
Neither this dispatcher nor the generic solver promises unbounded practical
reductions. Slots include propagators and irreducible numerators; exported
vector-key rules do not inherit the old 12-slot restriction.

Independent bounds remain: compact native powers are -64 through 63, legacy
packed ordering has a 34-coordinate ceiling, and artifact loading has separate
configurable resource budgets. None of these is a 12-slot rule-format bound.
