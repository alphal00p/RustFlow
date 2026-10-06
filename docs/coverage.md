# Implementation coverage

The [Python API and notebook status](python-notebook-status.md) records the
completed publication-wheel native gg→hg acceptance: all 16 fresh boundary
configurations, 4,360 transport coefficients, eight W/Z form factors and three
coherent observables, followed by independent 40-digit regeneration,
interruption/restart and nearby-point reuse. The calculation uses certified
supplied parent equations; independently deriving all parent equations with
RustRed is not established. External reference uncertainty caps the EW-square
comparison at 19 relative digits. The broader AMFlow/DiffExp parity goal remains
open; the milestones below retain their original scope and source provenance.

The latest [full release coverage](../reports/validation/2026-10-06-constraints-roots-cuts-integration/summary.json)
passes **651 top-level unit/integration tests across 103 targets**, one doctest
and six separate cache child checks, with 15 tests explicitly opt-in. Formatting
and strict all-target Python/stub-generation Clippy pass. The first run stopped
at an inherited test assertion expecting cache schema 5 instead of 7; after a
test-only repair, that target and all remaining targets passed. All 278 recorded
inputs were verified, with no production-source change between those runs.
This integrates rational endpoint constraints, recursive partial cut placements,
dense registered-root Python transport and balanced native series construction;
the newer constraint/slot Python bindings and immutable algebraic preparation
have separate gates. The frozen community publication wheel is unchanged.

The subsequent [Python interface integration](../reports/validation/2026-10-06-python-interface-integration/summary.json)
passes 80 targeted checks combining explicit physical deformation slots and
exact rational endpoint constraints. A further 26 overlapping checks validate
the corrected community re-export list, all affected Python adapters and native
constraints. Both gates pass formatting and strict all-target Python/stub Clippy.
The export regression executes the actual community wrapper statements and then
uses the exported classes for numerical and cache tests. It retains the initial
harness quoting error and the corrected reproducer's missing-export failure.
These are registered-module checks, not a rebuilt community publication wheel.

The latest isolated milestones add:

- [Exact correlated rational endpoint constraints](../reports/validation/2026-10-06-constrained-singular-endpoints/validation.json):
  118 top-level checks, six cache child checks, formatting and strict Clippy.
  Physical relations are explicit Gaussian-rational assertions; numerical small
  values are never promoted to exact constraints. Constrained root sheets,
  non-Gaussian coefficient fields and dimensional-sector selection remain open.
- [Recursive partial cut mass placements](../reports/validation/2026-10-06-partial-cut-boundaries/validation.json):
  32 native checks and a nonfactorized three-loop comparison with original
  AMFlow. Complete positive-energy final states and uniform virtual signs remain
  required; massive N-body leaves and interior phase-space poles remain open.
- [General dense registered-root Python transport](../reports/validation/2026-10-06-python-algebraic-transport/summary.json):
  43 release tests covering coupled analytic systems, both sheets, prescribed
  paths, endpoints, cancellation, restart and native generated stubs.
- [Balanced native series construction](../reports/performance/2026-10-06-native-series-construction/README.md):
  exact old/new output agreement in matched 48/61-master timings and a completed
  separate full gg→hg transport acceptance. This improves existing native series
  assembly without a Symbolica dependency patch.
- [Immutable exact algebraic preparation](../reports/performance/2026-10-06-exact-algebraic-preparation/report.json):
  120 affected integration tests and a full 16-configuration gg→hg transport
  acceptance pass. Matched checked transports improve by 1.12–1.17× across the
  tested 48/61-master first and nearby destinations, with identical non-timing
  outputs. Numerical coefficients, poles and residuals are still compiled at
  each requested precision; no supplied boundary or rounded coefficient is
  cached in the exact preparation object.

These isolated counts overlap and must not be added together. Their reports
retain exact revisions; none replaces the frozen notebook scientific run.

The [preparation/interface integration](../reports/validation/2026-10-06-prepared-interface-integration/summary.json)
subsequently passes **154 release tests across 24 targets**, formatting and strict
all-target Python/stub Clippy. All 283 recorded build/test inputs match the tested
commit before and after validation. The gate combines immutable algebraic
preparation with rational endpoint constraints, public Python exports, explicit
deformation slots, cut interfaces, prescribed continuation and native amplitude
assembly. One full crossed gg→hg scientific test is opt-in here; the isolated
preparation milestone above records its separate successful execution.

The [combined endpoint/cut/arithmetic gate](../reports/validation/2026-10-06-endpoint-cut-arithmetic-integration/summary.json)
passes **328 top-level release tests across 29 targets**, plus six successful
child-process checks, formatting and strict Python/stub-generation Clippy. Six
explicitly opt-in tests are ignored. All 260 recorded source/input files match
commit `ba16032` before and after validation. Supplied finite singular endpoints
now integrate with regular physical transport and a separate terminal cache;
native HEPKit cut diagrams share the Rust, Python and CLI projection adapters.
Bounded integer-polynomial residual arithmetic is available as an opt-in with
the existing ball path as fallback. At that revision, constrained endpoint
sectors, dimensional sector selection and partial cut placements remained open. This targeted gate
does not rebuild the community publication wheel or repeat its long scientific
acceptance.

The [Python prescribed-contour milestone](../reports/validation/2026-10-06-python-prescribed-transport/summary.json)
passes **58 release tests across seven targets**, formatting and strict
Python/stub-generation Clippy. `ContinuationPrescription` connects descriptive
Python transport and endpoint methods to the existing native planner. Native
exact-coordinate comparison is shared by source selection, Python admission and
cache hits; a compatible hit needs no new homotopy admission but still respects
the source policy and root germ. Both analytic logarithm sides, canonical root
sheets, typed failures, cancellation, stubs and terminal restart are checked.
The report preserves failed preliminary assertions and the exact-hit regression;
none is counted as passing evidence before the correction. This is embedded
Python validation, separate from the frozen community extension.
The subsequent [main integration gate](../reports/validation/2026-10-06-python-contour-integration/summary.json)
passes **62 release tests across nine targets**, including the cut Python API
and residual-arithmetic regressions, with formatting and strict Clippy. All 262
source/input hashes match the tested commit and remained unchanged during the
run. The longer scientific publication acceptance remains a separate revision.

The [combined Padé/cut release](../reports/validation/2026-10-06-pade-cut-integration/summary.json)
passes **574 release unit/integration tests** and one doctest, with 13 explicit
opt-in tests ignored, formatting and strict Python/stub-generation all-target
Clippy. All 248 recorded source/input files match the tested commit and remained
unchanged throughout the run. Optional [rational/epsilon Padé trials](pade-transport.md)
preserve source-domain checks and use the same accepted representation for saved
segments. [Mixed cut/virtual boundaries](cuts.md) now use shared region expansion,
native tensor projection, recursive vacuum evaluation and phase-space terminals.
This covers complete positive-energy final states with all uncut physical lines
deformed and a uniform virtual prescription; partial placements, massive N-body
leaves and interior phase-space poles remain unsupported. Neither the full
gg→hg publication acceptance nor the other opt-in scientific gates was rerun on
this source revision.

Release totals count each Cargo test binary once; six successful native-cache
child-process checks are recorded separately in each archive. Earlier reported
totals of 580 and 547 included those child runs; their corrected top-level counts
are 574 and 541. Raw logs and all passing outcomes are unchanged.

The [cache/arithmetic release](../reports/validation/2026-10-06-cache-arithmetic-release/summary.json)
passes **541 release unit/integration tests** and one doctest, with 13 explicit
opt-in tests ignored, formatting and strict Python/stub-generation all-target
Clippy. It integrates snapshot-local native Atom export reuse and Symbolica's
explicitly rounded Float operations. The numerical sources were verified
unchanged after the full run. The fresh
[ordinary transport comparison](../reports/performance/2026-10-06-ordinary-transport.json)
passes all three original-solver and independent refinement checks; ordinary
native transport remains slower than the original. This source revision is
separate from the completed publication-wheel gg→hg acceptance above.

The [combined linear/endpoint integration](../reports/validation/2026-10-06-endpoint-linear-integration/summary.json)
passes **74 targeted release tests**, formatting and strict Python/stub-generation
all-target Clippy. Recursive Euclidean FT now admits tested mixed linear families
through the public evaluation and projection interfaces. Registered-root systems
can use the shared rational Frobenius engine for local singular endpoints,
including infinity, resonant logarithms and common branch winding. These endpoint
expansions are standalone Rust functionality: they do not infer dimensional
sectors after epsilon expansion, establish global accuracy, or automatically
insert singular endpoints into the physical cache. The later cut release above
extends the separate automatic boundary workflow; upstream/performance parity
remains open.

The massless phase-space milestone passes **416 release unit/integration tests**,
one doctest, strict all-target Clippy and formatting. Nine long scientific tests
remain opt-in. Its [release report](../reports/validation/2026-10-05-massless-phase-space-release.json)
records the tested source and immutable library. The [N-body terminal](cuts.md)
uses one dimensional convolution formula, native HEPKit normalization, Symbolica
routing determinants and cut-aware reductions. Two- and three-loop finite-epsilon
values match independent automatic AMF discontinuities; a two-loop Laurent fit
matches the native three-body measure. General mixed-cut recursion remains open.

The [whole-segment verification milestone](../reports/validation/2026-10-05-whole-segment-residual-release.json)
passes **440 release unit/integration tests**, one doctest, strict all-target
Clippy and formatting, with nine long scientific tests still opt-in. Ordinary
rational, epsilon-hierarchy and registered-root charts retain all source
polynomial degrees when estimating differential defects over a trial disk.
Both step strategies now return the correct -999 in the
[previously false verified sparse-polynomial case](../reports/validation/2026-10-05-residual-alias-counterexample.json),
and the root representation has its own regression. Missing or inconclusive
majorants reject a trial. The optional bracket strategy charges every tested
candidate and commits only the selected state. Exact cache hits preserve source
provenance, and a captured erroneous old snapshot is rejected by source-fingerprint
compatibility before its values are admitted.

The newer [conditioned-transport release](../reports/validation/2026-10-05-conditioned-transport-release.json)
passes **469 release unit/integration tests**, one doctest, strict all-target
Clippy and formatting; 10 long scientific regressions remain opt-in. It adds
regular boundary-centered Möbius charts and directed rational/epsilon residual
enclosures against exact Gaussian-rational source coefficients. Symbolica owns
the polynomial and ball arithmetic. Registered-root residual assembly retains
its separately documented MPFR-estimate scope.

[Arithmetic conditioning checks](arithmetic-conditioning.md) now catch the
[coordinate/coefficient-cancellation failure](../reports/validation/2026-10-05-coordinate-roundoff-counterexample.json)
and request fresh higher-precision sources within a finite retry budget. Relative
componentwise perturbations preserve coordinate scaling and exact zero tails;
the original complex-mass profiles still pass. Values, physical defects, branch
state and saved segments use the same declared endpoint, including high-storage-
precision input coordinates. A stalled rounded step requests precision, and
nonfinite saved-interval coordinates are rejected. Cache evidence retains source
accuracy constraints and can fall back to another boundary after precision
exhaustion. Neither these local checks nor profile agreement proves accumulated
global accuracy; severe source cancellation can still exhaust the strict local
budget. Historical long reference comparisons retain their measured library
revisions and have not all been rerun on this release.

The preceding registered-root epsilon-shearing milestone passes **408 release unit/integration
tests**, one doctest, strict all-target Clippy and formatting. Nine long scientific
tests remain opt-in. The [release report](../reports/validation/2026-10-05-epsilon-algebraic-release.json)
records the exact tested sources and immutable compiled library. Rational and
registered-square-root physical matrices now share the diagonal rescaling adapter,
progressive cache and CLI. Exact root relations are normalized before epsilon
valuation and expansion; original norm domains, branch germs and boundary errors
remain enforced. New tests cover a common transformation across two physical
variables, complex paths, both root sheets, prescribed threshold continuation and
binary restart. General nondiagonal epsilon transformations and
epsilon-dependent radicands remain unsupported.

The preceding rational epsilon-shearing milestone passes **399 release unit/integration tests**, one
doctest, strict all-target Clippy and formatting. Nine long scientific tests are
explicitly opt-in. The [validation report](../reports/validation/2026-10-05-epsilon-shearing-release.json)
records the exact tested sources and dependency hashes. [Diagonal epsilon
rescaling](epsilon-shearing.md) now connects supplied rational systems, automatic
AMF/FT seeds, physical transport, target projection and the CLI through the same
growing cache. A native-IBP-certified massive bubble basis with an epsilon pole
matches an independent direct AMF evaluation after transport; original source
orders, reduction conditions, basis identity and uncertainty remain enforced.
Its registered-root extension is included in the new milestone above.

The complex-domain milestone passes 380 release unit/integration tests and one
doctest, strict all-target Clippy, and formatting. A final six-line correction
preserves native certified-zero precedence in the new mass-sheet guard; all
release library units, all six mass-sheet integration tests, the doctest,
Clippy and formatting pass after that correction. Other integration tests were
not repeated for this final change. The
[release report](../reports/validation/2026-10-05-complex-domains-release.json)
records both source snapshots and the preserved diagnostic failures.

Regular algebraic cache paths now accept exact rational-complex invariants and
masses, with exact root-sheet transitions and checked progressive restart.
Single-mass vacuum recursion admits Gaussian mass squares off the real mass cut
and reuses its massless child across masses and powers. Automatic AMF rejects
recognized complex-mass contours whose principal sheet is uncertified, before
reduction/cache lookup and again when evaluation options change. Compatible
conjugate three-loop cases pass analytic homogeneity and independent refinement
checks. These additions do not establish arbitrary complex-mass homotopy,
algebraic singular-endpoint support or feature parity. The
[comparison tables](benchmark-comparison.md) record current benchmark scope and
side-by-side reference timings.

The preceding physical-family milestone passes 362 release unit/integration
tests and one doctest, strict all-target Clippy, and formatting. Physical-family preparation
now reuses AMF's skipped initial reduction and basis refinement across every
coordinate, with native Gaussian factorization and preserved reduction domains.
Automatic seeds reject an unproved epsilon-dependent family pole bound before
backend calls or cache mutation. The
[release report](../reports/validation/2026-10-05-physical-preparation-options.json)
records the frozen source and library. The preceding
[boundary and contour milestone](../reports/validation/2026-10-05-boundary-contour-release.json)
added topology-independent single-mass vacuum loop peeling and native contour
certificates, including cold-cache derivative-side regressions. Public finite-epsilon tests compare the
new three-loop boundary route with FT under independent refinement. A separate
[public Laurent acceptance](../reports/validation/2026-10-05-single-mass-vacuum-laurent.json)
now passes both three-loop targets through epsilon zero at 20 digits, with a
disjoint sample grid and an independent native FT comparison. The longer
nonplanar physical continuation remains separate. Earlier milestones added prescribed continuation in
the growing cache and CLI, sharper inherited-error estimates and compatible
analytic singular-origin initialization. Separate production scientific runs pass the
full 108-master system, all 16 crossed Higgs+jet supplied-system cases, and the
13-master prescribed cache route. The [integrated report](../reports/validation/2026-10-04-prescribed-integration.json)
records their distinct source/binary provenance and retained failed attempts.
The [coherent full Higgs+jet matrix element](gg-hg-matrix-element.md) also matches
an independent oracle using native HEPKit, with 19 propagated relative digits.
Full feature parity, a uniform 20-digit amplitude result, general cut recursion,
and broad performance parity remain open. The chronological milestones below
retain their original scope and validation counts.

The tool is now named RustFlow; the Rust crate remains `symbolica-amflow`.
The [expanded parity goal](parity-plan.md) includes DiffExp, linear propagators,
cut phase-space integrals, a growing native boundary cache, and HEPKit integration.
These are separate from the completed original-paper acceptance below.

The first RustFlow transport milestone adds exact multivariate path pullbacks,
physical invariant derivatives including Gram variation, direct epsilon-coefficient
integration, independently checked intermediate checkpoints, and a binary
`RustFlowCache` passed by mutable reference. The regular physical-path example
uses 18 steps for its first destination, 1 from a newly cached nearby boundary,
and 0 for an exact repeated hit. See [the native API](rustflow.md).
The pinned live DiffExp polylogarithm example agrees to 30 digits using independently
constructed boundary series. `linear::PreparedLinearFlow` evaluates supported
rank-one linear denominators by quadratic deformation and Frobenius projection;
an HQET example matches its gamma formula to 20 digits under refinement.
`solve_integrals` now dispatches these linear families automatically, including
Laurent reconstruction. [FT recursion for linear families](linear-propagators.md)
also uses the automatic and projected interfaces; independent eikonal directions
and a rank-two two-loop linear form pass 20-digit analytic checks, and the
eikonal case has a live original AMFlow FT reference. AMF's rank-one contour
restriction remains intact. `solve_integral_combinations` applies exact numerator
weights across partial-fraction families at each epsilon sample before fitting,
and includes any additional poles from those weights.

The actual Higgs+jet application grid supplies a closed two-master nonplanar
subsystem. Its Rust transport agrees with the pinned live DiffExp oracle to
20 absolute digits at two destinations, reuses the first destination for the
second, and survives binary cache restart. The input's 24-digit evidence and
inherited errors are retained despite its much longer displayed mantissas.
This is a subsystem acceptance, not a full amplitude calculation.

The [point-bank validation report](../reports/validation/2026-10-04-rustflow-pointbank.json)
records 163 passing release tests, strict all-target Clippy and formatting,
including the live-oracle Higgs+jet subsystem and uncertainty regressions.
The original 20-digit two-loop acceptance is retained as prior evidence; this
milestone did not repeat its hour-long full Laurent run. Full DiffExp benchmarks,
general cut boundaries, and full Higgs+jet amplitude evaluation remain incomplete.

The native HEPKit milestone introduced the DOT/model CLI and delegated tensor
projection to HEPKit. The subsequent cut/boundary milestone adds native Linnet
connectivity and strongly connected components, reverse-unitarity IBP metadata,
positive-energy two-body phase space, and partial Frobenius coefficient matching.
The full [equal-mass banana system](banana-equal.md) reaches both reference points
from its analytic infinity boundary through ε⁴. Native graph families can now
prepare common physical equations and generate verified AMF seeds directly into
the growing cache. Exact reduction assumptions and original matrix domains are
retained even when a stationary coordinate cancels a pole from the pullback.
The [cut/boundary validation report](../reports/validation/2026-10-04-rustflow-cuts-banana.json)
records 206 passing release tests plus one documentation test, strict Clippy and
formatting. General mixed real/virtual cut boundaries remain incomplete.

The [branch/target milestone](../reports/validation/2026-10-04-rustflow-branches-targets.json)
passes 228 release unit/integration tests plus one documentation test, strict
all-target Clippy and formatting. It adds [regular square-root transport](algebraic.md)
with winding and per-call sheet tests, uncertainty-aware target projection from
cached masters, and faster [33,000-point cache selection](cache-selection.md).
The complete [unequal-mass banana](banana-unequal.md) checks all 75 coefficients
against original DiffExp, independent routes, increased precision/order, and
analytic tadpole products. These tests do not establish full parity.

The [multivariate algebraic milestone](../reports/validation/2026-10-04-algebraic-kinematic.json)
passes 239 release unit/integration tests plus one documentation test, Clippy and
formatting. It adds exact physical dlog pullback, faster root charts and shared
native reduction-graph expansion. The complete
[13-master planar one-loop five-point benchmark](fivepoint-planar.md) then passes
all 65 coefficient comparisons against live original DiffExp and ancillary data,
with three independently evaluated precision/order profiles. Its explicit
production regression took 630.50 seconds; it is opt-in for routine test runs.

The [algebraic cache and contour release gate](../reports/validation/2026-10-04-algebraic-cache-release.json)
passes 264 unit/integration tests and one documentation test, with strict Clippy
and formatting. The full five-point regression and a native root-isolation
performance reproducer remain explicit opt-in tests.

The [canonical-system release gate](../reports/validation/2026-10-04-canonical-release.json)
passes 283 unit/integration tests and one documentation test, strict Clippy and
formatting. [Canonical dlog systems](canonical-systems.md) now restrict their
letters to each path before assembling the connection and use the same growing
[physical cache](canonical-cache.md), whose current binary schema is 5. The CLI accepts
registered roots and explicit germs. Exact polynomial grouping and shared source
guards preserve removable letter, radicand and coordinate holes. The
[full 13-master cache experiment](full13-canonical-cache.md) checks every epsilon
coefficient, inherited accuracy, nearby reuse and restart against refined native
transport and original DiffExp. A [connected two-loop sector](fivepoint-zmz-sector.md)
adds 13 coupled/dependent masters with a live oracle and independent refinement;
it does not complete the 75-master parent family. The two longer scientific
regressions also [passed on the committed build](../reports/validation/2026-10-04-canonical-fivepoint.json),
separately from this routine gate (668.28 s one-loop; 291.68 s two-loop sector).

The [cache retry and MPL release gate](../reports/validation/2026-10-04-fallback-mpl-release.json)
passes 295 unit/integration tests and one documentation test, Clippy and formatting;
its report distinguishes concurrent RustRed source changes from the tested build.
[Accuracy-aware source retries](accuracy-fallback.md) try another compatible
boundary when the nearest source cannot meet accuracy. The CLI accepts canonical
letters and retry limits, and ordinary ODEs accept exact complex literals.
All [MPL notebook profiles](diffexp-mpl.md) now pass against original DiffExp and
independent precision/order refinement, including the weight-20 example.
The [full Higgs+jet supplied systems](gg-hg-coverage-inventory.md) pass a nearby
physical-point check of all 545 coefficients across 48 planar and 61 nonplanar
masters. The [crossed supplied-system check](../reports/validation/2026-10-04-gg-hg-crossed-form-factors.json)
now covers all 4,360 W/Z coefficients and eight form factors. Propagated
form-factor bounds support 19–20 relative digits in those separate W/Z contexts.
A [common physical-point comparison](../reports/validation/2026-10-04-gg-hg-coherent-form-factors.json)
now checks both masses at the same exact kinematics; its form-factor bounds support
18–20 relative digits. Native matrix-element contraction and uniform 20-digit
form-factor accuracy remain incomplete.

The [updated RustRed bridge gate](../reports/validation/2026-10-04-rustred-bridge-release.json)
passes 297 unit/integration tests plus one documentation test, strict Clippy and
formatting. Both reduction backends now receive explicit cut semantics.
Continuation budget exhaustion reports a typed limit error with rejection
counters. The report distinguishes the concurrent dependency revisions checked
by release tests and Clippy.

The current [prescribed-contour interface](prescribed-contours.md) accepts
explicit polynomial prescriptions. It does not infer them from a process.
Registered roots now participate in the same [progressive cache](algebraic-cache.md),
with explicit local sheet choices, inherited uncertainty and binary restart.
The separate [prescribed cached transport](prescribed-cache.md) entrypoint accepts
exact affine paths with polynomial prescriptions and an explicit homotopy-admission
policy. It retains actual root germs at verified physical checkpoints, excludes
complex detours from the physical bank, and uses schema 5 for continuation identity.
For connections with an exactly absent epsilon-zero term, inherited uncertainty
uses the finite Dyson series of the retained epsilon hierarchy. The [full prescribed-cache trial](../reports/validation/2026-10-04-prescribed-cache.json)
passes all 65 coefficients of the 13-master five-point system at a 20-digit target
and retains 427 physical bank entries. This includes independent transport
refinement and unchanged source-error evidence. General algebraic singular
endpoints and the full Higgs+jet amplitude remain incomplete.

The [108-master nonplanar benchmark](diffexp-nonplanar108.md) checks all 540
coefficients against original DiffExp and an independently refined native run.
Its initializer handles the compatible analytic, log-free sector of a pure-epsilon
regular singular origin whose registered roots remain nonzero; it does not supply
a general Puiseux or logarithmic algebraic endpoint solver.
The original full [74-master MZZ](../reports/diffexp/fivepoint-mzz-original.json)
and [86-master ZZZ](../reports/diffexp/fivepoint-zzz-original.json) PH1→PH6
benchmarks pass all 370 and 430 coefficient checks at their recorded 15-digit
settings. The [full native MZZ comparison](../reports/diffexp/fivepoint-mzz-validation.json)
also passes all 370 coefficients at 20 absolute digits using independent 60/56
and 70/64 working-precision/order profiles, with maximum refinement difference
1.9624e−48. This is supplied-boundary differential transport; its low-level
refinement estimate does not certify propagation of arbitrary boundary errors.
The [full 75-master ZMZ comparison](../reports/diffexp/fivepoint-zmz-validation.json)
now passes all 375 coefficients at 20 absolute digits, including independent
60/56 and 70/64 precision/order refinement. The maximum refinement difference
is 5.5571e−47 and all three endpoint root sheets agree. The first profile took
9308.38 seconds cumulatively, including capped/discarded work; the independent
high profile took 12157.19 seconds as a fresh process. These results use the
immutable scientific library recorded in the report; its new opt-in current-library
regression has not been rerun. The
[full 86-master ZZZ comparison](../reports/diffexp/fivepoint-zzz-validation.json)
now passes all 430 coefficients at 20 absolute digits, including independent
60/56 and 70/64 precision/order refinement and all three root sheets. The first
profile took 10122.59 seconds including capped/resumed work; fresh high transport
took 14042.05 seconds (14052.49 in-program total). Maximum refinement difference
is 1.5174e−46. The complete output was independently checked after the original
monitor ended and its unchanged child was adopted under the original cap. The
report preserves the unavailable exit status and roughly 345-second RSS polling
gap. This is supplied-boundary acceptance on the named immutable library. The
[separate 128-digit PH1→PH2 original profile](../reports/diffexp/fivepoint-zzz-128-original-failure.json)
ended with a native Mathematica memory failure after 1811.93 seconds and 36 of
108 segments under its declared 16 GiB per-process address-space cap. Neither
watchdog fired and no endpoint result was produced. The [native 128-digit
baseline](../reports/diffexp/fivepoint-zzz-128-native-launch.json) is now running
at 280 working digits/order 230 under fixed four-hour/16 GiB limits; its separate
320-digit/order-280 refinement remains unrun. Finite source data carry about
132 absolute accuracy digits. The staged MZZ offline regression has not yet been rerun against the
latest production library.

The original paper's mandatory numerical acceptance gate passes: all four
two-loop targets at s=30, t=-10/3, m²=1 have automatically generated Laurent
coefficients from ε⁻⁴ through ε⁰ with 20 verified decimal digits and the physical
+i0 prescription. Each coefficient agrees with the pinned upstream reference
within its recorded component precision. The 27-integral differential system
closes, with seven leading boundary coefficients generated recursively; closure
establishes a spanning basis, without certifying independence.

The final fit uses 31 exact samples ε=j/6200 and is checked against 27 samples
ε=j/2700. Working precision increases from 60 to 80 decimal digits and outer
series order from 80 to 96. All four results report one refinement and 20 verified
digits. Full coefficients, observed fit changes, exact grids, reference precision,
and source, binary, and log hashes are preserved in the
[two-loop acceptance report](../reports/validation/2026-10-04-paper-two-loop-acceptance.json).
The observed changes are consistency estimates, not rigorous error bounds or a
claim of accuracy beyond 20 digits.

The acceptance run exited successfully after 3553.007 seconds. It started with
exact preparation cached from the earlier one-sample run and used four additional
workers to prepare exact systems for the second grid. All numerical boundary
construction, continuation, and fitting ran in the main four-worker evaluation;
the helper supplied no numerical values. This timing includes a warm cache and
concurrent preparation and should not be interpreted as a cold four-worker
benchmark.

The earlier fresh evaluation of all four targets at ε=1/2700 completed in
368.906 seconds and passed a separate 60→80-digit, order 80→112 stability check.
Its first target also agrees with the earlier 38-integral-basis calculation.
Full sample values and the cross-basis comparison remain in the
[four-target sample report](../reports/validation/2026-10-04-paper-four-samples.json).

| Upstream operation | Rust implementation | Present limitation |
|---|---|---|
| Family setup, exact IBP | `family`, `reduction` | RustRed finite exact search; compiled arity registry (default 1–16) |
| Auxiliary differential equations | `differential_system` | Stable derivative closure required; residuals are not certified independent masters |
| SkipReduction | `differential_system_skip_initial` | Retains redundant target directions |
| RefineBasis | `refine` | Bounded D-factorizing swaps, reports if complete |
| Block structure | `DifferentialSystem::blocks` | SCC eigenvalue and recurrence LU blocks; full connection matrix |
| Ordinary propagation | `ode` | Pole-bounded Taylor steps with tail checks |
| Singular endpoints/infinity | `frobenius` | Rational indicial roots; integer shearing, Jordan-chain projector refinement, resonant logarithms |
| Auxiliary-mass placement | `MassMode`, `IntegralFamily::deform` | Intrinsic-mass groups, propagators, branches, loops, all lines, explicit positions; ordinary quadratic lines |
| Region enumeration/expansion | `regions` | Ordinary rank-one quadratic propagators; bounded routing budget |
| Tensor reduction | HEPKit/Spenso via `tensor` | Repeated single hard vector: rank 32; general rank 20 with native projector/output budgets |
| Massless bubble subloops | `bubble`, native adapter | Tensor rank at most eight, sum of bubble powers at most 32; transfer square must be external or an active denominator |
| Partial fractions/ISP completion | `integrand` | Affine scalar-product denominators, bounded term budget |
| Automatic boundary recursion | `recursive` | Depth 32; adaptive boundary order at most 16 in inverse mass; half-power numerator series |
| Feynman trick | `ft` | Recursive parameter DEs, regular rational reference points, Gaussian tensor terminals; Euclidean sectors only |
| Epsilon reconstruction | `epsilon`, `engine` | Conservative -2L pole bound; independent sample/precision refinement |
| Cache restart | `cache` | Reduction, native RHS-search checkpoints, and prepared AMF/FT systems; numerical memo is in-memory |

Broader multiloop and complex-kinematics coverage remains limited; the 316-integral
three-loop AMFlow 2.0 example is a later performance benchmark. Pole-disk contour planning is
implemented and tested on both logarithm branches. Exact complex-mass tadpoles
and a massless bubble with the exact external invariant `s = -2 + i` are tested
against analytic gamma formulas. The bubble exercises `KinematicPoint`
substitution and sampled IBP reduction, with 20-digit agreement under working
precision 60→80 digits and series order 80→112. This remains limited coverage
rather than general multiloop complex-kinematics validation.

A backend search budget must never be interpreted as a proof that a residual is
a master. The implementation checks repeated reduction and differential closure
and requires sufficient consistent boundary constraints. It returns typed
errors for incomplete closure, missing table rules, cycles, exhausted budgets,
unsupported normalization, nonfinite arithmetic, or failure of precision
refinement. Native RustRed calls are checked for cancellation before and after
the call; they cannot currently be interrupted from within this adapter.

## Numerical checks

The tests cover physical-region massless and massive bubbles, an on-shell
massless triangle, Euclidean massive bubbles, a non-vacuum single-mass sunrise compared against the
pinned upstream result with independent precision refinement, independent top
sectors, single-mass and two-mass two-loop sunsets with automatically generated region
boundaries, a three-loop single-mass banana vacuum with recursive FT boundaries, factorized
tadpoles, equivalent momentum routings and denominator permutations, linear
numerator boundaries, recursive FT sunsets, supplied-boundary analytic ODEs,
resonant logarithms, successive nilpotent balances, and Laurent fitting. Full
Laurent accuracy is tested with increased
precision/order and a second sample grid. For these automatic AMF examples, reference fixture digits are never used as
solver boundary data. DiffExp supplied-system regressions explicitly identify
their upstream boundary inputs separately from endpoint comparison values.

The test suite requires no Mathematica installation. Upstream regeneration is
optional and would require Mathematica and an upstream-supported reducer.

FT checks include public `solve_integrals` Laurent reconstruction for a Euclidean
massless bubble, analytic gamma-function comparison, and a fresh evaluator
restarting from the symbolic system cache with reducer calls disabled. Changing
cache-relevant options requires a new preparation. Numerical memo keys include
requested accuracy: tightening the tolerance at unchanged working precision and
series order triggers fresh continuation, while an identical request reuses the
value.

The final [release-check report](../reports/validation/2026-10-04-release-checks.json)
records 124 passing release tests, including doctests, strict all-target Clippy,
and formatting checks. It distinguishes the archived AMF acceptance source from
the later FT memoization fix; that fix leaves the AMF implementation unchanged.

## Runtime limits and exceptional samples

Both native adapters follow RustRed’s compiled arity registry, defaulting to
1 through 16 scalar-product slots. Its default
adapter search uses depth 2 and at most 4096 concrete integrals; differential
closure uses at most 12 reduction rounds. RustRed's compact integral indices
range from -64 through 63, including generated seed indices. Exact rational epsilon specialization
is enabled by default for high-level multiloop AMF evaluation. It is covered by
a test comparing the resulting dimensional sectors with an analytic sunset.
Enabling basis refinement retains symbolic epsilon before factorization, even
when sampled reduction is otherwise enabled. Refined `new_at_epsilon` flows
therefore remain reusable at other epsilon values.
Recursive multiloop boundary families follow the sampled-reduction option too;
their memo keys include the exact epsilon sample. Boundary calls at exceptional
samples where dimensional exponent classes collide retain symbolic epsilon.

Nilpotent normalization first searches small subsets of Jordan-chain
projectors, then applies exact Moser projector refinement. Large chain counts
skip subset enumeration. At most 32 balances are performed, and a balance must
decrease the pole-order/rank pair. Unsupported indicial roots and an unsuccessful
normalization return explicit errors. Automatic AMF
rejects epsilon samples that create extra resonances between symbolic exponent
classes. This guards against an observed endpoint-sector ambiguity at D=1.

Native factorized replay is the default coefficient backend. It preserves source
and denominator nonzero conditions, performs exact substitution before searching
the next RHS frontier, and can checkpoint completed work within a search round. A larger
`max_targets` budget can resume that checkpoint. Finite-depth search can still
produce expensive intermediate rules; a checkpoint or a completed reduction
round does not imply that the differential basis has closed.

The opt-in `parametric_rules` backend option reuses guarded formulas with symbolic
positive powers and fixed numerator indices. Each specialization checks its
domain, preserves nonzero conditions, and must decrease the native integral
order. Uncovered cases fall back to concrete searches, and all RHS integrals
still require closure. See the [reduction audit](reduction-search-audit.md) for
the paper-example search bottleneck and the scope of prototype measurements.

The separate opt-in `symmetry_rules` option proposes bounded loop translations
and external reflections from rational branch displacements and Gram entries.
RustRed verifies each map, including masses and the auxiliary mass, before
compiling finite numerator transport. Only complete, strictly descending rules
are used; unsupported maps and exhausted discovery/transport limits fall back
to IBPs. Enabled runs use distinct backend and checkpoint identities.
`SymmetryReduction` events expose discovery and application counts.

For large native frontiers, coefficient expansion can cost more than the extra
IBP searches needed by a conservative dependency traversal. The
`max_exact_frontier` option controls this tradeoff (512 by default). All reported
final reductions still undergo exact substitution and a check that every
surviving residual was searched. Taylor transport checks the differential-equation
defect at step endpoints and midpoints as well as its final series terms, and
its truncation tolerance tightens when guard precision increases.

The native adapter also reduces eligible massless two-point subloops by exact
beta-function identities before requesting IBPs. These rational identities use
the same integral order as the sector solver, and each proposed rule is checked
to strictly decrease that order. Unsupported routings, complex coefficients,
exceptional dimensions, and nondecreasing rules fall back to IBPs. Set
`bubble_subloops: false` to disable this optimization. Checkpoint keys distinguish
the changed ordering; legacy partial checkpoints are reused only if no already
searched sector changes ordering.

`max_sector_batch` bounds the number of targets sharing one native exact
elimination (unbounded by default). Smaller batches can reduce intermediate
algebra at the cost of repeated seed searches. They retain the same strict
integral order and final exact substitution. `SectorReduced` progress events
report seeds, generated rows, exact-trace size, and total/exact milliseconds.
The acceptance example accepts `--case-batch N` for profiling this tradeoff.

Native search checkpoints use a versioned binary codec with authenticated
payloads and atomic replacement. Legacy JSON checkpoints remain readable;
new writes use `.bin` files. An incompatible or corrupt binary checkpoint is
rejected rather than silently replaced by an older JSON snapshot. Completed
reduction and differential-system caches retain their JSON containers.
With native checkpoints enabled, completed batches are saved every 60 seconds
by default, at search completion, and when cancellation is observed. Configure
`checkpoint_interval` to change this cadence. Only completed searches enter the
saved visited set; unprocessed targets and new right-hand-side candidates remain
pending. Resuming a partial round may search a conservative superset before the
next exact substitution.
Short intermediate rounds share that interval; rewriting a large checkpoint
after every small frontier change can otherwise dominate the actual search.
The first exact substitution at each active-line level also forces a snapshot,
so a long coefficient replay cannot lose the preceding searches merely because
they took less than the periodic interval. `SubstitutionPlan` reports the actual
retained nodes, terminal count, sector threshold, and substitution direction.

Successful native calls also save a family identity bank. A new target batch
can reuse its reachable raw equations when family, exact epsilon, search depth,
ordering options, and dependency versions agree. Every surviving residual leaf
is searched again in that new batch; only equations are reused as completed
work. Target-specific checkpoints take priority. Banks preserve nonzero
conditions and use the same authenticated atomic codec. Concurrent writers may
replace a bank with a smaller valid subset, affecting performance only.

Increasing native search depth by one can also resume a shallower checkpoint.
Exact identities are retained, but every surviving residual is searched again
at the new depth. This matters even for a derivative-closed system: a redundant
basis can have spurious indicial modes. The physical paper subsector with lines
1, 3, 5, and 7 needs depth three in the current backend; depth two is rejected
by dimensional-sector validation.

`native_workers` enables bounded concurrent native batches (default one,
maximum 64 per reduction call). Results are merged in the same order as serial
search, and checkpoints retain only merged work. In-flight work is bounded by
one wave of batches. Account for both this setting and `FlowOptions::workers`
when budgeting total CPU and memory; the acceptance example exposes it as
`--native-workers N`.

Large reduction tables reuse Symbolica state headers during serialization and
reuse imported state maps during loading; entries remain readable in the
existing format. Conservative frontier traversal borrows graph keys instead
of allocating a key for every edge. Final exact substitution groups incoming
rational coefficients in balanced sums to reduce intermediate polynomial
growth; this changes evaluation order, not the reduction identities.
When at most `RustRedBackend::max_backward_frontier` terminal integrals remain
(128 by default), substitution runs from the leaves toward the targets so
apparent poles can cancel within each identity. Child expansions are released
after their last parent uses them. Wider frontiers use forward accumulation
to avoid storing a large residual map at every graph node. Raising this limit
may reduce polynomial growth and runtime when local cancellations are strong,
at the cost of retaining more intermediate coefficient maps. Lowering it can
reduce memory use; zero selects forward accumulation for every nonempty
terminal set. Both strategies perform exact algebra and use the same cycle
checks. The three paper examples expose `--max-backward-frontier N`.
Nondefault values distinguish complete reduction/system cache identities;
native partial checkpoints remain reusable because the saved equations and
integral order are unchanged. Performance depends on the particular reduction
graph and cancellation pattern.
Small-basis substitution shares dependencies across a batch of targets, and
completed frontier expansions are reused when producing the final table.
`DifferentialClosure` events report basis size and remaining derivative targets.

Successive native searches at the same sector threshold reuse each target's
exact weighted frontier, including searched residuals. Appended identities are
applied to those coefficients; changed existing rules invalidate the cache.
Entering a lower sector also clears it, and the original targets are expanded
again to recover discarded lower-sector contributions. The complete dependency
graph is checked for cycles and sector increases before reuse, including paths
that previously cancelled. This cache is in memory only and leaves checkpoint
formats and final closure checks unchanged. Regressions compare incremental and
fresh expansions across partial rounds, exact cancellations, changed parameter
sets, rule replacements, wide frontiers, and cycles. See the
[controlled replay measurements](performance.md#incremental-frontier-replay)
for the measured scope of this optimization.

Before coefficient propagation, the native adapter installs zero identities for
all reachable leaves whose sectors RustRed has already certified as scaleless.
This includes lower sectors that have not reached the search queue yet. The
certificate domain conditions are retained, and full dependency-graph validation
still runs first. The same pruning applies when final substitution resumes from
a completed checkpoint; uncertified sectors remain subject to ordinary search.

Pole finding normalizes each factor and scales its variable using a power-of-two
Fujiwara bound before calling Symbolica's root finder. Newton corrections and
polynomial reconstruction check the returned root set. Regressions include a
degree-12 factor from the paper's boundary reduction and roots at very large and
very small scales; unresolved clusters return a numerical error.
Ordinary Taylor propagation clears denominators exactly per matrix row and
recurs over finite polynomial coefficients, omitting structural zeros. Tests
compare these coefficients against the independent rational-series recurrence
at complex centers. Tail and differential-equation defect checks are unchanged.
The pinned Symbolica `main` revision also checks root corrections and relative
coefficient backward error, fixing the coefficient-scale-dependent convergence
report reproduced in `repros/symbolica-root-scaling`. Scaling and root-set
validation remain in the AMFlow adapter. Cache keys include the manifest and
lockfile fingerprint, so dependency updates invalidate previous prepared systems,
reduction tables, and native search checkpoints.

Vacuum tensor projection groups pairings under permutations of repeated hard
vectors. At rank eight, four copies of each of two hard vectors need a
three-dimensional invariant system instead of the full 105-dimensional pairing
matrix. General distinct vectors retain the full projection. Analytic rank-eight
contractions and mixed-vector contractions test the grouping. Frobenius
recurrences traverse nonzero series-matrix entries without allocating a fresh
numerical zero for every comparison.
Single-vector projection uses the exact angular-moment formula with a memoized
Wick pairing sum, bounded to 100000 states. This extends the boundary tensors
beyond rank eight without constructing the full labelled-pairing matrix.

Boundary orders are selected by the rank of coefficient constraints in the
Frobenius solution. The planner requests all overlapping regions at each selected
power, retains uncomputed regions as unknown, and uses known-zero power and
logarithm constraints. It avoids expanding every component to the first nonzero
term of every fundamental solution. `BoundaryPlan` progress events report the
number of selected region series, coefficient count, and largest half-order.
Boundary matching still checks every computed coefficient, and numerical
refinement remains required to establish accuracy.

Factorized boundaries collect exact integral products before applying the term
budget. A regression from the paper's six-line boundary family at half-order
eight reduces over 10000 monomial contributions to 75 distinct products,
removing a demonstrated boundary-budget failure without discarding terms.

`RecursiveBoundary::evaluate_many(family, targets, epsilon, precision)` evaluates
related boundary powers together. Region construction first collects the factors
needed by all selected coefficients, removes products with a scaleless factor,
and deduplicates requests within identical family definitions. Analytic and
user-supplied terminal values retain priority. Remaining requests share one
`PreparedFlow` per auxiliary-mass placement; vacuum targets choosing different
massive lines remain in separate groups. Prepared flows retain the complete
ordered target set in their cache key. Their target reduction coefficients are
multiplied into the endpoint Frobenius series before selecting the physical
limit, so poles in these coefficients are retained.

This follows the grouping in `BoundaryIntegrals` and `ReduceBoundary`, followed
by recursive evaluation of boundary master lists in `AMFSystemsSetup`, in the
[pinned upstream AMFlow.m](https://gitlab.com/multiloop-pku/amflow/-/blob/26005517a288086c4cb4d1b26d829691bc088485/AMFlow.m).
The batching regression checks preparation/backend-call savings for distinct
powers, differing vacuum mass choices, independent tadpole-product values,
and stability between 60 and 80 working decimal digits. Other checks cover
custom terminal providers, scaleless requests, duplicate targets, and reuse of
canonical numerical memo entries after denominator permutations. A connected
two-mass sunset batch with powers `111`, `211`, and `121` checks recursive
siblings against scalar solves and an independent beta-integral identity,
including a fresh evaluation at higher precision and expansion order.

Batching currently requires identical family definitions; it does not merge
families through a new momentum-routing or denominator-completion transformation.
The numerical memo remains per canonical integral, epsilon sample, and precision.
Recursion limits count nested batches rather than the number of integrals in a
batch, and encountering an active ancestor integral still returns a cycle error.
