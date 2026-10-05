# Implementation and reference comparisons

Snapshot: 2026-10-05. RustFlow is not feature-complete, and broad performance
parity is not established. The interfaces are topology-independent; current
algorithmic coverage and reduction budgets still restrict supported evaluations.
Detailed limitations and chronological validation are in [coverage](coverage.md)
and the [unified-flow design](unified-flow.md).

## Implementation status

| Part of the plan | Implemented | Still outstanding |
|---|---|---|
| Native family and reduction input | HEPKit graphs/kinematics, native tensor owners, exact families, RustRed and scoped table backends, derivative closure, reduction domains | Reducer search/compiled arity limits remain explicit; a target-only table need not cover derivatives or recursive families |
| Shared numerical engine | Rational Taylor transport, rational Frobenius matching with resonances/logarithms, direct epsilon hierarchies, algebraic root charts, whole-segment Taylor defects | Arithmetic cancellation and accumulated error accounting; general algebraic Puiseux/logarithmic singular endpoints and degenerate Gram charts |
| Automatic AMF boundaries | Regions, factorization, recursive boundaries, gamma terminals, generic single-mass loop peeling including complex masses off the mass cut; separate unequal-mass two-body and massless N-body cut terminals | General linear/cut recursion and automatic Higgs+jet initialization |
| AMFlow 2.0 options | Euclidean FT, configurable dimension, skipped initial reduction, basis refinement, symbolic restart | Broader validation/performance; FT physical contours remain outside its supported automatic scope |
| Physical transport and cache | Progressive binary `RustFlowCache`, compatible-source search, checked intermediate points, inherited errors, supplied/AMF seeds, regular complex paths with explicit root germs | General causal-path inference, arbitrary algebraic coordinate constants and non-affine cached routes |
| DiffExp controls | Exact multivariate pullback, canonical dlogs, supplied and compatible analytic-origin boundaries, explicit prescribed detours, exact diagonal epsilon rescaling of rational and registered-root systems, optional checked bracket proposals | General nondiagonal epsilon regularization, Möbius/Padé controls, full singular algebraic boundary matching |
| Higgs+jet amplitude | Supplied full master systems, W/Z crossings, native HEPKit contraction, HEFT square/interference | Uniform 20-digit EW square, fresh automatic AMF boundaries, matched full-amplitude timing |

AMF direct evaluation and AMF-generated seeds for physical transport share the
same continuation owners. No benchmark topology is selected by a hard-coded
family identifier. Successful two- and three-loop examples do not imply that
every loop count, kinematic point or mass pattern is already evaluable.

## AMFlow comparisons

All times are seconds. Rows specify the measured phase. Precision, checking,
initialization and resource settings differ unless the linked report states
otherwise; these columns are not general speedup ratios.
Each report records the library revision actually measured. Historical numerical
runs have not all been repeated after every subsequent feature change.
The first four rows compare the bundled standalone C++ DE solver. Their error
test is scaled: absolute error divided by the larger of one and the compared
value magnitudes must be below `1e-20`.

| Case | RustFlow | Original AMFlow | Numerical result and scope |
|---|---:|---:|---|
| [Analytic logarithmic chain](../reports/performance/2026-10-05-whole-segment-checked.json) | 0.04116 transport; 0.04809 process | 0.001907 transport; 1.72208 process | Current checked release, five-run medians; 20 digits and independent refinement |
| [Upstream 12-master DE](../reports/performance/2026-10-05-whole-segment-checked.json) | 0.57962 transport; 0.69189 process | 0.10217 transport; 4.71912 process | Current checked release, five-run medians; all 12 at 20 digits and independent refinement |
| [Paper-derived 27-master DE](../reports/performance/2026-10-05-whole-segment-checked.json) | 4.34252 transport; 4.54451 process | 0.06195 transport; 14.38687 process | Current checked release, five-run medians; all 27 at 20 digits, deterministic test boundary |
| [Resonant singular matching](../reports/performance/2026-10-05-whole-segment-checked.json) | 0.00734 matching; 0.01435 process | 1.14274 matching; 1.14877 process | Current release, five-run medians; 20 digits, analytic/refinement checks; different internal precision requirements |
| [Automatic massless bubble](../reports/validation/2026-10-04-live-upstream-oracles.json) | 0.1868 cross-check including refinement | 26.7388 `SolveIntegrals` | 20 digits against the original and analytic gamma formula; different measured work |
| [Automatic vacuum sunset, mass squares (1,1,0)](../reports/validation/2026-10-04-live-upstream-oracles.json) | 0.1179 cross-check including refinement | 9.0585 `SolveIntegrals` | 20 digits against the original and analytic gamma formula; different measured work |
| [Automatic connected three-loop single-mass vacuum, two targets at epsilon=1/10](../reports/performance/single-mass-full-workflow.json) | 2.814 process | 175.447 process | Both pass 35 relative digits against 40-digit original data; fresh caches, one worker, different bases/extra-order controls |
| [Required paper example: all four two-loop targets through epsilon zero](../reports/validation/2026-10-04-paper-two-loop-acceptance.json) | 3553.007 complete acceptance | [2428.977 process](../reports/performance/2026-10-05-original-paper-full-workflow.json) | All 20 complex coefficients agree within recorded accuracy at the requested 20 digits; different caches, workers, bases and precision/sample profiles |

The first original four-target run failed during Wolfram child-thread
initialization. An unchanged copied subsystem succeeded after explicitly setting
Wolfram's internal thread counts to one. A subsequent
[30-minute run](../reports/performance/paper-original-censored.json) completed
all 23 boundary systems and eight of fourteen top-system samples before its cap
(1805.707 seconds including cleanup). A separate fresh run with a declared
one-hour cap completed successfully in 2428.977 seconds, including all recursive
boundaries, reduction, sample evaluations and fitting. The original uses one
worker, a 51-integral top basis, 23 boundary systems, 14 samples, 180 working
digits and order 360 plus 50 extra terms. The historical native acceptance uses
a 27-integral basis, warm exact caches, four workers plus four helpers, 31 fit
samples checked against 27 validation samples, and 60/80 working digits with
orders 80/96. This is a complete-workflow comparison with different work and
settings, not a matched speed ratio. Observed coefficient differences near
`1e-38` do not promote the original's configured 20-digit accuracy claim.
The native three-loop **Laurent** AMF/FT
comparison took 120.115/52.336 seconds at different profiles; those are separate
from the finite-epsilon original comparison above.

## DiffExp comparisons

The large systems use supplied boundaries. Their checks validate numerical
transport and reference agreement, not automatic AMF construction of those
boundaries. "Low/high" denotes independent precision/order profiles, not repeated
timing trials. Most rows have one timing observation per profile.

| Case | RustFlow | Original DiffExp | Numerical result and scope |
|---|---:|---:|---|
| [Equal-mass three-loop banana](../reports/diffexp/banana-equal-validation.json) | 2.486 initialization + 0.335 + 0.583 transport | 43.309 initialization + 20.805 + 33.859 transport | All 20 coefficients at each of two endpoints, 20 absolute digits |
| [Unequal-mass three-loop banana](../reports/diffexp/banana-unequal-validation.json) | 16.78 low / 25.63 high; 8.04 alternate, in-program totals | 350.10 alternate process; direct run censored at 1599.51 | All 75 coefficients, 20 absolute digits; path/profile differences |
| [MPL G(1,0,1;4)](../fixtures/diffexp/mpl-101-4.json) | Not separately timed | 0.2845 / 0.2879 transport | 20-digit cross-check; original order refinement and notebook comparison |
| [Trailing-zero MPL](../reports/diffexp/mpl-additional-validation.json) | 0.5121 low / 0.3559 high in-program total | 0.7683 / 0.9261 helper | 20 digits; different precision/order settings |
| [Complex-letter MPL](../reports/diffexp/mpl-additional-validation.json) | 0.08505 low / 0.06067 high in-program total | 0.4433 / 0.5102 helper | 20 digits; different precision/order settings |
| [Weight-20 MPL](../reports/diffexp/mpl-additional-validation.json) | 5.9746 low / 10.1817 high in-program total | 33.6691 helper | 20 digits; different precision/order settings |
| [Full 13-master planar one-loop five-point](../reports/diffexp/fivepoint-planar-validation.json) | 59.55 low transport; 79.65 total | 89.44 transport; 92.03 process | All 65 coefficients, 20 absolute digits, independent refinement |
| [Connected 13-master ZMZ sector](../reports/diffexp/fivepoint-zmz-sector-validation.json) | 66.16 low / 195.98 high transport | 96.99 transport | All 65 coefficients, 20 absolute digits; a sector of the 75-master system |
| [Full 108-master nonplanar five-point](../reports/diffexp/nonplanar108-validation.json) | 233.89 low / 684.63 high transport | 200.73 accurate transport | All 540 coefficients, 20 absolute digits; original fast 44.08-second profile only passes 14 digits |
| [Full 74-master MZZ five-point](../reports/diffexp/fivepoint-mzz-validation.json) | 3524.86 low / 4541.57 high transport | 1147.09 transport | All 370 coefficients, 20 absolute digits, independent refinement |
| [Full 75-master ZMZ five-point](../reports/diffexp/fivepoint-zmz-validation.json) | 9308.38 cumulative low / 12157.19 fresh high process | 1847.95 transport; 1854.56 process | All 375 coefficients, 20 absolute digits, independent refinement |
| [Full 86-master ZZZ five-point](../reports/diffexp/fivepoint-zzz-first-endpoint.json) | 10122.59 cumulative first endpoint | 2453.68 transport | All 430 first comparisons pass 20 absolute digits; independent high profile pending |
| [ZZZ 86-master, separate 128-digit notebook profile](../reports/diffexp/fivepoint-zzz-128-original-failure.json) | [Baseline running](../reports/diffexp/fivepoint-zzz-128-native-launch.json); no endpoint result | Mathematica memory failure after 1811.93 s; 36/108 segments | Acceptance pending; original failed under declared 16 GiB address-space cap, no wall/RSS watchdog censor |

The 75- and 86-master cumulative times include capped predecessors, discarded
work and checkpoint restoration. They are not uninterrupted transport timings.
The higher profiles start independently from the supplied original boundary.
Original MZZ/ZMZ/ZZZ notebook runs request 15 digits; their observed reference
agreement is sufficient for the stated 20-digit comparisons, but does not change
their configured accuracy target. The full 75-master refinement now passes; the full 86-master refinement is
still pending. These scientific runs use the immutable library revisions named
in their reports; the prepared current-library full-system regressions remain
separate, opt-in validations.

A [controlled step-selection experiment](../reports/performance/step-bracket-leg100.json)
on one supplied 75-master segment reduced native transport from 729.713 to
581.216 seconds by using fewer Taylor charts. All 375 coefficients agree within
20 absolute digits and all root sheets match. This is one pair of measurements;
the strategy is now an opt-in control. This historical measurement predates
the new whole-segment residual checks and is not a timing of the corrected
solver. The first four AMFlow rows above now measure the corrected release against
its frozen archive; their extra checks and trial counts are recorded.
Neither change is included in the full-system timings above.

## Higgs+jet reuse and matrix elements

These [nearby-point measurements](../reports/performance/gg-hg-matched-nearby.json)
are medians of three repetitions on the same CPU, with all coefficients compared
at 20 absolute digits and 24-digit source evidence. Native fixed-precision and
checked-cache interfaces do different work; the latter independently refines
transport and propagates source uncertainty.

| Case | RustFlow | Original application/DiffExp | Scope |
|---|---:|---:|---|
| 48 planar masters, nearby W point | 0.976 fixed; 4.533 checked | 1.060 public driver | 240 coefficients |
| 61 nonplanar masters, nearby W point | 2.046 fixed; 11.543 checked | 2.376 public driver | 305 coefficients |
| Exact cached repeat, planar / nonplanar | 0.00190 / 0.01746 | 0.0602 / 0.2018 | Original driver reloads matrices; inner selection is much shorter |
| [Coherent EW squared matrix element](../reports/validation/2026-10-04-gg-hg-squared-matrix-element.json) | No comparable isolated timing | No comparable isolated timing | Native HEPKit and independent original arithmetic agree; propagated supplied-input errors support 19 relative digits |
| [HEFT square and EW interference](../reports/validation/2026-10-05-gg-hg-heft-interference.json) | No comparable isolated timing | No comparable isolated timing | 30-digit arithmetic check; interference has 20 conditional relative digits after propagating EW input errors |

There is no established full-amplitude end-to-end speed ratio. Source accuracy,
AMF initialization, uncertainty checks, serialization, master transport and
amplitude contraction must be included consistently before making that claim.
