# RustRed assessment for native finite-density flow

At RustRed revision `78969aab524b7d6a2eec36f59b01e9af1e04cc04`, the
guarded source API is sufficient for the validated massive two-loop flow.
No RustRed correctness defect has been demonstrated by this work. Practical
four-loop closure is still unproved: the saved E7 pilots encounter bounded
discovery and frontier limits, not a proof that the native algorithm cannot
close their sources.

RustRed is now vendored in this repository. The first local change combines
exact coefficients of repeated `NoApplicableRule` remainders before returning
a recursive reduction. Previously, different reduction paths could emit the
same unresolved label separately, even when their coefficients cancelled.
That was a correct uncollected expression but could give RustFlow unnecessary
provisional labels. The change retains every accumulated nonzero condition and
keeps work-limit, pole, and other failure records separate. It does not change
source discovery, ordering, rule certification, or persistence. Its native
replayed-diamond regressions and controlled E7 comparisons are recorded
separately from the historical pilots below.

RustFlow supplies physical distribution identities, their admitted index
domains, deformation, and retained coefficient conditions. RustRed owns rule
discovery, ordering, application, exact replay, and persisted proof programs.
The native rule program must reconstruct every target and every derivative
of the proposed spanning basis before RustFlow calls it closed. An unresolved
remainder is never sufficient evidence for a master integral.

The massive sunset closes with an 11-member double-cut basis after RustFlow
requests constant frontier sectors and refines encountered finite guard
intervals. The complete amplitude, including its nonzero vacuum, both single
cuts, and double cut, passes 40 independent fixed-dimensional reference
comparisons and 30 precision/order/path refinement comparisons. These are
two targets at epsilon 4/5, not generic four-loop acceptance. See the
[saved comparison](../reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-frontier-reference-comparison.json)
and [guard refinement description](finite-density-guard-refinement.md).

| Observed limitation | Current owner and evidence | Appropriate next step |
| --- | --- | --- |
| Broad nonpositive index boxes consume native domain budgets | RustFlow chooses the requested boxes; RustRed reports `DomainBudget` explicitly. All retained discovery gaps in the original depth-3 E7 pilots have that reason. Exact zero faces remove those gaps in the single-cut first round and greatly reduce them in the double-cut first round. | Continue bounded derivative closure and targeted guard refinement before recommending an upstream algorithm change. |
| Constant provisional sectors were not searched | RustFlow originally requested only derivative sectors. Requesting the whole frontier reduced the massive double-cut basis from 64 to 11 while preserving exact replay and numerical results. | Keep the validated generic frontier policy; do not hand-select masters. |
| Extra storage coordinates cost search work | RustFlow freezes unused native coordinates to zero; native guarded construction still enumerates with storage arity N. | Consider a native guarded physical-arity/fixed-tail API as a performance extension. |
| Massless origin and endpoint admission | These require physical contour, regulator-order, and boundary arguments outside a symbolic reducer. | Retain explicit physical admission certificates and finite-basis audits in RustFlow. Algebraic closure alone does not admit numerical evaluation. |

The depth-3 E7 single cut has physical arity 16, 84 source rows, 1650 certified
native rules and 5 applications. Its frontier reaches 147 above the configured
128 limit after 21.611 seconds. The double cut has physical arity 18 in storage
20, 328 source rows, 874 rules and one application; after three rounds its
frontier is 88 and closure remains unresolved (108.605 seconds). Both use
depth 3, 8192 domains, no adaptive guard passes, and a 180-second process cap.
These are algebraic source pilots with no contour or numerical admission.
Their factory does not attach sealed massless-origin lower-contact zero
evidence. The later typed comparison described below tests that stronger source
context independently.
The [owner report](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-final-owner-gates.json)
contains all six pilot attempts, exact commands, resources, and source hashes.

The optional native extension should accept physical arity and fixed storage
coordinates directly in guarded source construction, drawing on the existing
`SourceSystem::new_with_fixed`/active-arity infrastructure. It must bind that
information into ordering, replay, persistence and context identity, reject
tail escapes at every application boundary, and avoid any ordinary scaleless
sector inference. This is a proposed API/performance improvement, not a
required correctness repair. The [capacity audit](finite-density-native-capacity.md)
records the distinction from the ordinary application's capacity registry and
the current symbolic-power/export limits.
At the pinned revision, the relevant native locations are
`crates/rustred-core/src/solver/guarded/model.rs:87`
(`SourceSystem::new`), `solver/source.rs:45` and `:200` (fixed coordinates and
active arity), `solver/guarded/replay.rs:71` and
`solver/guarded/persistence.rs:99` (the current physical-arity-equals-N checks).

A small executable reproducer is
[finite_density_capacity.rs](../tests/finite_density_capacity.rs). It derives
an actual two-factor weighted surface source, embeds only its labels and
guards into capacity four, and compares reduction, replay and decoded-program
application. It also rejects invalid tails, domain widening, dummy coefficient
dependence and cache mismatches. The physical 7-to-12 and 9-to-12 comparisons
add exact equality of closed bases, matrices, weights and conditions plus
numerical agreement; see the
[capacity comparison](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-physical-comparison.json).

For a larger upstream performance reproducer,
[finite_density_generic_closure.rs](../tests/finite_density_generic_closure.rs)
generates the genuine four-loop source corpus from the graph input. The
[archive map](../reports/validation/2026-10-09-finite-density-native-assembly/four-loop-source-pilot-archives.json)
locates lossless `sources.json.gz` and `result.json.gz` files, with hashes of
both original bytes and archives. Every source row, guard, nonzero condition,
native gap and round program is retained. No supplied numerical oracle enters
source generation or reduction.

The opt-in zero-face experiment changes only discovery presentation:
requested ordinary index zero receives the exact face `[0,0]`, and a requested
negative index receives `(-infinity,-1]`; positive boxes are unchanged. Each
new box is a subset of the previous nonpositive box and still contains its
requested point. Other axes remain symbolic, physical admission is intersected
afterward, and no new source equation, zero identity or elimination is added.
The same native final replay/closure audit remains mandatory. The experiment
remains off by default.

The [controlled comparison](../reports/validation/2026-10-09-finite-density-native-assembly/four-loop-zero-face-comparison.json)
uses exactly the same mathematical source data, target list and discovery
budgets. The new executable enables Python, so its dependency-feature digest
changes; that provenance difference is recorded separately from source equality.
The single-cut first round produces 5410 rules and 88 applications, with 17
remaining `UnprovedDescent` domains and no `DomainBudget` gaps. The double-cut
first round produces 5237 rules and 75 applications, with 3 `UnprovedDescent`
and 445 `DomainBudget` gaps. Both remain unresolved: their frontiers, 324 and
197, exceed the unchanged 128 limit. Their closure times are 64.246 and 57.809
seconds. These times compare different completed search work to the three-round
baseline and do not establish a speedup.

The [larger-frontier follow-up](../reports/validation/2026-10-09-finite-density-native-assembly/four-loop-zero-face-expanded-comparison.json)
retains the three-round/depth-3/8192-domain/180-second controls, raises the
frontier cap to 1024 and requested-label cap to 4096, and still uses no adaptive
guard passes. The single cut completes two rounds in 158.025 seconds, applies
7204 rules cumulatively, and then exceeds the frontier cap at 1869. Its final
program has 3191 rules and reports 5919 `DomainBudget`, 15 `SearchExhausted`,
and 3 `UnprovedDescent` gaps. The double cut reaches its first frontier of 197
and starts searching derivatives, but hits the 180-second process cap before
that next discovery completes. Its final rule/gap counts are unavailable; the
first-round exact program is preserved. Neither case establishes closure.

The small `UnprovedDescent` witnesses say that a recentered source proof misses
the requested domain; their finite width-two intervals are candidates for the
existing native guard-refinement policy. This is a failure to certify that
particular candidate/domain, not evidence that no valid reduction exists or
that a different ordering is necessary. The
[extracted witnesses](../reports/validation/2026-10-09-finite-density-native-assembly/four-loop-zero-face-recenter-witnesses.json)
record exact domains and source/program identity. They have not been minimized
to a smaller source corpus or rerun as isolated boxes. A next bounded comparison
can enable the existing native guard refinement while preserving all unresolved
parent gaps and the final exact closure audit. The next physical comparison
must distinguish physical source support from native search coverage rather
than interpreting a legacy formal-source frontier as a backend obstruction.

The [typed origin comparison](../reports/validation/2026-10-09-finite-density-native-assembly/E7-origin-pilot-comparison.json)
now admits E7's sealed dimensional-origin evidence and repeats the original
three-round, frontier-128, depth-3, 8192-domain profile with ordinary zero faces
and adaptive guard refinement disabled. Both public occupied-flow preparations
fail explicitly before numerical transport (test exit 101). The single cut
ends at frontier 147 in 22.168 seconds and the double cut at frontier 88 in
109.7 seconds. All saved requested and frontier arrays in all three rounds are
exactly equal to their formal-source counterparts. The source identity differs
because it now binds the physical origin proof. Exact native programs and
metadata are retained in the
[lossless archives](../reports/validation/2026-10-09-finite-density-native-assembly/E7-origin-pilot-archives.json).
This comparison still does not add independently certified zero domains for
unrestricted virtual polynomial directions. It demonstrates neither a closed
four-loop basis nor a RustRed correctness defect.

### First-residual coverage isolated from whole-family closure

The historical E7 origin-only programs were inspected and then reloaded through native exact replay. For both first residuals, no saved rule domain contains the target, although 20 original Lorentz-source guards admit it. The nearest applicable target patterns are excluded by concrete raised ordinary-index guards (or a separate compact-loop cut guard), rather than coefficient poles or failed descent at application.

Bounded rediscovery using the unchanged historical source corpus finds native rules on exact-point and one-ray boxes for both targets. More usefully, the double-cut ordinary-zero-face box finishes its initial discovery at depth 3 / 8192 domains with 3652 rules and no discovery gaps, and covers the target. The corresponding broad sign box still misses it with 6251 domain-budget gaps. The single-cut broad box does cover its target when isolated from the larger request queue. These are concrete frontend domain-presentation and queue-coverage effects, not evidence that RustRed lacks the necessary row algebra.

The resulting double-cut reduction still leaves 59 explicit child residual terms; this is not a closed spanning basis or an amplitude result. Complete source-bound programs, exact requested boxes, conditions, nearest guard/proof records and 11 measured controls are saved in [the historical coverage probe](../reports/validation/2026-10-09-finite-density-native-assembly/E7-origin-native-coverage-probe/README.md). No source row, native order, physical guard or elimination implementation was changed.

### Requested-index rays with the vendored residual fix

The [current checkpoint](../reports/validation/2026-10-09-finite-density-native-assembly/requested-rays-checkpoint.json) compares the same physical E7 sources with requested-index rays disabled and enabled. Both controls use the committed native residual aggregation fix, sealed free-virtual zero sectors, ordinary zero faces, depth 3, 8192 discovery domains, three closure rounds, frontier cap 1024 and requested-label cap 4096. Adaptive guard passes are disabled. A ray starts at the actual requested nonzero integer; a zero coordinate is fixed. This changes discovery domains, while every resulting rule retains its exact native guard and proof.

The single-cut control changes from frontiers 318/1726 in 150.24 seconds to 106/538/1608 in 41.72 seconds. Both exceed the frontier budget. The double-cut control without rays times out at 180 seconds after frontier 172; with rays it reaches frontiers 68/308/899 in 123.69 seconds and exhausts the three-round budget. These are bounded coverage improvements, not closed connections. All programs, failures and resource records are preserved losslessly; timings are from concurrent processes on a shared host.

Seventy-five nonignored regression tests pass, together with complete massless fixed-D and Laurent runs (70 independent-reference comparisons and 54 refinement comparisons). A massive finest-profile regression passes all ten vacuum/cut/total reference comparisons and ten historical comparisons. Its final double-cut spanning basis has ten elements after exact residual cancellation, versus eleven historically; the evaluated target values agree. No master-minimality claim follows from that change.

Subsequent checkpoints implemented target-priority scheduling while preserving explicit `DomainBudget` gaps for unvisited domains and exact source replay. Broader four-loop numerical acceptance remains outstanding.

### Three-loop active-target and direct-zero checkpoint

The 2026-10-10 [canonical checkpoint](../reports/validation/2026-10-09-finite-density-native-assembly/active-zero-canonical-checkpoint.json) passes 167 selected tests: 78 native guarded tests and 89 RustFlow tests, with one native test explicitly ignored. It binds 1,811 source files and ten release artifacts. An earlier failure in two new queue tests was corrected by canonicalizing their source-replayed fixture RHS order; the failed logs remain available.

Native application now expands harder integrals first, according to the existing replayed order, so descending paths coalesce before shared children expand. On the saved three-loop program, this removes repeated work and preserves every previously completed coefficient map. A bounded direct-zero prepass also finds previously hidden one-original-row Lorentz zeros. Both changes retain guards, nonzero conditions and exact proof replay. Neither change establishes three-loop closure.

The opt-in frontend closure follows the actual original target sums and their reachable derivatives. It retains historical requests for audit without turning them into additional targets, and excludes individually unresolved retired rows from final candidate data. Under this policy the complete massless two-loop fixed-dimension and Laurent checks pass 70 independent-reference comparisons and 54 refinement comparisons. Those numerical binaries precede the test-fixture-only correction; their exact source snapshot is retained separately.

The complete three-loop scalar-plus-raised run still stops in its first single-cut sector after eight rounds at the explicit 1,025-label frontier limit. The isolated double-cut control reaches the same limit after eight rounds. Neither run produces numerical values. An isolated derivative query is covered by point or priority-guided native discovery, but its reduction still leaves twenty uncovered children; local coverage is not connection closure. The user-requested moving-shell/fixed-spatial-support deformation is being investigated separately before choosing the next development route.
