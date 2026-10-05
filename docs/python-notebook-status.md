# Python API and Higgs-plus-jet notebook status

Snapshot: 2026-10-05, during final publication-wheel acceptance. The implementation
is not feature-complete. This page separates completed numerical evidence from
checks still running on the dependency graph intended for publication.

| Plan item | Implemented and checked | Remaining acceptance |
|---|---|---|
| Descriptive Python API | Optional PyO3 bindings in this crate, registered under `symbolica.community.hep.integration`; shared native HEPKit/Symbolica objects, scoped reductions, precision evidence, progress/cancellation and generated stubs | Advanced Rust Frobenius/infinity/endpoint operations, custom physical routes and cut-graph evaluation are not separate Python interfaces |
| Unified numerical evaluation | Automatic AMF, Euclidean FT, finite-epsilon sampling and Laurent reconstruction share the series machinery used for physical transport | Full AMFlow/DiffExp generality and performance parity remain open; see the [comparison table](benchmark-comparison.md) |
| Original AMFlow acceptance | All four two-loop paper targets through epsilon zero passed the requested 20-digit comparison and independent sampling/refinement | That completed gate does not establish full upstream feature coverage |
| Growing boundary cache | Binary persistence, compatible-source selection, intermediate points, source fingerprints and retained uncertainty; completed samples are stored separately | Final-wheel forced recomputation, interruption/resume and nearby-query acceptance are still running |
| Exact Higgs-plus-jet inputs | Native physical families, certified 48/61 canonical basis maps, exact normalization and kinematics-dependent form-factor projections | No claim of automatic reduction or boundary coverage for every possible multiloop family |
| Fresh physical boundaries and amplitude | The publication wheel completed all 16 native 30-digit configurations from empty numerical caches, all 4,360 transport coefficient comparisons, eight W/Z form factors and three coherent observables | Forced 40-digit regeneration, interruption/resume, independent stability and nearby-point checks remain |
| Independent Euclidean anchors | Publication wheel passed all 545 comparisons at 20 digits with fresh native boundaries | This prerequisite alone does not certify the physical calculation |
| Marimo notebook | `symbolica-community/examples/hep/gg_hg.py`, stage controls, precision/provenance tables and cache controls; initial HTML export, native SVG display and populated original notebook execution passed using a copied completed native cache. Actual Chromium transport, amplitude and restart clicks passed, with native diagram/form-factor panels inspected | Finish the complete 40-digit headless acceptance; boundary-generation and cancellation controls have native/controller evidence, not browser-click coverage |
| Packaging | Rust core release suite: 534 passed, 12 opt-in tests ignored; strict Python-feature Clippy and formatting passed. Original native wheel component suite: 206 passed, one inherited C++ export failure. Corrected tensor stubs: seven checks passed. Isolated tensor fix passes native and actual Pyodide checks with/without NumPy | Existing full browser wheel still contains the inherited tensor list-output bug; its fix has not been integrated into the frozen host. The whole browser suite has not been rerun with that fix |

The final publication run has now completed
[all 16 physical starting configurations](../reports/validation/2026-10-05-public-wheel/physical-boundary-validation.json)
at 30 verified digits, with 1,248 completed finite-epsilon samples. Fresh boundary
generation took 6,467.54 seconds. The
[cold physical calculation](../reports/validation/2026-10-05-public-wheel/cold-physics-validation.json)
also passed transport and native amplitude comparison, with the same 25/26/47
propagated relative-digit estimates as the earlier build. Its runtime, dependency
pins, model and numerical steering sources remain frozen while the later
restart/refinement stages continue.

## Completed cold calculation and precision scope

The [historical cold report](../reports/validation/2026-10-05-gg-hg-historical-cold/summary.json)
records native generation without reference-seed fallback, with two independent
sample sets per source and 1,248 completed finite-epsilon samples. All 4,360
transport coefficients passed the 20-digit mixed absolute/relative comparison.
The native propagated relative-digit estimates were 25 for the coherent EW
square, 26 for interference and 47 for the infinite-top HEFT square. Reference
uncertainties support only 19, 20 and 39 relative digits respectively; the EW
reference comparison therefore does not independently establish 20 digits.

| Stage | Historical build | Publication build |
|---|---:|---:|
| All 16 fresh native boundaries | 7,399.44 s | 6,467.54 s |
| All 16 physical transports | 405.37 s | 400.07 s |
| Native amplitude assembly | 50.75 s | 35.96 s |
| Binary reload and repeated transport | 75.39 s | 73.88 s |
| Warm transport, 16 exact hits and zero ODE steps | 69.44 s | 68.43 s |
| Warm amplitude assembly | 4.06 s | 4.03 s |

These observations used four boundary workers sharing a 256-sample-worker budget
on a busy cluster. Warm controller timings include persistence after each query;
they are not isolated cache-lookup timings. The
[cache audit](../reports/validation/2026-10-05-public-wheel/warm-cache-audit.json)
identified repeated serialization. A subsequent
[copied-bank measurement](../reports/validation/2026-10-05-public-wheel/populated-notebook-validation.json)
took 0.268 seconds for sixteen exact evaluations without saving, 56.221 seconds
for sixteen unchanged saves and 55.918 seconds for the original combined
controller path. Values, errors, provenance and cached evidence agreed, and the
original archived files stayed unchanged. This one sequence ran before amplitude
preparation in a fresh process on a busy host; it measures a persistence
bottleneck but establishes no optimization speedup. A matched full-application
reference timing is not yet available.

The [interactive notebook check](../reports/validation/2026-10-05-public-wheel/notebook-interactive-validation.json)
served the unchanged notebook with native Marimo and clicked transport, amplitude,
binary reload and repeated amplitude controls in Chromium. All completed without
browser errors. A separate panel capture confirmed the populated transport and
eight-row form-factor tables and native diagram display. Both checks used copies
of the committed cold bank; its originals, runtime and steering sources were
verified unchanged. Click-to-completion timings include UI refresh and persistence,
and are not isolated solver measurements.

## Publication and remaining work

RustFlow milestones and validation reports are pushed to its repository. The
community notebook/API changes and dependency pins are pushed in
[draft PR #19](https://github.com/symbolica-dev/symbolica-community/pull/19),
authored by `ValentinHirschi`. Its head is `8758b93d1e1a6fed2aace41493a91be569f98a98`;
merge and release remain pending the acceptance and upstream issues below. The
notebook's repository-relative destination is `examples/hep/gg_hg.py`; the local
publication checkout is `/common/dev/symbolica-community/loop-integration-publication`.

The inherited tensor list-output failure has a separate
[draft owner PR #126](https://github.com/alphal00p/gammaloop/pull/126).
[Focused native tests](../reports/validation/2026-10-05-public-wheel/tensor-sequence-owner-validation.json)
pass with and without NumPy on that isolated fix. Its
[actual Pyodide validation](../reports/validation/2026-10-05-public-wheel/tensor-sequence-browser-validation.json)
also passes both output types. The live wheel has not been changed; the isolated
tensor host is not a rerun of the full community browser suite.

The final wheel's forced cancellation after one fresh 40-digit sample passed;
regeneration resumed in a fresh session. Next are all sixteen refined sources,
independent stability and nearby-query checks. Update the draft's acceptance
evidence after those checks; the live calculation and its sources stay frozen.
Broader work remains on general linear/cut recursion, algebraic singular
endpoints, causal paths, nondiagonal epsilon transformations, advanced Python
solver access, and performance across the reference benchmark suite. Numerical
success for the demonstrated families is not a claim of arbitrary-family parity.

Evidence for this packaging milestone is in the
[native core report](../reports/validation/2026-10-05-public-owner-core.json),
[wheel component report](../reports/validation/2026-10-05-public-wheel/component-validation.json),
[Euclidean anchors](../reports/validation/2026-10-05-public-wheel/anchor-acceptance.json),
[initial notebook export](../reports/validation/2026-10-05-public-wheel/notebook-render-validation.json),
[native diagram display](../reports/validation/2026-10-05-public-wheel/notebook-native-display-validation.json),
[tensor stubs](../reports/validation/2026-10-05-public-wheel/tensor-stub-package-validation.json)
and [browser checks](../reports/validation/2026-10-05-public-wheel/browser-validation.json).
