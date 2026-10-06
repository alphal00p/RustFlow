# Python API and Higgs-plus-jet notebook status

Snapshot: 2026-10-06, after successful final publication-wheel acceptance. The
Python API and native Higgs-plus-jet numerical milestone have passed; the broader
AMFlow/DiffExp implementation is not feature-complete. This page distinguishes
the demonstrated calculation from remaining generality and packaging work.

| Plan item | Implemented and checked | Remaining acceptance |
|---|---|---|
| Descriptive Python API | Optional PyO3 bindings in this crate, registered under `symbolica.community.hep.integration`; shared native HEPKit/Symbolica objects, scoped reductions, precision evidence, progress/cancellation and generated stubs | Advanced Rust Frobenius/infinity/endpoint operations, custom physical routes and cut-graph evaluation are not separate Python interfaces |
| Unified numerical evaluation | Automatic AMF, Euclidean FT, finite-epsilon sampling and Laurent reconstruction share the series machinery used for physical transport | Full AMFlow/DiffExp generality and performance parity remain open; see the [comparison table](benchmark-comparison.md) |
| Original AMFlow acceptance | All four two-loop paper targets through epsilon zero passed the requested 20-digit comparison and independent sampling/refinement | That completed gate does not establish full upstream feature coverage |
| Growing boundary cache | Binary persistence, compatible-source selection, intermediate points, source fingerprints and retained uncertainty; completed samples are stored separately. Final-wheel forced recomputation, interruption/resume and nearby reuse passed | Serialization overhead remains a performance issue |
| Exact Higgs-plus-jet inputs | Native physical families, certified supplied 48/61 canonical basis maps and parent connections, exact normalization and kinematics-dependent form-factor projections | Top-level RustRed derivation of those parent connections is not demonstrated; general multiloop reduction coverage remains open |
| Fresh physical boundaries and amplitude | All 16 native 30-digit configurations from empty numerical caches, all 4,360 transport coefficients, eight W/Z form factors and three coherent observables passed. Independent 40-digit regeneration, interruption/resume, uncertainty consistency and nearby-point checks also passed | External EW-square comparison remains limited to 19 relative digits; no matched full-amplitude reference timing |
| Independent Euclidean anchors | Publication wheel passed all 545 comparisons at 20 digits with fresh native boundaries | This prerequisite alone does not certify the physical calculation |
| Marimo notebook | `symbolica-community/examples/hep/gg_hg.py`, stage controls, precision/provenance tables and cache controls; complete headless acceptance passed. HTML export, native SVG display and populated original notebook execution passed using a copied completed native cache. Actual Chromium transport, amplitude and restart clicks passed | Boundary-generation and cancellation controls have native/controller evidence, not browser-click coverage |
| Packaging | Rust core release suite: 534 passed, 12 opt-in tests ignored; strict Python-feature Clippy and formatting passed. Original native wheel component suite: 206 passed, one inherited C++ export failure. Corrected tensor stubs: seven checks passed. Isolated tensor fix passes native and actual Pyodide checks with/without NumPy | Existing full browser wheel still contains the inherited tensor list-output bug; its fix has not been integrated into the frozen host. The whole browser suite has not been rerun with that fix |

The final publication run completed
[all 16 physical starting configurations](../reports/validation/2026-10-05-public-wheel/physical-boundary-validation.json)
at 30 verified digits, with 1,248 completed finite-epsilon samples. Fresh boundary
generation took 6,467.54 seconds. The
[cold physical calculation](../reports/validation/2026-10-05-public-wheel/cold-physics-validation.json)
also passed transport and native amplitude comparison, with the same 25/26/47
propagated relative-digit estimates as the earlier build. The complete
[publication acceptance](../reports/validation/2026-10-06-gg-hg-publication-complete/summary.json)
subsequently passed independent 40-digit regeneration, all uncertainty checks,
interruption/resume and nearby reuse. The root process exited with code zero;
final runtime attestation and post-completion source hashes passed before any
presentation or compatibility follow-up was applied.

The [independent scope audit](../reports/validation/2026-10-06-notebook-finalization/independent-scope-audit.json)
distinguishes the exact supplied parent equations from computed numerical data.
Generic common-mass transformations use those certified equations to close the
parent flow. Native region enumeration, recursive boundary-family IBPs, analytic
terminals and finite-epsilon reconstruction then generate the boundary values.
This route admits supplied symbolic reduction results; it does not establish
that RustRed independently rederived all 48/61 parent master equations.

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
| Resumed independent 40-digit boundaries | 10,317.13 s | 5,503.79 s |
| Refined transport | 426.26 s | 412.67 s |
| Refined amplitude | 51.64 s | 36.88 s |
| Nearby transport, all 16 reuse accumulated points | 374.08 s | 363.54 s |
| Nearby amplitude | 4.23 s | 4.02 s |

The historical build subsequently
[completed its entire acceptance with exit code zero](../reports/validation/2026-10-05-gg-hg-historical-complete/summary.json).
All sixteen independent 40-digit boundary fits passed; all 4,360 transported
coefficients were consistent with the cold run's uncertainty estimates, and
the refined form factors and observables passed their checks. Resumed boundary
generation took 10,317.13 seconds after a 299.42-second intentionally interrupted
stage; refined transport and amplitude took 426.26 and 51.64 seconds. Nearby
transport reused accumulated physical points for all sixteen configurations,
inserted 32 points, and took 374.08 seconds, followed by 4.23 seconds for amplitude
assembly. Completed sample hashes survived the interruption and resume; final
runtime attestation passed. The publication wheel now independently passes the
same complete gates: 1,248 initial and 1,376 refined samples, all 4,360 stable
coefficients, all eight form factors and three observables, and 32 points
inserted by nearby transport. Each report retains its own runtime provenance.

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

A separate [cache-save optimization PR](https://github.com/alphal00p/RustFlow/pull/1)
reuses native Symbolica exports within each snapshot. Four new release tests,
forty existing cache integration tests, formatting and strict release Clippy pass.
A paired profile with the real canonical topology and synthetic numerical values
measured 6.23–20.47× faster payload serialization with byte-identical output.
These are isolated serialization results, not a notebook speedup; the publication
runtime remains unchanged. The review also produced a standalone
[Symbolica export-race reproducer](../mre/symbolica-export-registry-race/README.md).
Its validated owner fix is proposed in
[Symbolica PR #52](https://github.com/symbolica-dev/symbolica/pull/52) and is not
part of the tested runtime. Native import/export regressions, standalone
reproducers, generator reentrancy checks and strict Clippy passed in isolation.

The separate [native-arithmetic follow-up](https://github.com/alphal00p/RustFlow/pull/2)
uses official Symbolica's explicit rounded scalar operations in the complex
arithmetic layer. It preserves the working precision and operation order;
24,192 bit comparisons, 76 selected regressions, formatting and strict Clippy
passed. Its report makes no new scientific acceptance or performance claim,
and it has not changed the publication runtime.

The [interactive notebook check](../reports/validation/2026-10-05-public-wheel/notebook-interactive-validation.json)
served the unchanged notebook with native Marimo and clicked transport, amplitude,
binary reload and repeated amplitude controls in Chromium. All completed without
browser errors. A separate panel capture confirmed the populated transport and
eight-row form-factor tables and native diagram display. Both checks used copies
of the committed cold bank; its originals, runtime and steering sources were
verified unchanged. Click-to-completion timings include UI refresh and persistence,
and are not isolated solver measurements.

A documentation-only notebook follow-up was validated at community commit
`059ec1faf663e38ccb4def8bfeb368eb7b5edc8e`, on
`codex/gg-hg-notebook-finalization`. Its exact-input provenance, comparison
precision and explicit loop measure passed
[static Chromium rendering and notebook smoke checks](../reports/validation/2026-10-06-notebook-finalization/render-validation.json).
The [measure display](../reports/validation/2026-10-06-notebook-finalization/normalization.png)
was visually inspected. This empty-cache presentation check starts no numerical
work. Its wording and the separately tested Python compatibility changes were
applied after the full numerical run's final attestation and source rechecks.

## Publication and remaining work

RustFlow milestones and validation reports are pushed to its repository. The
community notebook/API changes and dependency pins are pushed in
[draft PR #19](https://github.com/symbolica-dev/symbolica-community/pull/19),
authored by `ValentinHirschi`. Its head is
`5c204ba34b836d448f9d166d546baecca8fda8e3`; the full numerical run used
`8758b93d1e1a6fed2aace41493a91be569f98a98`. Native dependencies, extension bytes,
the calculation controller and anchor runner are unchanged between these
revisions. The differences are presentation, polling compatibility and the
documented CI repairs. Merge and release remain pending upstream/CI work. The
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
regeneration resumed in a fresh session and preserved completed sample hashes.
All sixteen refined sources passed, using 1,376 new-generation samples. Refined
physical transport and observable stability passed, and all sixteen nearby
configurations reused accumulated physical points with preserved provenance.
The [complete report](../reports/validation/2026-10-06-gg-hg-publication-complete/acceptance.json)
records the final integrity attestation; its summary includes the authoritative
process exit rather than inferring success from a progress file.

Separate Python 3.10 checks found and fixed the headless runner's use of the
built-in timeout exception: `Future.result` raises a distinct exception on that
Python version. The follow-up preserves pending polling, worker exception
identity and final failure progress. It passed 72 focused checks on Python 3.10
and 70 on Python 3.14 with identical native extension bytes; the difference is
two tests for exception classes that become aliases on newer Python. A separate
conditional `tomli` dependency fixes the inherited Python 3.10 test-collection
failure. Both commits are pushed on `codex/python310-toml-tests`, ending at
`d36b9f4fe2b2b0766f313dab05237e90ea544640`, and are now included in the prepared
publication branch. The final applied checkout also passed all 70 focused
Python 3.14 controller, anchor and notebook checks. See the
[compatibility evidence](../reports/validation/2026-10-06-python310-compatibility/polling-validation.json)
and [applied-steering checks](../reports/validation/2026-10-06-gg-hg-publication-complete/prepared-steering-validation.json).

The full Python 3.10 suite before the polling fix recorded 1,095 passes,
36 failures and 172 skips. Its
[source-provenance audit](../reports/validation/2026-10-06-python310-compatibility/full-suite-provenance.json)
identifies 29 missing native Linnet dependency failures, inherited documentation
and accessor issues, an unrelated showcase symlink issue, the polling defect and
the known C++ export defect. These are being repaired or tracked separately;
the focused green checks do not mean the entire community suite is green.
The [subsequent documentation/accessor repair](../reports/validation/2026-10-06-python310-compatibility/ci-repair-validation.json)
passed 653 documentation examples and 32 focused checks, retaining only the same
29 missing-Linnet failures in that documentation run. API signatures were
unchanged. Verification with a matching native Linnet wheel remains separate.

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
