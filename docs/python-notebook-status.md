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
| Fresh physical boundaries and amplitude | An earlier frozen build completed all 16 native 30-digit configurations from empty numerical caches, all 4,360 transport coefficient comparisons, eight W/Z form factors and three coherent observables | Repeat on the publication wheel, followed by forced 40-digit regeneration and independent stability checks |
| Independent Euclidean anchors | Publication wheel passed all 545 comparisons at 20 digits with fresh native boundaries | This prerequisite alone does not certify the physical calculation |
| Marimo notebook | `symbolica-community/examples/hep/gg_hg.py`, stage controls, precision/provenance tables and cache controls; initial HTML export passed, and all three native amplitude diagrams rendered as SVG and in the notebook display container | Execute the populated notebook using a copy of the completed native cache, then finish the complete headless acceptance |
| Packaging | Rust core release suite: 534 passed, 12 opt-in tests ignored; strict Python-feature Clippy and formatting passed. Original native wheel component suite: 206 passed, one inherited C++ export failure. Corrected tensor stubs: seven checks passed | Existing browser tensor evaluation fails on list output; an isolated owner fix passes native NumPy/no-NumPy regression checks, but its browser rerun is pending |

The final publication run has now completed
[all 16 physical starting configurations](../reports/validation/2026-10-05-public-wheel/physical-boundary-validation.json)
at 30 verified digits, with 1,248 completed finite-epsilon samples. Fresh boundary
generation took 6,467.54 seconds; physical transport has started. Its runtime,
dependency pins, model and numerical steering sources remain frozen. The physical
amplitude and later restart/refinement stages have not yet passed on this wheel;
earlier build results are not relabeled as its validation.

## Completed cold calculation and precision scope

The [historical cold report](../reports/validation/2026-10-05-gg-hg-historical-cold/summary.json)
records native generation without reference-seed fallback, with two independent
sample sets per source and 1,248 completed finite-epsilon samples. All 4,360
transport coefficients passed the 20-digit mixed absolute/relative comparison.
The native propagated relative-digit estimates were 25 for the coherent EW
square, 26 for interference and 47 for the infinite-top HEFT square. Reference
uncertainties support only 19, 20 and 39 relative digits respectively; the EW
reference comparison therefore does not independently establish 20 digits.

| Historical stage | Measured wall time |
|---|---:|
| All 16 fresh native boundaries | 7,399.44 s |
| All 16 physical transports | 405.37 s |
| Native amplitude assembly | 50.75 s |
| Binary reload and repeated transport | 75.39 s |
| Warm transport, 16 exact hits and zero ODE steps | 69.44 s |
| Warm amplitude assembly | 4.06 s |

These observations used four boundary workers sharing a 256-sample-worker budget
on a busy cluster. Warm controller timings include persistence after each query;
they are not isolated cache-lookup timings. The
[cache audit](../reports/validation/2026-10-05-public-wheel/warm-cache-audit.json)
identifies repeated serialization but contains no measured cost attribution or
speedup. A matched full-application reference timing is not yet available.

## Publication and remaining work

RustFlow milestones and validation reports are pushed to its repository. The
community notebook/API changes are committed on the local publication branch;
community publication remains pending the complete notebook acceptance. The
notebook's repository-relative destination is `examples/hep/gg_hg.py`; the local
publication checkout is `/common/dev/symbolica-community/loop-integration-publication`.

The inherited tensor list-output failure has a separate
[draft owner PR #126](https://github.com/alphal00p/gammaloop/pull/126).
[Focused native tests](../reports/validation/2026-10-05-public-wheel/tensor-sequence-owner-validation.json)
pass with and without NumPy on that isolated fix. The live wheel has not been
changed, and a successful browser rerun is still required.

Next are the final wheel's cold amplitude, forced cancellation/restart and
40-digit regeneration, populated notebook execution, and copied-cache profiling.
After those checks, publish the community update with its exact validated pins.
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
