# Finite-density implementation and acceptance plan

This is an implementation plan and a record of exact input certificates. It is not
an end-to-end numerical acceptance report. The checkout audited before changes is
`b3a4843e8327835d1ec3ada5a6f32f1841bab2c2`. The supplied schema-3 oracle has been
preserved byte for byte in `fixtures/finite_density/oracle_results_supplied.json`.
Its SHA-256 is `d49acc63e40b9fdccf7fb3eddea6e71f211aba0067021e917925274952ad66d1`.

The concrete machine-readable contract is
[`finite-density-acceptance.json`](finite-density-acceptance.json). Predictions,
reference comparisons and precision successes are initially empty. Every missing
independent reference is explicitly marked; none may become a boundary value.

## Owners and evidence gates

| Gate | Existing source owners | Required evidence |
| --- | --- | --- |
| Native dependency integration | `Cargo.toml`, `Cargo.lock`, `src/reduction.rs` | One compatible Symbolica/RustRed/HEPKit dependency graph, guarded replay test, ordinary/cut/transport regression results |
| Physical input and measure | `src/hepkit.rs`, `src/family.rs`, `src/cuts.rs` | Exact supplied incidence/routing/charge certificate, Symbolica polynomial conversion including medium coordinates, distributional identities and explicit contour admission |
| One complete small graph | `src/cut_flow.rs`, `src/cut_regions.rs`, `src/regions.rs`, `src/integrand.rs`, `src/tensor.rs` | Zero-, one- and two-cut massive sunset, raised charged line, original polynomial held fixed under independent mass differentiation, surface contribution |
| Boundary extensions | `src/recursive.rs`, `src/boundary.rs`, `src/frobenius.rs` | Tadpole-only recursive policy propagated through children/cache identities; projected hard/soft expressions retained; weighted soft integration; logarithmic region matching; bounded recursion |
| Generic flow | `src/reduction.rs`, `src/physical_family.rs`, `src/epsilon.rs`, `src/frobenius.rs` | Guarded target/derivative closure and conditions, rational target projection before the endpoint limit, independent precision refinements |
| Four-loop acceptance | `examples/finite_density`, supplied oracle and independently generated references | Three certified distinct nonfactorized families, raised and numerator targets, a massless endpoint, at least ten stable digits, saved predictions before comparison |

Boundary changes necessary for the sunset belong before its evaluation gate. A
failed small example must identify the missing owner or mathematical condition;
it is not a reason to replace recursion with an analytic sunset formula.
`RustRedBackend::bubble_subloops` must be false on this acceptance path. Ordinary
callers keep their established shortcuts.

The baseline lock selects RustRed `acc92b0dad27b11fd194a4c284765fb6a93cbc94`,
Symbolica community `ed2374f1d880d52c3a7ca48cd7c22f4baad5c020` and GammaLoop/HEPKit
`635d3a1feb44067583a668c7d18f67405fe144e1`. The guarded RustRed candidate is
`78969aab524b7d6a2eec36f59b01e9af1e04cc04`. Integration must inspect source identity,
companion packages and the embedding host; this list does not certify that a
combined build has passed. Generated build/test reports record the final pins.

## Graph evidence and exact input isolation

HEPKit supplies graph, routing and tensor facilities through GammaLoop and Linnet;
its [project inventory](https://hepkit.org/#projects) identifies these owners.
For model-free scalar inputs, `linnet::half_edge::builder::HedgeGraphBuilder` can
build the supplied incidence directly. `GraphIntegral::new` in `src/hepkit.rs`
validates a full native model diagram when one is available. A routing matrix
alone is not a physical graph certificate.

The three four-loop certificates use seven, eight and nine physical edges,
respectively, and four, five and six vertices. The last is a triangular prism.
The certificates record directed incidence B, routing R, exact B R = 0, full
cycle-space dimension, chemical conservation and connected complements for every
retained occupied cut. Every maximal nonzero routing minor has determinant +1 or
-1. Physical edges exclude numerator completion slots. Full vertex-permutation
canonicalization preserves multiplicities; distinct vertex/edge counts and
canonical adjacency establish that these are different multigraphs. Exhaustive
exact rank tests rule out splitting each routing into independent denominator
blocks. Thus the nonfactorization evidence goes beyond connectivity.

Run the deterministic validation-only checker with:

```sh
python3 tools/finite_density/acceptance_manifest.py
```

`--write` regenerates its JSON artifacts. It uses Python rational arithmetic only
for finite graph certificates, not polynomial conversion or integral evaluation.
The oracle's algebraic expressions must be parsed with native Symbolica in Rust
validation. `examples/finite_density/oracle_definitions.json` includes exactly
108 input definitions and basis maps, with no constants, answer coefficients,
uncertainties or determinant tags. Its original per-target indexed coefficients
are input polynomial identities, not reference answers.

Mandatory supplied four-loop records are I37, I91 and I115. I37 includes a raised
neutral denominator; it does not provide a raised occupied-line reference. Each
family additionally requests a raised occupied line with a polynomial numerator.
The seven-edge family's supplemental reference has now been independently
derived by regulated tensor/beta integration and checked using native precision
and epsilon refinements; see the [E7 reference derivation](finite-density-e7-reference.md).
The other two supplemental references still require generation. I115 supplies a finite coefficient through epsilon order zero; higher orders
remain unspecified. Its rationally encoded central value remains approximate
when an uncertainty is supplied. The absolute criterion in the manifest applies
only to zero reference coefficients. The massive two-loop graph and multiple-chemical-potential
coverage require references beyond the massless oracle. A complete independent
massive reference has now been generated at D=12/5, with 12 components/totals
passing node and precision refinements. Complete saved native values for both
targets now pass all 40 vacuum/cut/total comparisons and 30 precision/order/start
refinements at epsilon=4/5, with 40 guard digits held fixed. Native guard and
frontier-sector discovery closes occupied bases of sizes 6, 6 and 11. The earlier
7/6/64 bases also pass at the matching high-guard profile; the historical
20-guard-digit endpoint failure remains recorded.
The [current validation report](../reports/validation/2026-10-09-finite-density-native-assembly/README.md)
preserves measured differences and failure evidence. Its independent Laurent
extension through orders [-2,0] now passes 108 component-coefficient refinement
comparisons and six analytic UV-residue checks, using six-sector subtraction and
exact rational interpolation; the empirical errors are not rigorous interval
bounds. Complete native Laurent comparison remains pending.
Concrete reference-only kernels and the E8 shared-pole obstruction are recorded
in [the remaining-reference investigation](finite-density-missing-references.md).
No E8 or prism supplemental reference values have been generated, and their
precision feasibility remains unestablished.

## Remaining acceptance evidence

Exact certificates and bounded algebraic pilots are useful gates, but do not
constitute a numerical finite-density evaluator. Before claiming completion, the
saved report must include the complete amplitude for the massive sunset and all
mandatory four-loop targets, contour/placement admissions, guarded replay and
nonzero conditions, integrated boundary provenance, endpoint cancellation tests,
independent reference errors, and measured stability under separate changes in
precision, boundary order, start scale and epsilon grid. Report runtime and peak
memory with a clear scope, including whether preparation and caches are counted.

No numerical oracle comparisons are claimed by the manifest generator. Broader
oracle coverage is optional only after the mandatory targets and the additional
massive, raised-occupied and medium tests pass.
