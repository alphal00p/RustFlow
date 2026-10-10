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
| Three-loop intermediate acceptance | `examples/finite_density/massless_three_loop_chain.json`, native full-amplitude owners | Complete vacuum and occupied-cut sum, scalar and raised occupied-line mixed-medium targets, independently refined fixed-dimension and Laurent comparisons before further four-loop numerical scale-up |
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
bounds. Complete native Laurent comparison now passes 30 independent-reference
checks and 24 refinements across five profiles, including a separate epsilon
grid change. Required four-loop numerical comparisons remain absent.
Concrete reference-only kernels and the E8 shared-pole obstruction are recorded
in [the remaining-reference investigation](finite-density-missing-references.md).
The E8 supplemental reference remains unavailable. The prism supplemental
reference now has a complete independently refined Laurent[-4,0] result: all 20
coefficient refinements pass, with largest empirical normalized absolute
uncertainty 5.03e-16. This is reference availability, not native AMF acceptance.

At the user's request, complete three-loop numerical validation now precedes
further four-loop numerical scale-up. The frozen intermediate input is a
nonfactorized three-vertex, five-edge massless graph with cuts `[]`, `[0]`, `[3]`
and `[0,3]`. It requests a scalar and a raised occupied line with the original
numerator `g1_2+u1*u2`, including its nonzero upper Fermi-surface contribution.
The [independent reference](finite-density-three-loop-reference.md) is available
at D=13/2 and through Laurent orders [-3,0], with arithmetic, quadrature and
exact-algebra checks. Both a source-only pilot and the subsequent production
schedule close the `[0]` cut family, while complete preparation and numerical
validation remain pending; no three-loop amplitude comparison is claimed. The earlier saved
[rule-reuse controls](../reports/validation/2026-10-09-finite-density-native-assembly/three-loop-rule-union-controls.json)
and [residual-search diagnostic](../reports/validation/2026-10-09-finite-density-native-assembly/three-loop-residual-batch-probe/summary.json)
record this limitation. This intermediate gate does not replace any original
four-loop family, raised target or mandatory oracle record.

The subsequent native source-projection checkpoint retains explicit measure-zero
domains during original-row replay and versions guarded persistence accordingly.
Its initially exposed massless two-loop index-guard failure is resolved by
[bounded native point refinement](finite-density-conditional-points.md). The
massless fixed-dimension/Laurent and massive regression comparisons pass again;
the three-loop frontier still exceeds its configured limit. The separate
[fixed-spatial-occupation moving-shell assessment](finite-density-moving-shell-exploration.md)
does not establish a multi-loop advantage, so the production deformation remains
fixed-shell while reduction closure is investigated.

The subsequent [native ordering controls](../reports/validation/2026-10-10-occupation-order-experiment/README.md)
do not establish closure with occupation degree moved earlier in the order.
The [reciprocal occupied-energy exploration](../reports/validation/2026-10-10-massless-reciprocal-energy-exploration/README.md)
also gives no measured reduction benefit from widening the admitted energy
powers alone. Exact shifted sources shorten two small reductions within the
original polynomial domains and reduce early frontier growth, but the bounded
three-loop pilot still fails to close. Neither experiment changes production
ordering, massless admission or the three-loop numerical acceptance gate.

Further [correlated Lorentz-source controls](../reports/validation/2026-10-10-global-boost-source-experiment/README.md)
derive an exact polynomial boost identity and a narrowly guarded simple-cut
surface recurrence. Each improves selected reductions, but their separate
bounded full pilots remain unclosed. The narrow recurrence and the earlier
polynomial temporal source shifts retain their complementary improvements in
[combined small controls](../reports/validation/2026-10-10-combined-polynomial-source-experiment/README.md).
They remain exploratory source additions. An isolated
[direct-match search lookahead](../reports/validation/2026-10-10-native-direct-hit-lookahead/README.md)
gives no improvement; its motivating example already has a first indirect
winner. Production search remains unchanged while bounded source selection is
investigated with exact replay and explicit condition-coverage checks.

The [normalized-boost attribution](../reports/validation/2026-10-10-normalized-boost-reciprocal-experiment/README.md)
finds shorter small reductions after adding inverse occupied-energy powers, but
the [polynomial raw Ward controls](../reports/validation/2026-10-10-polynomial-raw-ward-attribution/README.md)
reproduce all 19 coefficient and condition maps without that domain extension.
Keeping the energy factor explicit in the narrowly guarded surface identity,
together with exact polynomial source shifts, gives 143 right-hand-side terms
instead of 192. Both 62-source orders retain this small-control result.
Production massless admission remains polynomial.

The subsequent [complete native closure pilots](../reports/validation/2026-10-10-native-source-portfolio-experiment/README.md)
close the `[0]` singleton family on the same seven basis integrals, with exact
agreement of the final original-target and basis-derivative coefficient rows.
Unchanged RustRed needs 16 rounds and 3,553 rules; an isolated source-selection
prototype needs 15 rounds and 3,324 rules. Their retained condition lists differ
and remain explicit. Both runs include their timed-out initial segments and
verified native continuations; shared-host elapsed times are not a comparative
performance result. These are algebraic closure results for one cut, not the
assembled three-loop amplitude.

Since unchanged RustRed closes the polynomial corpus, production integration
uses the optional `polynomial-closure-v1` source policy and keeps native search
unchanged. The [corrected integration gate](../reports/validation/2026-10-10-polynomial-closure-integration-v2/polynomial-closure-v2-root-gates.json)
passes 105 tests, including exact agreement of the public factory with the
complete frozen corpus. The unchanged native binary and sources authenticate
reuse of the earlier 84 passing RustRed tests. The initial attempt's two invalid
nonunit-routing fixtures and their correction remain recorded separately; no
production admission rule was relaxed. Full native three-loop preparation,
boundaries, transport and Laurent refinement remain pending. Those numerical
gates still precede four-loop numerical scale-up.

The integrated polynomial policy also passes the massless fixed-dimension and
Laurent two-loop regressions and the massive two-loop profile. Four bounded
production three-loop preparation controls fail before generating predictions:
raising conditional-point passes and total allocation exposes further coupled
index conditions; reducing the per-residual allocation to one leaves a growing
unreduced frontier. Their [frozen report](../reports/validation/2026-10-10-polynomial-closure-integration-v2/README.md)
preserves every failure. The successful standalone pilot additionally searches
exact points after ordinary ray misses. An optional production schedule now
transfers this existing native mechanism. Its [integration checkpoint](../reports/validation/2026-10-10-requested-ray-point-integration/README.md)
passes 115 tests and retains authenticated reuse of the 84 native RustRed tests.
The massless two-loop comparisons pass at fixed dimension (40 references and
30 refinements) and in the Laurent expansion (30 references and 24 refinements).
The massive profile passes ten independent references, ten historical comparisons
and two assembly checks. No new elimination backend or source-admission
extension is introduced by that scheduling change.

Production preparation now closes the three-loop `[0]` cut on seven basis
integrals and 3,553 native rules, matching all 16 rounds of the standalone pilot.
The `[3]` and `[0,3]` cuts still exceed the frontier cap with the new schedule.
Separate controls using the existing residual schedule also fail: the former
encounters a vanishing index condition and the latter exceeds the frontier cap.
The complete three-loop amplitude has no prediction yet. The `[3]` cut energy
is a difference of completion factors, so the current exact-`c*E` source policy
does not add its shifted temporal and boost presentations. General polynomial
reindexing for that case is being assessed separately before changing the policy.

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
