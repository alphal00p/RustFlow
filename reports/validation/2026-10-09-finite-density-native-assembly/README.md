# Native finite-density assembly work, 2026-10-09

This report covers work after foundation checkpoint `530b7c8`. The foundation
snapshot and its input/oracle provenance remain in
[`../2026-10-09-finite-density-foundations/README.md`](../2026-10-09-finite-density-foundations/README.md).
The selected complete massive-sunset amplitude passes at D=12/5 and through
Laurent order zero. Mandatory four-loop acceptance remains pending; no numerical
AMF-to-supplied-oracle comparison has been performed.

## Current checkpoint: complete massive amplitude and Laurent coefficients pass

The complete massive sunset now passes independent comparison at epsilon=4/5
(D=12/5), for both the scalar target and the raised original medium-numerator
target. Its vacuum, two single cuts, two-cut contribution and assembled total
are compared independently. The raised reference includes its nonzero moving
Fermi-surface terms. All **40 reference comparisons and 30 independent
refinements** pass. Predictions were saved before the comparison utility read
the independent reference.

The configurations `(digits, series order, occupied start scale)` are
`(18,60,8)`, `(28,60,8)`, `(28,80,8)` and `(28,80,12)`, with 40 guard digits and
native frontier-sector search enabled throughout. The criterion is relative
agreement within 1e-12 for magnitudes at least 1e-20 and absolute agreement
within 1e-25 below that threshold. The largest relative reference discrepancy
is 9.44882265757e-36; the largest relative refinement difference is 9.83297e-50.
These observed discrepancies do not provide rigorous interval error bounds;
the independent reference has its own empirical quadrature/precision estimates.
See [`full-sunset-frontier-reference-comparison.json`](full-sunset-frontier-reference-comparison.json)
for exact decimals, criteria and hashes of all input artifacts.

The occupied source-closed bases have sizes 6, 6 and 11. Native search of
provisional sectors, including auxiliary-constant sectors, uses the same source
corpus and guarded replay; no minimality or complete integer-domain coverage is
claimed. The four-profile prebuilt native process passed in 153.2576642698
seconds with peak RSS 52372 KiB, excluding compilation and Nix startup. Its
predictions and source checkpoints are in `full-sunset-frontier-sectors/`, with
resources in
[`full-sunset-frontier-sectors-resources.json`](full-sunset-frontier-sectors-resources.json).
The supporting source-bound regression batch passes 103 selected gates; see
[`frontier-sector-regressions.json`](frontier-sector-regressions.json).

A controlled rerun retained the old 7/6/64 bases and the `(18,60,8)` profile,
but increased guard digits from 20 to 40. It passed in 134.7241758639 seconds
with peak RSS 193684 KiB. This isolates inadequate guard precision in the
earlier endpoint failure from the later basis-size improvement; it does not
license ignoring an uncancelled endpoint diagnostic. The old and new saved
predictions agree for all ten sector/target values at their matching profile,
with largest relative difference 2.29744838567709e-46. This is a numerical
agreement test, not an exact master-basis transformation. Both predictions,
their input definitions and their closed native programs are hashed in
[`full-sunset-basis-comparison.json`](full-sunset-basis-comparison.json).

Complete Laurent refinement now also passes for both targets and orders
`[-2,-1,0]`: **30 independent-reference comparisons and 24 refinements** across
five saved profiles. The first four independently change digits, series order
and occupied start scale as above, using epsilon-grid denominator 1000. The
fifth repeats `(28,80,12)` with denominator 2000; guard precision and source
search policy stay fixed. The largest relative reference discrepancy is
1.429870667895067e-18 and the largest relative refinement difference is
9.249528747539574e-45. Laurent comparisons require relative agreement within 1e-12 for magnitudes
at least 1e-18, absolute agreement within 1e-20 below that threshold, and
imaginary-zero agreement within 1e-20. They retain the empirical-reference
qualification of the fixed-dimension gate. See
[`full-sunset-laurent-reference-comparison.json`](full-sunset-laurent-reference-comparison.json)
for all decimal differences and prediction/reference hashes.

The initial bounded process reached its 1800-second limit after saving three
profiles. Those predictions were preserved; the final two were resumed with
the exact original prebuilt binary and completed in 1238.984469 seconds with peak RSS 55296 KiB. The
run logs, resource records and timeout remain separate evidence rather than
being reported as a single uninterrupted success. This comparison covers the
complete sum, including the nonzero vacuum and all occupied contributions with
the original raised numerator and moving support.

All mandatory four-loop native predictions and AMF-to-oracle comparisons remain
absent. Earlier failures below are retained as historical evidence and are not
the current massive-small-graph acceptance result.

## Historical guard-refinement checkpoint: closure passed, endpoint failed

The guard-refinement checkpoint passes 91 selected regression gates, including
33 finite-density unit gates and the refreshed input, normalization, measure,
boundary and ordinary-flow regressions. The Python-feature all-tests type check
also passes; this does not claim a Python host import. Exact commands, source
snapshot and counts are in
[`guard-refinement-regressions.json`](guard-refinement-regressions.json).

Native requested-target and derivative closure now passes for **all three
occupied sectors** of the selected massive sunset. The full-sample preparation
retains source-replayed closed bases of sizes 7 for cut `[0]`, 6 for cut `[1]`,
and 64 for cut `[0,1]`. Their original guarded programs and closed checkpoint
records are in `full-sunset-guard-refinement/cut-{0,1,01}/`. These are finite
closed generating sets, not claims of minimal master counts or complete
coverage of every integer index domain. Previously unresolved parent domains
remain recorded. No source equation or physical support assumption was added
by splitting search domains.

The complete sample was then attempted at epsilon=4/5, 18 working decimal
digits plus 20 guard digits, series order 60 and occupied start scale 8.
Guard refinement used at most 3 passes, 256 added domains and interval width 2.
The N9 boundary supplied 51 region series and 53 coefficients, and transport
reached the physical endpoint. Endpoint extraction returned
`Numerical("uncancelled physical endpoint divergence")`; the test failed
(exit 101). This reported numerical failure does not establish divergence of
the original positive-mass integral. No complete prediction was written or
accepted, and no complete independent-reference or oracle comparison occurred.

The test harness reported 121.18 seconds. The enclosing prebuilt-process
measurement was 121.2337249289 seconds with peak RSS 184372 KiB, excluding
compilation and Nix startup. See
[`full-sunset-guard-refinement-resources.json`](full-sunset-guard-refinement-resources.json)
and its adjacent log. The earlier 59.2705-second failure below remains an
unchanged historical closure failure. At this checkpoint the endpoint failure
became the blocker; the later successful high-guard results above supersede it.

## Normalization adapter

`src/finite_density/normalization.rs` provides exact native mixed-measure,
quadratic-index, polynomial Wick and unexpanded MSbar conversion functions.
Virtual native loops use d^Dk/(i*pi^(D/2)); independent occupied loops use
d^Dq/pi^(D/2) theta(q0) theta(mu-q0) C_n(q²-m²). The future-shell convention is
q_M=(-i*P_E0,-P_Espatial), so g_E=-g_M and u_E=+i*u_M.

For L loops and k occupied coordinates the measure-only conversion is
(2*pi)^k*(4*pi)^(-L*D/2). Each original quadratic index contributes (-1)^n for
ordinary and required-cut slots alike; no extra per-cut sign is applied.
Routing determinants remain separate geometry data. `OccupiedCutFamily`
coefficients already include index and Wick phases. Compact-shell Euclidean seed
values already have their physical signs and spatial measures and must not be
converted a second time.

## Verification status

| Check | Status and scope |
| --- | --- |
| Joint release type check | Current frontier-sector checkpoint passes `cargo check --locked --release --features python --tests`; see its source-bound regression report |
| Exact sunset cut phases and odd medium numerator map | Passed |
| Unexpanded MSbar identity for L=1..4 and k=0..L | Passed at 60-digit working precision with a 50-digit comparison |
| Massive one-loop complete vacuum plus occupied assembly | Passed with independent 40/60-digit settings, raised power and support-threshold checks |
| Native massive sunset mixed and fully occupied leading boundary coefficients | Passed against independent Euclidean seed products, with a 35-digit comparison |
| Native occupied weighted closure | All selected sunset cut sectors close with exact replay; current full-sample basis sizes 6, 6 and 11 |
| Massive sunset vacuum and first single-cut endpoints | Passed both original targets at D=12/5 against independent quadrature: 16 comparisons across four configurations, plus 12 independent refinement checks |
| Current native input regression | Refreshed frontier-sector binary passes 9/9; direct wall time 0.1008 seconds; binary hash and log retained |
| Complete massive sunset at the physical endpoint | Passed at epsilon=4/5: all vacuum/cut/total values for both targets, 40 reference comparisons and 30 refinements |
| Complete massive sunset Laurent coefficients | PASS: five profiles, 30 independent-reference checks and 24 refinements for both targets through order zero |
| Four-loop predictions and Laurent stability | Not performed |
| AMF prediction records compared against the supplied oracle | 0 |
| Independent analytic references compared against the supplied oracle | I37: three exact coefficient identities; reference-only check |

All four normalization integration tests passed with no ignored tests. Direct
execution of the refreshed linked binary took 0.3690 seconds and maximum resident memory
12328 KiB; compilation is excluded. See `normalization-regressions.log` and
`normalization-regressions.json` for exact output and measurement scope.

The sunset boundary test goes through the native region expansion, hard tensor
projection and `IntegratedOccupiedBoundary` owner. Its hard factor uses recursive
ordinary AMF with bubble shortcuts disabled and the tadpole-only terminal policy.
The test then applies the whole original target phase and the native mixed measure.
At D=3, charged masses 1/2 and chemical potential 1, the expected leading
Euclidean coefficients are -1/(64*pi²) for the mixed region and +1/(64*pi²) for
the fully occupied region. These are coefficients at large auxiliary mass,
not numerical values of the full graph. Their references come from independent
Gaussian and compact one-shell integrations, not the supplied oracle.

The one-loop test has a nonzero massive vacuum term. It independently assembles
the dimensionally continued D=3 answer -max(m,mu)/(4*pi), and checks cancellation
of the raised-power vacuum and occupied pieces above threshold. It does not
establish contour admission or dimensional continuation for a multiloop weighted
family.

## Reproduction

```sh
nix develop --command cargo check --locked --release --tests
nix develop --command cargo test --locked --release --test finite_density_normalization
```

Replace pending rows only from actual command output. Compilation time is not an
evaluation benchmark. Independent working-precision checks here do not satisfy
the requested four-loop Laurent stability gate, which also requires boundary,
start-scale and epsilon-grid refinements.

Remaining mandatory work includes complete weighted differential closure and
transport, recursive retained soft-pole boundaries, finite-density endpoint
projection, full vacuum-plus-occupied assembly, certified multiloop common-contour
admission and all requested independent references and numerical acceptance runs.

## Independent massive reference

`tests/finite_density_reference.rs` supplies a separate, explicitly invoked
reference generator at D=12/5 (epsilon=4/5), where both selected massive targets
are UV convergent. It uses direct Schwinger vacuum integration, reference-only
Feynman parameters and compact radial/angular quadrature. It imports only native
arithmetic from RustFlow; it does not call AMF, reduction, boundary or compact
distribution owners and reads no oracle or feature predictions.

The derivation in
[`docs/finite-density-reference.md`](../../../docs/finite-density-reference.md)
has been independently reviewed by two agents. It retains the original raised
momentum numerator, differentiates moving support, and reports nonzero upper
surfaces separately. Six vacuum sectors and exact fifth-power radial/angular
maps remove the endpoint singularities. Gauss-Legendre node and working-precision
refinements are required before a reference is marked stable. This supplies
finite-dimension validation points, not Laurent coefficients at epsilon=0.

Explicit execution passed both tests, with no ignored tests remaining in that
run. The 12 reported components/totals passed all refinement gates. The largest
48-to-64-node absolute change was 1.069e-28; the largest 40-to-60-digit change
was 2.760e-45. Total direct runtime was 5.0437 seconds and maximum resident memory
12292 KiB, excluding compilation. The estimate comes from independent node and
precision changes, not a rigorous interval bound. Results and source digests are
in `independent-reference/independent-reference.json`; invocation/resource data
are in `independent-reference-regressions.json`. No feature predictions have
been read or compared by the generator. Reproduction:

```sh
RUSTFLOW_DENSITY_REFERENCE_REPORT=reports/validation/2026-10-09-finite-density-native-assembly/independent-reference \
  nix develop --command cargo test --locked --release --test finite_density_reference \
  generate_independent_massive_sunset_reference -- --ignored --nocapture
```

The prism supplemental target was also corrected to raise its charged slot 1;
slot 0 is neutral. The deterministic acceptance generator now records the exact
raised charged slot for every supplemental graph target. The supplied oracle
targets and preserved oracle file remain unchanged. The regenerated exact
manifest checker and all nine native input regression tests passed with the
updated supplemental definition. This rerun also checks all 108 exact original
numerator/indexed-term identities; it does not compare numerical oracle answers.

## Native closure and full-amplitude attempt

Native `OccupiedCutFamily` constructs the oriented cut coordinates and exact
target phases. Weighted preparation derives its Lorentz and distribution
identities, invokes guarded native discovery/replay, and builds a differential
connection only after target and derivative closure is established. Strictly
positive shell masses justify the recorded lower-energy support zeros. Source
conditions retain their original denominator domains; the numerical connection
checks them after exact epsilon specialization.

The shared connection owner, integrated occupied boundary owner, symbolic
power/log matching and ordinary zero-cut owner are connected in Rust. Earlier
regression results in `connection-regressions.json` record 89 passing assertion
gates against `compiled-source-hashes.json`; subsequent support-aware sources
have a separate `support-source-hashes.json` snapshot. A passing diagnostic
harness does not mean its weighted closure result was successful: each native
case retains an explicit `closed` or `unresolved` outcome.

The first complete-sunset attempt is preserved in `full-sunset-attempt-1/`.
Its two single-cut sectors source-closed with 17 and 13 masters, respectively.
The fully occupied N9 sector exhausted the provisional-frontier budget after
seven rounds, with 267 provisional integrals. The process returned failure
(exit 101), taking 59.2705 seconds and peak resident memory 33848 KiB. Compilation
and Nix startup are excluded from these figures. The failure occurred during
preparation, before any numerical vacuum/cut values or full prediction were
written; no independent-reference comparison was made.

A separate focused native vacuum/single-cut test passed (1/1). It evaluates both
original targets at epsilon=4/5 using configurations `(digits, series order,
occupied start scale)` of `(18,60,8)`, `(28,60,8)`, `(28,80,8)` and `(28,80,12)`.
Its measured process runtime was 11.3105 seconds and peak resident memory
30768 KiB, excluding compilation and Nix startup. Original predictions and the
native execution log are preserved in `single-sunset-attempt-1/`.

Only after those predictions were saved, the independent quadrature snapshots
were read by `tools/finite_density/compare_massive_reference.py`. The raised
single-cut reference includes both its bulk and its nonzero moving-support
surface. All 16 reference comparisons and all 12 independent changes of working
precision, series order or start scale passed. For values with magnitude at least
1e-20 the criterion is relative error at most 1e-12; below that magnitude it is
absolute error at most 1e-25. The predicted imaginary parts also pass the 1e-25
absolute zero check. The largest relative reference discrepancy is
2.3987682841272801e-30 and the largest relative refinement difference is
2.3987597043329849e-30.

| Contribution, final `(28,80,12)` configuration | Absolute difference from reference | Relative difference |
| --- | --- | --- |
| Scalar vacuum | 3.064155723e-43 | 1.606204298e-41 |
| Raised original numerator, vacuum | 1.582366217e-43 | 1.606168351e-41 |
| Scalar cut 0 | 1.332405863e-38 | 1.270643363e-36 |
| Raised original numerator, cut 0 | 2.130787669e-38 | 8.579690196e-36 |

`single-sunset-reference-comparison.json` retains the full saved decimal strings,
signed real and imaginary differences, complex norms, criteria and SHA256 digests
of every compared artifact. Agreement below the observed quadrature refinement
changes is an observed discrepancy, not a proven accuracy bound. The original
generator's zero-comparison metadata remains unchanged because it records its
independent generation stage. Reproduce the separate comparison with:

```sh
python3 tools/finite_density/compare_massive_reference.py
```

These four historical unique sector/target comparisons did not complete the
amplitude. The current checkpoint above now compares both single cuts, the
two-cut term and the complete sum. At that earlier snapshot, complete Laurent refinement and all
mandatory four-loop runs remain open. No AMF
prediction has been compared against the supplied oracle, and no completed
generic evaluator is claimed by this stage.

The new routing-dependent heavy-edge contour argument is documented in
[`finite-density-contours.md`](../../../docs/finite-density-contours.md).
It is a sufficient uncut-amplitude continuation certificate, separate from shell
threshold, convergence, reduction and numerical acceptance requirements.

Later independent four-loop reference generation and a separate supplied-I37
reference cross-check are recorded in
[`independent-e7-reference/README.md`](independent-e7-reference/README.md).
That cross-check compares one supplied record's three coefficients against an
already saved independent analytic reference; all three agree exactly. It is
reference versus reference, with zero AMF predictions read and zero AMF-to-oracle
comparisons. Earlier zero-comparison statements describe their historical stage.

The mathematical massless-channel checker is recorded in
[`massless-channels.json`](massless-channels.json). All 3 E7 and 7 prism occupied
cut sets pass the sufficient positive-eta criterion; 12 of 29 cut sets of the
eight-edge family pass, leaving 17 unresolved by this certificate. The
[common-contour derivation](../../../docs/finite-density-massless-contours.md)
and [dimensional compact-distribution construction](../../../docs/finite-density-massless-distributions.md)
are mathematical preparation for flowing graphs. They do not admit massless
flowing amplitudes or justify their eta=0 endpoint projector. The separately
admitted polynomial-only massless compact terminal has its own dimensional
origin prescription and tests in `frontier-sector-boundary-gates.json`.
