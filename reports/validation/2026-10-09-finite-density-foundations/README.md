# Finite-density foundation validation, 2026-10-09

Status: implementation foundations, **not completion of the requested generic
finite-density evaluator**. This report records exact input admission, guarded
infrastructure, compact seeds and boundary-owner extensions separately from the
missing assembled multiloop calculations. No numerical oracle predictions or
comparisons have been performed.

## Input provenance and isolation

- RustFlow baseline: `b3a4843e8327835d1ec3ada5a6f32f1841bab2c2`.
- Dependency integration checkpoint: `e189dbe`; final feature source
  identities and test logs are recorded in the companion dependency reports.
- Guarded RustRed candidate integrated: `78969aab524b7d6a2eec36f59b01e9af1e04cc04`.
- Symbolica community: `ed2374f1d880d52c3a7ca48cd7c22f4baad5c020`.
- HEPKit/GammaLoop: `635d3a1feb44067583a668c7d18f67405fe144e1`.
- Oracle: `fixtures/finite_density/oracle_results_supplied.json`, schema 3,
  preserved without modifications.
- Oracle SHA-256:
  `d49acc63e40b9fdccf7fb3eddea6e71f211aba0067021e917925274952ad66d1`.
- Complete supplied inventory: 95 four-loop records, 13 lower-loop records and
  nine constants, including all reference coefficients, uncertainties and tags.
- Definition-only export: `examples/finite_density/oracle_definitions.json`.
  Reference constants, answer coefficients, uncertainty fields and determinant
  tags are excluded; algebraic indexed-target coefficients remain input data.
- Numerical oracle records actually compared: **0**. No prediction artifact exists.

## Verification evidence

| Check | Result |
| --- | --- |
| Independent exact manifest checker | Passed: three distinct nonfactorized four-loop graph certificates, one massive two-loop graph, two independently charged cycles, 108 definition-only records and original oracle hash |
| Native Symbolica original-definition and target-conversion identities | Passed: all 108 records, routing squares, legacy dotted numerators, inverse maps and indexed target coefficients |
| Complete reference-expression parsing in native Symbolica | Passed: constants, supplied Laurent coefficients and uncertainty magnitudes parse as validation-only data |
| Native graph/cut/input failure paths | Passed: all nine integration tests, including all108 identities, graph complements, nonunit rational routing, consistent edge reversal, wrong incidence/charges, independent species, zero-power medium offsets and independent mass derivatives; see `input-regressions.json` |
| Guarded native dependency integration | Pending final command results from dependency owner |
| Compact seed distribution/threshold regressions | Pending final command results |
| Ordinary/cut/transport and boundary-owner regressions | Passed: 62 tests across nine targets, recorded in `guarded-regressions.json`; additional targeted gates remain pending |
| Python feature compilation | Passed; actual Python host import and runtime invocation were not run |
| CLI preparation smoke | Passed on the massive sunset and all three four-loop inputs; saved `*-preparation.json` reports explicitly mark `algebraic_preparation_only` and `numerical_evaluation.available=false` |

Pending rows must be replaced by measured results before any corresponding check
is reported as passing. The graph checker is an exact finite combinatorial
calculation, not a numerical multiloop benchmark. Its assertions include signed
incidence conservation, charge cycles, all retained connected cut complements,
nonzero routing determinants, vertex/edge deletion connectivity, exhaustive
vertex-permutation graph canonicalization and exact row-rank partition tests
excluding independent denominator blocks.

## Implemented capabilities and limits

The native input path accepts user-defined graph incidence, arbitrary exact
rational routings consistent with that incidence, independent line masses,
conserved nonbranching species charges, polynomial scalar and medium numerators,
uniform powers, Laurent orders and requested precision. Preparation uses native
Linnet graph connectivity and Symbolica affine basis conversion. It preserves
zero-power chemical assignments, independent mass coordinates and fixed physical
input expansion coefficients. Preparation reports are explicitly labelled as
algebraic; there is no benchmark-name evaluation dispatch. Loop charges are
currently integers and edge charges are unit signed values with at most one
species per edge. Rational routings requiring fractional loop charges after a
basis change are not supported; unrestricted arbitrary-routing coverage remains
incomplete.

The guarded RustRed integration keeps occupation semantics distinct from required
cuts, with guarded-source replay infrastructure. It does not demonstrate generic
four-loop weighted derivative closure. Compact polynomial shell seeds account
for upper-surface terms when raising a physical line. They handle below-support
zeros and admitted threshold limits or derived singularity diagnostics; they do
not evaluate a retained virtual denominator or supply arbitrary massless
infrared continuation. Their energy insertion is real E^r; an original Euclidean
q0^r numerator requires i^r on the shell q0=iE, independently of the full
amplitude normalization.

The selected massive sunset has a sufficient analytic admission for real
nonnegative fixed-shell auxiliary mass: its one-cut bubble Feynman-parameter
polynomial and fully occupied neutral denominator are positive throughout the
open domain |m1-m2| < M with positive masses. The physical point (1/2,1/2,1)
lies inside this domain and away from occupied thresholds. This supports
independent physical-mass derivatives and the inherited common contour for
that small graph; it does not establish generic four-loop contour admission or
provide an integrated boundary or amplitude value.

Boundary-owner changes provide an explicit tadpole-only recursive policy,
projected hard/soft factors with coordinate maps, and logarithmic region matching
with unknown coefficient/depth/rank failures. These are prerequisites for a
weighted evaluator, not integrated finite-density boundaries by themselves.

## Mandatory results still absent

| Required result | Missing work |
| --- | --- |
| Complete massive small nonfactorized graph | Evaluate and assemble its vacuum, mixed and fully occupied terms, including a raised charged line, medium numerator and nonzero surface term |
| Generic finite-density reduction | Derive/admit sufficient weighted source identities and establish target and differential closure with replay and nonzero guards |
| Occupied recursive boundaries | Integrate retained nonpolynomial soft factors through a terminating weighted recursion and establish their numerical accuracy |
| Common-contour continuation | Resolve branch/pinch admissions and retain the regulated complete-amplitude prescription through deformation and endpoint limits |
| Whole-amplitude normalization | Validate vacuum and occupied normalization together on a massive example with a nonzero vacuum contribution |
| Finite-density endpoint projection | Exercise rational target cancellation, uncancelled divergence and insufficient order through a complete finite-density evaluation path |
| Three distinct four-loop families | Produce predictions for the certified seven-, eight- and nine-edge graphs, including raised and polynomial targets and massless endpoints |
| Independent supplemental references | Generate and certify raised-occupied, massive, medium and multiple-chemical-potential references absent from the supplied oracle |
| Ten stable requested Laurent digits | Independently refine precision, boundary order/start scale and epsilon grid and record errors before comparison |

The mandatory supplied oracle subset is I37, I91 and I115, plus lower-loop L2_1.
I37 has symbolic supplied coefficients. I91's finite uncertainty is 1/1000;
I115's finite uncertainty is 1/10^10. These uncertainties constrain reference
agreement without weakening the separate ten-digit prediction-stability gate.
Higher unsupplied Laurent orders remain unspecified. The acceptance manifest
records each supplemental reference as required but not generated.

## Reproduction commands

From the repository root, with Python available for the definition-only checker:

```sh
python3 tools/finite_density/acceptance_manifest.py
nix develop --command cargo test --locked --release --test finite_density_input
nix develop --command cargo run --locked --release --bin rustflow -- \
  finite-density-prepare examples/finite_density/massive_two_loop_sunset.json
```

The CLI prints an algebraic preparation report. The Python host exposes
`prepare_finite_density(input_json, cut_budget=65536)` with the same JSON report.
Neither command produces a finite-density multiloop prediction.

Actual command logs, test counts, runtime and memory measurements belong beside
this report. Their measurement scope must state whether compilation, reduction
preparation and cache reads are included. No runtime or accuracy is inferred
from the existence of a parsed fixture or a bounded algebraic pilot.
