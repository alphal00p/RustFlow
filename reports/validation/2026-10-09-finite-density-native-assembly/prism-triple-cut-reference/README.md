# Independent raised-prism reference

This directory covers only the supplemental target with powers
`[1,2,1,1,1,1,1,1,1]` and fixed original shifted-Euclidean numerator
`g1_2^2+g1_3*g2_4` from `examples/finite_density/triangular_prism.json`, at mu=1.
It does not cover the other prism target `g2_3^3`. No native AMF prediction or
supplied oracle value was read by the reference generator.

The complete reference sums all three pair cuts and the triple cut. Vacuum
and singleton zeros use the separately documented common regulated dimensional
continuation. Pair contributions use the five virtual identities independently
verified by native exact source replay; the optional external origin of those
candidate identities is retained in `../prism-virtual-exact-reference/`.
Kira and Wolfram are not required to generate or compare these references.
The triple-cut calculation retains the original mass derivative, shell
Jacobian, and moving upper Fermi surface as 42 separate exact monomials.

The derivation is `docs/finite-density-prism-barnes-reference.md`. The exact
monomials, primary mathematical references, Gamma-contour continuation,
per-sample pole separation, finite-part removal of the analytic regulator,
and uniform Stirling-tail checks are retained. Failed diagonal-contour
optimization is preserved under `diagnostics/failed-diagonal-contour-translation/`
and excluded from reference acceptance.

## Recorded refinements

Complete finite-epsilon refinement at epsilon=1/101 passes to about 7.11e-17
relative; a separately chosen rational analytic-index direction changes the
result by about 2.99e-34 relative. The epsilon=1/10000 refined pilot changes by
2.79e-26 relative. These are empirical checks, not rigorous error enclosures.

The Laurent controller fits epsilon^6 times the saved amplitude using exact
rational interpolation weights, retaining every coefficient without forcing
small values to zero. It requests orders [-4,0], while keeping the additional
working-pole coefficients as diagnostics. All five complete profiles are:

| Profile | Digits | Epsilon samples | Barnes step / cutoff |
| --- | ---: | --- | --- |
| `laurent-baseline` | 70 | i/10000, i=1..12 | 1/30 / 14 |
| `laurent-nearby-grid` | 70 | i/500000, i=1..12 | 1/30 / 14 |
| `laurent-fine-baseline` | 70 | i/1000000, i=1..12 | 1/30 / 14 |
| `laurent-fine-precision` | 80 | i/1000000, i=1..12 | 1/30 / 14 |
| `laurent-fine-quadrature` | 80 | i/1000000, i=1..12 | 1/40 / 18 |

All five profiles passed. `reference.json` is the published independent
reference; `laurent-reference-validation.json` records all 20 passing raw
coefficient refinement comparisons. `laurent-refinement-pending.json` is a
preserved historical partial-run report and is superseded by the final report.
The complete gate also checks ten-versus-twelve-node interpolation and the
imaginary part of the finest coefficients. Primary acceptance uses the raw
Euclidean coefficients: relative change <=1e-12, or absolute change <=1e-25
when both magnitudes are below 1e-20. Rescaled coefficients are retained with
their actual empirical differences; they are not independently rounded to
exact zero.

Reproduce the saved-output comparison and reference export:

```sh
python3 tools/finite_density/compare_prism_reference_grids.py \
  reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference \
  --output reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference/laurent-reference-validation.json \
  --reference-output reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference/reference.json
```

`generator/` contains the standalone reference implementation. It uses native
arithmetic utilities and Symbolica Gamma functions, with no AMF computation.
`standalone-linkage.json` records its actual source/library linkage. Each
sample records the executable and exact contour-input hashes before launch.
Contour input archives are checked byte for byte after decompression.
Whole-controller resource records include its executable-hashing memory;
per-sample resource files report the numerical process separately. Compilation
and Nix startup are excluded from timed numerical runs.

No full prism native closure, native numerical acceptance, or oracle comparison
is implied by any reference result in this directory.

The nearby epsilon grid was added because the initial coarse-to-fine change in the rescaled finite coefficient, 1.43e-12, slightly exceeds the planned 1e-12 absolute reference uncertainty. Final uncertainty is assessed from the nearby/fine grid, fine precision, fine quadrature, and ten/twelve-node changes; the larger coarse-grid change is retained as a convergence diagnostic. No tolerance was relaxed.


The normalized residue is `4.33925042047248508077` (empirical absolute
uncertainty `2.88e-21`) and finite coefficient `17.03999268944775288649`
(uncertainty `5.03e-16`). The raw finite coefficient is
`4.32431526349508047544e-8` (uncertainty `1.66e-24`). Small higher-pole
coefficients remain as fitted values in the artifact. These uncertainties
measure observed refinement changes, not rigorous integration-error bounds.
The final gate verifies 60 sample pole/tail checks, 2520 original-monomial
regulator-residue cancellations and 120 lossless contour archive round trips.
`generation-resources.json` aggregates the actual completed controller and
numerical-process resource records without including compilation or Nix startup.


`evidence-replay.json` checks all 120 exact pole/tail reports by regenerating
them from the 60 archived numerical contours. Their epsilon and contour hashes
are bound to the actual numerical command provenance. This closes a missing
explicit comparator binding; no numerical value or tolerance was changed.
