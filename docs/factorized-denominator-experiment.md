# Factorized denominator disk bounds

This experiment strengthens the local rational ODE defect certificate. At the
time of these measurements it ran as an isolated candidate. The implementation
has since been adopted into the publication branch; the measurements below
retain their original dependency and source provenance. They do not establish
full native notebook acceptance for the updated public dependency graph.

Symbolica factors each exact row denominator over its integer polynomial ring.
The complete product, including constant content and multiplicities, must
reconstruct the original denominator exactly. A compilation-local cache also
supplies factors to the existing pole scan. Existing exact Gaussian
specialization and native polynomial shifts produce the ball coefficients.

On each trial disk, each factor receives the existing directed triangle lower
bound. An inconclusive factor cannot certify a product. Positive bounds are
multiplied downward with their multiplicities. The implementation uses the
stronger of this product and the expanded polynomial bound, preserving the
expanded residual numerator, error budget, conditioning check, and source
domain. Denominator shifts and bounds are shared across epsilon channels.

For example, `(1-z)^10` on the disk of radius `1/3` has a strictly positive
factor bound although the expanded triangle bound is inconclusive. Conversely,
the expanded bound for `(1-z)(1+z)` at radius `1/2` is stronger and is retained.
Tests also cover true poles, two inconclusive factors, constant and complex
content, specialization, empty products, and distinct coupled epsilon rows.

## Captured nonplanar ODE replay

The following runs use the same captured native 61-component source and stored
dyadic starting values at epsilon `1/3700`, with 90 working digits, an unchanged
80-digit local gate, and a 2,000-predicate limit. The comparison threshold was
fixed at 50 mixed-scale digits before comparing outputs. All 61 components
passed exact rational squared-norm comparisons.

| Implementation | Series order | Accepted charts | Predicates | Transport seconds | Speed relative to baseline |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline | 96 | 263 | 1,008 | 911.483 | 1.00 |
| Candidate | 96 | 147 | 433 | 501.896 | 1.82 |
| Candidate | 160 | 69 | 133 | 385.971 | 2.36 |

Compilation and captured-boundary setup took 0.425, 0.407, and 0.381 seconds,
respectively. The maximum absolute / mixed differences from the baseline were
`7.30e-80 / 6.36e-86` at order 96 and `1.31e-79 / 8.51e-84` at order 160.
These are measured agreement values, not physical accuracy certificates.
The replay excludes infinity initialization, endpoint matching, and Laurent
fitting. Scientific jobs ran concurrently, so these are diagnostic timings,
not isolated comparisons against an upstream implementation.

The local replay inputs, exact output dyadics, probe build commands, and hashes
are recorded under `target/gg-hg-basis-work/factorized-denominator-*` and
`target/factorized-denominator-*-replay.*`. The accompanying validation report
preserves comparison results and provenance. The native runner and each linked
library were copied before execution to avoid artifact collisions with the
baseline worktree.

## Validation status

The candidate passed 137 release library tests (two ignored), including nine
new regressions, and strict release Clippy for the library and all tests.
All 397 integration tests across 89 targets also passed (ten ignored), including
the unequal-banana precision and independent-contour check. Production options
and the Python host remain unchanged.

## Fresh native physical boundary

The candidate generated `W-NP_EW1-1` from an empty numerical bank at 30 requested
digits, with 60 guard digits, order 96, 64 epsilon workers, and a 2,000-predicate
limit. Only three compatible exact IBP checkpoint files were copied; their
originals remained unchanged. All 78 evaluations across the independent
37- and 41-point grids completed, and the resulting boundary passed all 305
published source coefficients at the reference's recorded 24-digit cap.
The maximum mixed difference from that reference is `2.42e-31`.

Native generation took 1,731.489 seconds, compared with 2,476.472 seconds for
the baseline at the same settings (observed factor 1.43). The concurrent
workloads differed. Exact warm reload took 1.821 seconds and produced identical
values, errors and provenance without evaluating any samples.

The 305 candidate coefficients also passed a separate comparison with the
validated baseline at 30 mixed-scale digits and within their combined recorded
errors; the maximum mixed difference was `6.90e-81`. These checks validate one
complete nonplanar source boundary, not the full sixteen-source amplitude
workflow. No numerical reference boundary was supplied to the evaluator.

See the [native source report](../reports/validation/2026-10-05-native-nonplanar-physical-factorized-source.json),
[published-source comparison](../reports/validation/2026-10-05-native-nonplanar-physical-factorized-comparison.json),
and [baseline comparison](../reports/validation/2026-10-05-native-nonplanar-physical-factorized-baseline-comparison.json).
