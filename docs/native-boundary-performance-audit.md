# Native boundary performance audit, 2026-10-05

This read-only audit concerns runtime `6e86afd`, the release-mode
`gg_hg_boundaries` runner, and the native common-mass PL48/NP61 calculations.
It does not establish performance parity with upstream codes.

The [recorded planar acceptance](../reports/validation/2026-10-05-native-planar-boundary.json)
passes all 240 canonical coefficients through epsilon power four at the
Euclidean anchor, with 20 verified digits and a reference accuracy cap of 40
digits. Reference data were loaded only after native generation. The successful
16-worker continuation took 396.010 seconds after a cancelled 4-worker run of
356.556 seconds; all 16 completed sample files from the cancelled run remained
byte-identical. This is a resumed timing, not a cold benchmark. A subsequent
verified-boundary restart took 0.451 seconds and evaluated no epsilon samples.

Subsequent [NP acceptance](../reports/validation/2026-10-05-native-nonplanar-boundary.json)
passed all 305 coefficients with 20 verified digits. Its resumed 16-worker run
took 1761.449 seconds after a 422.618-second interrupted 4-worker run; four
completed sample files remained byte-identical. The warm reload took 1.860
seconds. A [first physical planar source](../reports/validation/2026-10-05-native-planar-physical-source.json)
passed 30-digit native refinement in 1653.328 seconds with 16 workers, and all
240 coefficients agree with the independent source at its 24-digit evidence
cap. The sixteen physical configurations, transport, and coherent amplitude
remain separate acceptance work; see the
[native acceptance runner](gg-hg-native-boundaries.md) and
[coverage inventory](gg-hg-coverage-inventory.md).

## Observations

The first eight PL samples in the 4-worker run give these mean stage intervals.
Events from concurrent workers are interleaved; these are averages of matched
phase cohorts, not instrumented CPU profiles.

| PL48 phase | Mean seconds per sample |
|---|---:|
| Compile connection | 0.072 |
| Infinity expansion, including initial shared preparation | 11.679 |
| Recursive boundary construction and matching | 20.874 |
| Main transport | 40.814 |
| Endpoint expansion | 6.579 |

Exact infinity preparation ran once in 3.484 seconds. Later numerical infinity
and endpoint evaluations took approximately 9.8 and 5.8 seconds. Main PL
transport used 43 accepted Taylor steps. Its 12-master one-loop boundary child
used 103 steps and approximately 11.9 seconds per sample. The first three NP
main transports took approximately 218 seconds each; the fourth took 282
seconds. These NP observations precede full-fit validation.

The PL process used approximately 3.84 CPU cores and 548 MiB RSS with four
workers. A snapshot of the resumed 16-worker process showed 14.79 CPU cores,
2.86 GiB RSS, and 18 threads. This does not suggest sustained serialization by
the preparation cache. The first run wrote less than 2 MiB of progress/sample
data at the resource snapshot; numerical work dominated.

Local evidence: `target/gg-hg-planar-native-fit.log`,
`target/gg-hg-nonplanar-native-fit.log`, and the `progress.jsonl` files below
`target/gg-hg-native-acceptance/{pl20,np20}/runs/`. The durable planar report
above records run identifiers, build digests, accuracy, restart, and timing
scope; local progress logs are not committed scientific reference data.

## Sampling budget

The [fit controller](../src/engine.rs), `fit_samples_refined_leading`, uses the
generic two-loop pole envelope through epsilon power four, currently `-4..=4`.
It requires agreement of independent sample grids, precision, and order.
The following table uses the native runner's default 40 guard digits and initial
order 80. The notebook instead uses 60 guard digits and initial order 96, so its
30-digit seeds start at 90 working digits and refine to at least 110.

| Requested digits | First two grids | Working digits | Series orders | Samples if all four attempts run |
|---|---|---|---|---:|
| 20 | 31 + 35 = 66 | 60, 80 | 80, 96 | 148 |
| 30 | 37 + 41 = 78 | 70, 90 | 80, 96 | 172 |

Thus sixteen 30-digit configurations with this envelope require at least 1,248
finite-epsilon evaluations. Precision-error hints can increase working
precision further. Worker count is excluded from both
[sample checkpoint settings](../src/sample_checkpoint.rs) and
[symbolic system keys](../src/cache.rs), permitting the validated 4-to-16-worker
restart. Higher precision/order samples have distinct identities and cannot
substitute for independent validation.

## Candidates for measured follow-up

- [Frobenius evaluation](../src/frobenius.rs) still reparses exact rational
  entries and constructs the same recurrence matrix for each generalized
  eigenvector. At dimension 48 and order 80 this constructs 8.85 million matrix
  slots per endpoint evaluation. Cache sparse exact coefficient arrays and
  hoist matrices/resonance decisions per eigenvalue/order, retaining the
  existing arithmetic and coupled recurrence.
- [Recursive boundaries](../src/recursive.rs) repeat exact region enumeration,
  powers, factorization, and scaleless classification across samples. Those
  symbolic pieces can be cached by family, basis, deformation, and required
  order; numerical rank/admission decisions must remain sample-specific.
  Child-flow lookup/build/insert also permits duplicate concurrent misses:
  the initial PL log contains paired 3-master preparations near 35.500,
  36.114, and 36.834 seconds.
- [Laurent fitting](../src/epsilon.rs) reconstructs and factors the same sample
  matrix separately for every output. A shared factorization with multiple
  right-hand sides avoids this repetition without reducing sample evidence.
- [Taylor transport](../src/ode.rs) uses a pole-radius/3 proposal and halving.
  Its existing bracketed policy is an appropriate controlled comparison because
  every proposed step retains tail, endpoint/midpoint defect, whole-segment,
  conditioning, and domain checks. Current progress events do not expose
  rejected-predicate counts, so accepted-step counts alone do not justify
  changing tolerances or geometric limits.

These are source-supported candidates, not measured speedup claims. Retain the
guard digits, pole-cancellation checks, and independent refinement until any
alternative has equivalent accuracy evidence.
