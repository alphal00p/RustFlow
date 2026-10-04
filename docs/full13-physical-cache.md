# Growing physical cache on the full 13-integral planar fixture

The full planar one-loop five-point fixture passed progressive cache reuse at two nearby physical points. Each request contains all 13 basis integrals through epsilon order 4: 65 complex coefficients. The exact destinations are `PH1 + (PH6-PH1)/1000` and `PH1 + (PH6-PH1)/500`. The native cache admitted both paths using its source-domain and real nonzero-radicand checks.

The starting values are supplied upstream PH1 reference data. Their recorded absolute accuracy is 132 digits; the test preserves a conservative `1e-130` error on every inexact input, including finite-accuracy zeros, and advertises only 128 source digits. It does not construct this boundary with AMF.

The second request selected the first destination as its starting point. The bank retained 16 verified entries, including independently checked intermediate points. Binary save/load preserved the entries, and repeating the second destination returned the cached value without transport. Both endpoints were also evaluated directly from the original PH1 source in a separate, more precise run.

| Check | Maximum absolute difference across all 65 coefficients |
| --- | --- |
| First cached endpoint versus refined direct-source transport | `1.881e-80` |
| Second cached endpoint versus refined direct-source transport | `2.898e-80` |
| Second cached endpoint versus fresh original DiffExp | `1.5454e-36` |

These are observed comparison differences. The cache conservatively advertises 48 digits at the first point and 45 at the second after propagating input uncertainty. The acceptance target is 20 digits. The error estimates use finite-precision refinement and weighted amplification, not directed-interval certification.

The fresh original calculation used the unchanged, hash-pinned DiffExp package and ancillary files, with working precision 150, expansion order 60, accuracy goal 25 and epsilon order 4. Its exact destination was the nearby `1/500` point, not the published PH6 endpoint. The reported total error estimate was about `5.62e-35`, consistent with the observed cross-implementation discrepancy. Boundary values, precision metadata, output values and source hashes are preserved in the [report](../reports/performance/2026-10-04-full13-physical-cache.json).

## Scope of the timings

These measurements use frozen standalone prototype modules on CPU 31 with one worker. Other independent validation runs occupied CPUs 28–30. They predate later shared source-domain and canonical-form cache changes; they are not current production throughput measurements.

| Operation | Observed wall time |
| --- | ---: |
| Dense physical canonical-system construction | 8.00 s |
| Cache identity and domain construction | 26.50 s |
| First cached query, with numerical refinement | 64.31 s |
| Second cached query, with numerical refinement | 52.03 s |
| Save 16 entries | 10.40 s |
| Load 16 entries | 36.28 s |
| Exact repeated query | 1.56 s |

The cached requests start with 60 working digits/order 64 and refine precision/order internally. The independent direct-source requests start with 90 working digits/order 96 and request 40 verified digits. Their different settings prevent a meaningful query speed ratio. Original DiffExp took 4.69 seconds for its transport and 8.01 seconds including process setup, with different settings and work; this is a correctness comparison, not a performance-parity claim.

The preserved baseline spent 424.02 seconds constructing the same cache identity and reached its 600-second limit before returning an endpoint. Profiling identified Symbolica's general-expression `AtomField` zero tests inside root-coefficient grouping. The candidate instead uses native exact `Q` polynomials and `MultivariatePolynomial::to_polynomial_in`, then applies the same root-power parity reduction. This reduced the completed identity phase to 26.50 seconds. The benchmark changes only that grouping operation; it keeps every domain condition and fallible quotient inversion.

Focused exact tests cover four inverse identities, a four-root sum, nonunit rejection, signed powers, reordered root declarations, an exact complex radicand and a `10^-1000` coefficient. Further query profiling still finds material repeated source-domain work. Canonical-form reuse is therefore a separate follow-up measurement.

The report embeds full numerical results, the fresh original result, source and binary hashes, compiler commands, resource limits, baseline evidence and the exact precision policy. Raw workspace artifacts remain under `target/fivepoint-cache-exact-coefficients`; the timed-out baseline remains under `target/fivepoint-cache-patch`. No amplitude, automatic-boundary or Monte Carlo throughput coverage is claimed.
