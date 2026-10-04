# Full 13-integral canonical cache validation

The canonical-form cache passed the same two nearby physical-point checks as the
[dense-system experiment](full13-physical-cache.md), including every one of the
13 planar one-loop five-point basis integrals through epsilon order four. The
exact destinations are `PH1+(PH6-PH1)/1000` and `PH1+(PH6-PH1)/500`. Both paths
pass the existing real-radicand, sheet and source-domain checks.

The second request starts from the cached first destination. Sixteen verified
entries survive binary save/load, and an exact repeat performs no transport.
All 65 coefficients at each destination agree with independent higher-precision,
higher-order transport from PH1: the largest absolute differences are
`1.881e-80` and `2.898e-80`. The second endpoint also agrees with the completed
fresh original DiffExp calculation within `1.546e-36`, consistent with its
reported error estimate of about `5.62e-35`. This reuses the original partial-point
oracle from the dense experiment; it is not another original run.

PH1 is a supplied upstream boundary with recorded 132-digit absolute accuracy.
Every inexact component, including finite-accuracy zeros, carries a conservative
`1e-130` absolute error. The source advertises 128 digits; the two transported
points advertise 48 and 45 digits after propagating uncertainty. The much smaller
observed refinement differences do not imply 80 verified digits. The requested
accuracy is 20 digits, and the checks provide consistency/error estimates rather
than directed-interval certificates. Automatic AMF boundary construction, full
amplitudes and threshold-crossing cache reuse are outside this experiment.

The frozen standalone implementation uses `CanonicalAlgebraicSystem` and the
same `RustFlowCache` owner and binary codec as dense systems. It retains ordered
letters, matrices, registered roots and original source-domain guards. Its source
snapshot includes both the native exact-polynomial improvement and the shared
source-domain owner; comparisons with the earlier dense experiment cannot be
attributed to one change.

| Operation | Observed wall time |
| --- | ---: |
| Canonical-system construction | 0.1020 s |
| Cache identity construction | 0.00577 s |
| First query, including numerical refinement | 7.2768 s |
| Nearby second query from the first point | 3.3417 s |
| Direct PH1 to second point, matching numerical settings | 10.2585 s |
| Save / load 16 entries | 0.0435 / 0.0855 s |
| Exact repeated query after load | 0.00488 s |

Both the progressive and matching-settings direct requests use 20 requested
digits, 40 guard digits and initial series order 64. They have different path
lengths and source uncertainties: the direct request begins with 128 source
digits, while the progressive second request begins with 48. The matching direct
result differs from the progressive result by at most `7.411e-81`. Separate
validation runs request 40 digits with 50 guard digits and initial order 96.

The complete run, including those independent checks, took 67.03 seconds on
CPU 31 with one worker. Peak polled resident memory was 97.6 MB, and the binary
bank is 547,122 bytes. Other validation jobs occupied CPUs 28–30. These are
descriptive measurements, not amplitude throughput or cross-implementation
performance-parity claims. Further profiling was deferred in favor of broader
coverage.

The [report](../reports/performance/2026-10-04-full13-canonical-cache.json)
preserves every value and error estimate, exact coordinates, compiler commands,
source and binary hashes, settings, raw oracle results and resource limits.
Frozen raw artifacts remain under `target/fivepoint-canonical-cache`.
