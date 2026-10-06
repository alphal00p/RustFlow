# Bounded integer source residuals

This milestone adds opt-in `ResidualArithmetic::AdaptiveInteger`; ball arithmetic
remains the default. The implementation and its contracts are described in
[`docs/adaptive-integer-residuals.md`](../../../docs/adaptive-integer-residuals.md).
[`report.json`](report.json) records the measured complete-controller timings,
validation, implementation/dependency digests, and archive checksums.

The nearby fixed-boundary transports improve by 3.05–6.40 times in these runs.
The longer paths improve by 1.03–1.12 times because conservative ball fallback is
frequent. These are measurements against the existing Rust ball controller on a
shared, unpinned cluster host, not external-code performance parity. Exact DE
compilation, cache lookup, automatic boundary generation, registered-root gg→hg,
and amplitude assembly are excluded. The stored nearby boundary vectors have no
independent uncertainty cap, so 20-digit agreement/refinement is conditional on
those values and is not evidence of 20-digit physical-integral accuracy.

`evidence.tar.gz` contains three distinct records:

- `frozen-adaptive-default`: the earlier isolated candidate, its source, fixtures,
  three alternating timing repetitions, fresh precision/order comparisons, and
  the 95-test/19-target gate. That candidate selected adaptive arithmetic by
  default; this is not the final policy.
- `final-opt-in`: the final source manifest, changed source files, Cargo inputs,
  native and Python/stubgen strict Clippy logs, 41 passing ODE unit tests, 36
  passing public/CLI tests, and a final-source nearby/adverse/refinement run.
  Two benchmark tests are ignored during ordinary unit tests; one was explicitly
  executed for the final comparison. The final public tests separately witness
  option routing and independent epsilon refinement.
- `earlier-harness-failure`: a retained experiment whose first harness incorrectly
  assumed the old-ball degree-64 case must succeed at 201 bits. It correctly
  requests at least 244 bits. The corrected harness compares typed failures and
  also succeeds at 397 bits; production arithmetic was not changed to hide this.

The test-only timing wrapper selects each arithmetic owner explicitly through
the same continuation controller, independently of the default option. The
final routing tests verify that public `ball` and `adaptive_integer` options
select the intended owners for rational/epsilon and identity/Möbius charts.

## Reproduction

Use the repository's native development environment and a valid Symbolica
license. Extract into a new output directory so benchmark result files never
overwrite the retained archive:

```bash
mkdir -p target/adaptive-evidence
tar -xzf reports/performance/2026-10-06-adaptive-integer-residual/evidence.tar.gz \
    -C target/adaptive-evidence
cargo test --release --locked --lib ode::
cargo test --release --locked --test residual_arithmetic --test conditioning \
    --test continuation_budget --test diffexp --test mobius --test pade --bin rustflow
cargo clippy --release --locked --all-targets -- -D warnings
cargo clippy --release --locked --all-targets --features python,python_stubgen -- -D warnings
```

For the long-path constructor and complete-controller experiment:

```bash
ADAPTIVE_RESIDUAL_BENCH_DIR="$PWD/target/adaptive-evidence/frozen-adaptive-default/adaptive-benchmark" \
    cargo test --release --locked --lib \
    ode::integer_residual::bench::adaptive_residual_performance_and_full_transport \
    -- --exact --ignored --nocapture
```

Replace `adaptive-benchmark` with `adaptive-nearby` to measure nearby transports.
For fresh 201-bit/order-80 and 397-bit/order-112 comparisons plus adverse sources,
use the same environment variable with the test name
`ode::integer_residual::bench::adaptive_residual_refinement_and_adverse_cases`.
The exact stored input and endpoint arrays, source equations, mode counters,
precision requests, and individual timings are retained in the archive.

No community-host Python runtime rebuild is claimed by this milestone. The
Python constructor/getter and generated-stub integration were compiled with
strict all-target Clippy; existing host integration can consume the new option
when the library revision is updated.
