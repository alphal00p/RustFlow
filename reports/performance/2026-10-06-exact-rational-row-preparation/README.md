# Reusing exact native rational-row preparation

This change separates exact rational-row preparation from numerical compilation.
Symbolica continues to own rational conversion, denominator GCD/division/products,
factorization, coefficient extraction and root solving. Preparation retains exact
numerators, denominators, common-denominator rows and factor certificates; numerical
compilation uses fresh parameter values and working precision on every call.

The algebraic preparation object reuses these exact rows for its residuals, kernel
and root equations across the checked physical-transport profiles. The compatibility
rational compiler still prepares rows on each call; pure rational physical routes
and transformed local charts do not yet share those preparations across profiles.
There is no rounded-value cache, new numerical algorithm, dependency patch, public
option or cache-format change. See [the API notes](../../../docs/algebraic-preparation.md).

## Complete transport measurements

The baseline is commit `20a2af0b538de003ad32159dd47dd26f915afd5d`, which already
includes balanced native-series construction and immutable algebraic normalization.
The table shows final-source medians from three alternating frozen-binary pairs on CPU 35 of
the shared cluster host. Compiler/test activity continued on other CPUs. These
compare two versions of this library, not upstream DiffExp or AMFlow runtimes.

| Supplied gg→hg system | Operation | Baseline (s) | Prepared rows (s) | Ratio |
|---|---|---:|---:|---:|
| Planar, 48 masters | Fixed precision integration | 1.023 | 1.021 | 1.002× |
| Planar, 48 masters | Checked first destination | 4.775 | 4.600 | 1.038× |
| Planar, 48 masters | Checked nearby reuse | 4.772 | 4.613 | 1.035× |
| Nonplanar, 61 masters | Fixed precision integration | 2.496 | 2.472 | 1.010× |
| Nonplanar, 61 masters | Checked first destination | 12.391 | 11.466 | 1.081× |
| Nonplanar, 61 masters | Checked nearby reuse | 12.379 | 11.566 | 1.070× |

All non-timing outputs are exactly identical in all three pairs: complete fixed,
checked, nearby and independent direct-second values, propagated uncertainties,
source-capped verified digits, cache insertions, precision/order and step diagnostics.
The supplied source accuracy cap remains 24 digits; first and nearby destinations
report 21 and 20 verified digits at 282 bits/order 96 with one accepted step and
no rejection. Independent precision/order refinement and input uncertainty
propagation remain active.

The fixed integration kernel is unchanged. Its small timing variations, including
first-pair ratios below one, are retained as noise rather than claimed preparation
gains. The improvement applies to repeated checked compilation. It does not establish
upstream performance parity, automatic boundary-generation performance or full
amplitude runtime. The earlier diagnostic profile attributed 8.48% of inclusive
samples to rational-row compilation, motivating this limited follow-up.

Each timing invocation sets `ROOT_PROFILE_REPEATS=1`, yielding one actual run per
system. The copied benchmark harness still emits a hardcoded `settings.repeats=3`;
the archive records all six invocations and the actual run arrays, and the comparison
uses those arrays. Raw paired timings from both pre-layout and final-source
measurements are retained in `report.json` and the archive; the table uses only
the final-source pairs.

## Numerical and ownership checks

The new tests compare against a frozen copy of the previous row compiler, checking
complex coefficients, pole ordering, denominator factors and content, cleared rows
and exact source coefficients at 80, 201 and 400 bits. A separate check uses 700-bit
stored parameters to establish that exact source specialization keeps their full
stored dyadics, and that a later higher-precision compilation recomputes numerical
values rather than promoting low-precision values. Changed parameters, mutated
caller rows, zero specialization, excluded denominator domains and typed cancellation
are covered. All 37 affected private ODE tests passed on the pre-layout source, including mapped
source residuals, conditioning and the three new preparation tests.

Numerical operations retain their previous order. Factor first-occurrence ordering
is preserved even when the first numerical specialization makes a factor constant.
Numeric and exact-source sparse-zero decisions remain distinct. Each compiled
system keeps its original exact row expressions and the current caller's parameter
snapshot for coordinate pullbacks. Preparation retains exact expressions for its
lifetime and increases retained exact metadata; it is not a global cache. Native
algebra retains its existing resource limits and cancellation granularity.

The pre-layout source passed 138 integrations across 25 binaries and the complete
standalone gg→hg gate in 649.87 seconds. Strict Clippy then identified a large
`PreparedConnection` enum variant. The only subsequent source change boxes its
prepared algebraic object, adding one allocation per prepared physical route.
The final source passed 41 checks across 10 affected physical/cache test binaries and
all six checked benchmark invocations, followed by strict native and Python/stubgen
all-target Clippy. The standalone full gg→hg test directly uses compiled systems
and does not exercise `PreparedConnection`; its pre-layout result is retained
without claiming a second run for that layout change.

The archive preserves both source manifests, the pre-layout enum source, the Clippy
failure and correction, all numerical and formatting logs, standalone harness,
exact output comparisons and source/dependency checksums.
The full gg→hg gate retains all 4,360 reference/transport coefficients, 4,376
precanonical projection rows and eight W/Z form-factor checks. Existing observable
accuracy limitations remain: uniform 20-significant-digit form-factor accuracy is
not newly claimed. Gate counts and terminal outcomes are recorded in `report.json`.

No factorization, convolution, graph or tensor owner is reimplemented. The official
Symbolica community dependency remains unchanged. The normal implementation source
fingerprint changes, so old persisted boundaries still face the existing
source-sensitive compatibility check.

## Reproduction

Use the native development environment with a valid Symbolica license and a private
`CARGO_TARGET_DIR`. Extract `evidence.tar.gz`; copy the harness example and support
file to an isolated checkout's examples directory. Build the baseline commit above
and this candidate, retaining both executables for the archived alternating driver.

```bash
cargo build --release --locked --example registered_root_profile
ROOT_PROFILE_REPEATS=1 taskset -c 35 \
  target/release/examples/registered_root_profile target/rational-row-values
cargo test --release --locked --lib ode::
cargo test --release --locked --test gg_hg \
  complete_crossed_gg_hg_systems_refine_and_match_original -- --exact --ignored
cargo clippy --release --locked --all-targets -- -D warnings
cargo clippy --release --locked --all-targets --features python_stubgen -- -D warnings
```
