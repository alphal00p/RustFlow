# Native series construction and registered-root transport

This change builds fixed-precision Taylor series by balanced additions of native
Symbolica `Series` monomials. It rounds each supplied coefficient exactly as
before, skips exact rounded zeros, and retains the original absolute remainder.
The monomials have disjoint powers, so balancing does not reorder overlapping
floating-point sums. No series algebra, convolution, graph, or tensor owner is
reimplemented, and no Symbolica change is required.

The following are medians of three alternating baseline/candidate binary runs
on CPU 35 of the same shared cluster host. These measurements compare RustFlow
against its preceding construction, not against Mathematica/DiffExp. Other
compiler/test workloads were active on other CPUs; the host was not reserved.

| Supplied gg→hg system | Operation | Before (s) | After (s) | Ratio |
|---|---|---:|---:|---:|
| Planar, 48 masters | Fixed precision integration | 1.247 | 1.032 | 1.208× |
| Planar, 48 masters | Checked first destination | 6.160 | 5.464 | 1.127× |
| Planar, 48 masters | Checked nearby reuse | 6.209 | 5.425 | 1.145× |
| Nonplanar, 61 masters | Fixed precision integration | 2.902 | 2.502 | 1.160× |
| Nonplanar, 61 masters | Checked first destination | 15.963 | 14.770 | 1.081× |
| Nonplanar, 61 masters | Checked nearby reuse | 16.213 | 14.817 | 1.094× |

Every non-timing output agrees exactly between the frozen binaries: fixed,
checked, nearby and independent direct-second endpoint values; propagated
uncertainty estimates; source-capped verified digits; cache insertion counts;
and selected precision/order/step diagnostics. Checked transports retain the
existing independent precision/order comparison and boundary uncertainty
propagation. The supplied source cap is 24 digits; first and nearby destinations
report 21 and 20 verified digits respectively. Both select 282 working bits and
order 96 with one accepted step and no rejected step.

The initial non-interleaved batches are also retained. They measured the
nonplanar checked operation at 16.093 s before and 16.118 s after (0.998×), while
the alternating comparison measured 1.081×. This variation is a reason to keep
the raw timing samples and avoid a universal speedup claim. Full amplitude
runtime, fresh automatic boundary generation and external-code performance
parity are outside this measurement.

## Profiling and native ownership

The diagnostic instrumented run is separate from all reported uninstrumented
timings. Its inclusive perf samples attributed 50.63% to local chart building,
24.80% to algebraic compilation, and 9.18% to `fixed_series`. MPFR complex
multiplication/rounding and native solution convolution remain substantial
costs. These are overlapping inclusive samples, not additive fractions.
The perf record contains 28,408 samples and reports 12 lost samples.

At official Symbolica community commit
`c3408e4ba1d3bdd4ea55678fad50e27009be13d4`, native `Series` exposes monomials and
addition but no coefficient-vector/iterator constructor or polynomial-to-series
conversion. Balanced native additions therefore improve construction without a
local dependency patch. The first experiment kept zero monomials and regressed
sparse final-only series; that negative result is preserved. Filtering exact
zeros after the unchanged rounding operation removed that regression.
Builder-only dense order-64/96 speed ratios are approximately 6.4×/8.4× at
201 bits. They are not complete-transport ratios.

An independent follow-up can reuse immutable exact root normalization across
the repeated compilations of one pulled-back path. Fresh numerical coefficient,
pole, root and source compilation must still occur at each working precision.
That follow-up is not part of this change or these measurements.

## Validation and evidence

The archive records the unchanged baseline helper, complete candidate Rust
sources/Cargo inputs, exact source manifests, copied standalone benchmark
harnesses, raw JSON values/timings, build/test logs, temporary instrumentation,
perf symbol report and the earlier sparse-construction regression. Large
executables, raw `perf.data` and redundant binary cache files remain local;
checksums identify the profiling executable and raw recording.

The focused native-series tests compare the previous linear construction across
56 dense, sparse, leading-zero and all-zero cases, two working precisions and
orders up to 512. They check exact coefficients, precision, unknown remainder
and native products. The broader validation covers rational/algebraic transport,
epsilon hierarchies, conditioning, branches, Möbius/Padé charts, saved/cache
points and singular endpoint expansions. The separately executed full gg→hg
acceptance checks all 4,360 coefficients in 16 configurations under two
precision/order profiles, 4,376 precanonical projection rows and the eight W/Z
form factors. The separate native amplitude test also passes. These preserve
the existing projection accuracy limits; uniform 20-significant-digit accuracy
for every observable is not newly claimed.

See [`report.json`](report.json) for exact gate counts, source/dependency hashes
and the archive checksum. No public options, Python bindings or binary cache
schema change. The existing source-sensitive implementation fingerprint changes,
so old cache snapshots remain subject to the normal compatibility rejection.

## Reproduction

Use the native development environment and a valid Symbolica license. Extract
`evidence.tar.gz` into a fresh output directory. For the full-path harness, copy
`harness/examples/registered_root_profile.rs` and its `support` directory into
an isolated checkout's `examples` directory, then build and run:

```bash
cargo build --release --locked --example registered_root_profile
ROOT_PROFILE_REPEATS=3 taskset -c 35 \
  target/release/examples/registered_root_profile target/series-comparison-values
```

Build the baseline at commit `7e01651d1c50010ab9b140829db8835c9f9f5c09` and the
candidate with this change; retain both executable files for alternating runs.
The archived `series_builder_probe.rs` expects the baseline `src/fixed_series.rs`
and should be run from the baseline checkout. Its candidate constructor is
self-contained in the probe and uses only native `Series` operations.

```bash
cargo test --release --locked --lib fixed_series::tests::
cargo test --release --locked --test algebraic --test algebraic_cache \
  --test algebraic_endpoints --test canonical_algebraic --test complex_ode \
  --test conditioning --test diffexp --test epsilon_algebraic --test gg_hg \
  --test gg_hg_amplitude --test mobius --test mobius_cache --test pade
cargo test --release --locked --test gg_hg \
  complete_crossed_gg_hg_systems_refine_and_match_original -- --exact --ignored
cargo clippy --release --locked --all-targets -- -D warnings
cargo clippy --release --locked --all-targets --features python,python_stubgen -- -D warnings
```
