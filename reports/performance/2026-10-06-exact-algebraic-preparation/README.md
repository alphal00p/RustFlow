# Immutable exact preparation for registered-root transport

`AlgebraicSystem::prepare()` now retains an immutable exact snapshot of a pulled-back
system, its native root normalization and original source-domain guards. Regular
and prescribed physical routes prepare once and reuse this exact data for route
admission, both refinement profiles, and boundary uncertainty propagation.
Every numerical compilation remains fresh at the requested precision: coefficients,
source residuals, root equations and pole approximations are never promoted from
an earlier profile. See [the API notes](../../../docs/algebraic-preparation.md).

This is an extraction of existing Symbolica-owned algebra. It adds no independent
root, polynomial, matrix, graph or tensor algorithm and requires no dependency
patch. The official Symbolica community revision stays unchanged. Exact preparation
retains its source and normalized expressions in memory for its lifetime; it is
neither a global cache nor a serialized boundary record. Cancellation checks bracket
native phases; an individual native algebra call is not preempted.

## Measurements

These are medians of three alternating frozen baseline/candidate binary pairs on
CPU 35 of a shared cluster host. The baseline already contains the preceding
balanced native-series construction milestone. Other compiler/test jobs were
active on other CPUs. These compare library versions, not upstream DiffExp runtime.

| Supplied gg→hg system | Operation | Before (s) | Prepared (s) | Ratio |
|---|---|---:|---:|---:|
| Planar, 48 masters | Fixed precision integration | 1.040 | 1.038 | 1.001× |
| Planar, 48 masters | Checked first destination | 5.459 | 4.871 | 1.121× |
| Planar, 48 masters | Checked nearby reuse | 5.445 | 4.875 | 1.117× |
| Nonplanar, 61 masters | Fixed precision integration | 2.488 | 2.587 | 0.961× |
| Nonplanar, 61 masters | Checked first destination | 14.833 | 12.797 | 1.159× |
| Nonplanar, 61 masters | Checked nearby reuse | 14.903 | 12.737 | 1.170× |

Every non-timing output is exactly identical in all three pairs: complete fixed,
checked, nearby and independent direct-second values; uncertainty estimates;
source-capped verified digits; cache insertions; selected bits/order; and accepted
and rejected steps. The supplied source has a 24-digit accuracy cap. Checked first
and nearby destinations retain 21 and 20 verified digits, respectively, at 282 bits
and order 96, with one accepted step and no rejection. These are supplied-boundary
transports with the existing source uncertainty propagation and independent
precision/order comparison. They do not measure fresh boundary construction or
full amplitude assembly, and do not establish external-code performance parity.

A single fixed-precision integration cannot benefit from sharing preparation
across profiles; its unchanged/slightly slower measurements are retained explicitly.
The timing driver sets `ROOT_PROFILE_REPEATS=1` for each of six invocations. The
archived harness's JSON `settings.repeats` field remains hardcoded to 3; the actual
`runs` arrays each contain one sample, and the driver records the three alternating
pairs. The comparison script uses those actual samples, not that stale field.

## Validation and evidence

Four new public regressions check mutation-independent snapshots, both root sheets
in a coupled epsilon hierarchy, fresh high-precision roots, original domain holes,
validation, and typed cancellation. The affected integration gate passes 120 tests
across 20 binaries. Native and Python/stubgen strict all-target Clippy and formatting
checks pass. The explicit long gg→hg gate is recorded separately in `report.json`.
It retains all 4,360 transport/reference coefficients, 4,376 precanonical projection
rows and eight W/Z form-factor checks. Existing observable accuracy limits remain;
uniform 20-significant-digit accuracy is not newly claimed.

The archive retains the frozen source manifest, source/Cargo inputs, standalone
benchmark harness, raw output values/timings, binary checksums, and validation logs.
Executables and redundant binary caches remain local. An initial test-only compile
error used a nonexistent convenience method, then was fixed to use native Float
square root. A wider gate exposed a pre-existing cache-schema test expecting 5 while
this checkout writes 6; the test now reads the current version (at least 6) and rejects every older version. Both failure
logs and corrected passing logs are retained. An accidentally misdirected dependency
build was stopped before this crate compiled and restarted in the isolated target;
its log is retained separately and is not validation evidence.

This optimization changes neither cache serialization nor compatibility criteria.
The implementation source digest changes normally, so prior binary caches remain
subject to source-sensitive compatibility checks. Singular endpoints retain their
raw exact pullback and existing Frobenius owner. Public Python/CLI behavior does not
require new options.

## Reproduction

Use the native development environment and a valid Symbolica license. Extract
`evidence.tar.gz`; copy its harness example and support file into an isolated
checkout's `examples` directory. Build the baseline at
`3d29512a0021d94e50dc4daf07854c99040c7267`, and the candidate with this change:

```bash
cargo build --release --locked --example registered_root_profile
ROOT_PROFILE_REPEATS=1 taskset -c 35 \
  target/release/examples/registered_root_profile target/preparation-values
cargo test --release --locked --test algebraic_preparation
cargo test --release --locked --test gg_hg \
  complete_crossed_gg_hg_systems_refine_and_match_original -- --exact --ignored
cargo clippy --release --locked --all-targets -- -D warnings
cargo clippy --release --locked --all-targets --features python_stubgen -- -D warnings
```

Preserve each binary for the archived alternating driver. Use a private
`CARGO_TARGET_DIR` when sourcing an environment that sets a publication target.
