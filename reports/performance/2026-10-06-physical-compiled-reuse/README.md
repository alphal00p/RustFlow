# Initial supplied-boundary physical transport

The algebraic physical route now retains the numerical compilation paired with
its final successful precision/order refinement and uses it for uncertainty
admission. Previously that same immutable source was compiled again at the same
precision after branch and checkpoint classification. Rational routes retain
their existing behavior. Cancellation is checked before uncertainty admission;
the uncertainty proofs and transactional cache admission are unchanged.

This is an object-lifetime change within one call. It adds no global cache,
threading, numerical algorithm, public option, dependency patch or persisted
format. It uses the existing native Symbolica/Numerica owners. The normal source
fingerprint changes, so old binary cache identities are still checked normally.

## Matched native measurements

Baseline commit: `699691c8c48717d4b1a0f3f637844a7d95656888`.
Both frozen release executables ran on CPU 43 of the shared AMD EPYC 9754 host,
with one runtime worker, requested accuracy 20 digits, 60 guard digits and starting
series order 96. Compilation and other workloads continued on the host. These
are one matched pair, not statistical evidence of a small improvement.

The inputs are the actual W permutation-2 planar (48 masters) and nonplanar
(61 masters) initial physical boundaries exported from the independently refined
native publication run. Their exact dyadic values, error arrays, 40-digit source
caps, original provenance and 415-bit storage were retained. The importer checks
the native system's exact start, dimension, range and root germ, then explicitly
admits them as supplied evidence under the current identity. It does not load a
foreign binary cache or import already transported destination values.

| Initial transport | Before (s) | Reused final object (s) |
|---|---:|---:|
| Planar | 8.0154 | 8.0442 |
| Nonplanar | 19.7090 | 19.5211 |

There is **no material speedup demonstrated by this pair**. Every non-timing
scientific output is exactly equal: all stored coefficient and error dyadics,
branch germs, cache coordinates and provenance, verified/input digits,
diagnostics, and complete saved Taylor-segment coefficient hashes. Each result
has one accepted step, no rejected steps, two inserted cache boundaries,
349 working bits/order 128 and 37 verified output digits capped by 40 input
digits. The 48 focused integration tests also pass, covering fresh preparation,
source exclusions, refinement, algebraic endpoints, branch-aware persistence,
source floors/caps, automatic error-proof retries and transactional cancellation.

## Measured remaining cost and a transport-only profile

A separate sampling run of the frozen candidate recorded 2,568 `cpu-clock:u`
samples at 99 Hz with no lost samples. It used 16 KiB DWARF callchains and
non-inlined inclusive symbol reports. Its complete scientific output is again
exactly equal to the unprofiled candidate. Inclusive samples include algebraic
local-chart construction (58.26%), exact-source residual-chart construction
(9.50%), conditioned segment evaluation (8.37%), and prepared numerical row
compilation (2.22%). These overlapping callchain percentages are not additive
stage timers. Self time is dominated by native MPFR/GMP arithmetic and allocation.

The same candidate was therefore tried with **requested accuracy still 20**, the
same 40-digit supplied source and all existing checks, but 40 initial guard
digits and starting series order 64. The solver still performs independent
precision/order refinement, residual and conditioning checks, uncertainty
propagation and source-capped admission.

| Initial transport | Guard 60/order 96 (s) | Guard 40/order 64 (s) | Reduction |
|---|---:|---:|---:|
| Planar | 8.0442 | 4.8545 | 39.65% |
| Nonplanar | 19.5211 | 12.2843 | 37.07% |

Both lower-profile results finish at 282 bits/order 96 with one accepted step,
no rejections and 37 verified output digits. Their source evidence, branch
germs, cache coordinates and caps are unchanged. Exact rational comparison
checks all 545 final complex coefficients and every stored cache boundary:
the Euclidean difference is within the sum of the two independently admitted
error bounds, and each profile's error and cross-profile difference meet the
requested 20-digit scaled tolerance. The comparison uses squared rational
inequalities, without conversion to machine floating point. Progressive cache
checkpoints correctly retain their own 37-digit input cap; the final boundary
retains the supplied 40-digit input cap.

This is one native run of each profile on two systems. No notebook defaults are
changed here. All 16 configurations and the final amplitude need their own
validation before adopting the smaller starting profile. These timings exclude
seed decoding/system setup and do not measure cold AMF boundary generation,
full notebook or browser/WASM execution.

## Reproduction and evidence

`report.json` records source/dependency digests, frozen executable hashes, input
hashes, timing and validation results. `evidence.tar.gz` includes raw outputs,
source manifests, the measured and final benchmark harnesses, comparison scripts,
compiler/check logs, the exact supplied seed export and textual sample reports.
The large native binaries and raw `perf.data` remain in the recorded local
evidence directory and are identified by hashes rather than committed.

The final harness adds guard/order environment controls and fixes two redundant
borrows found by strict Clippy. Those harness-only edits do not change production
source; the original measured harness is preserved as well. The benchmark example
is archived rather than installed as a public API.

With the native development environment and a valid Symbolica license, extract
the evidence, copy `initial_transport_profile.rs` into an isolated checkout's
`examples` directory, and decompress `boundaries.json.gz`.

```bash
cargo build --release --locked --example initial_transport_profile
taskset -c 43 target/release/examples/initial_transport_profile boundaries.json high.json
INITIAL_TRANSPORT_GUARD=40 INITIAL_TRANSPORT_ORDER=64 taskset -c 43 \
  target/release/examples/initial_transport_profile boundaries.json low.json
python3 compare_profiles.py high.json low.json comparison.json
cargo test --release --locked --test algebraic_cache --test algebraic_preparation \
  --test physical_transport --test prescribed_cache --test algebraic_endpoints \
  --test fundamental_boundary_errors
cargo clippy --release --locked --all-targets -- -D warnings
cargo clippy --release --locked --all-targets --features python_stubgen -- -D warnings
cargo fmt --check
```
