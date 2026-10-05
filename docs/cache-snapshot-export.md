# Native Atom export reuse during cache saves

`RustFlowCache::save` memoizes Symbolica's native `Atom::export` result within
one snapshot. The same memo covers both system representations, basis and
normalization data, continuation prescriptions, coordinate values and root
germs. Repeated zero entries in dense storage of sparse canonical matrices no
longer repeat Symbolica's state-export work.

The stored data still contains one self-contained native export for each Atom.
The boundary-cache wire format, version, payload checksum, source fingerprints,
validation, numerical precision and atomic replacement protocol are unchanged.
The memo is discarded before encoding the completed snapshot; it does not
survive saves or grow with the lifetime of a cache. It uses additional temporary
memory proportional to the distinct expressions and their native exports.
Symbolica owns expression equality, hashing, export and import.

A native export can include registered polynomial and finite-field state in
addition to the expression's symbols. If unrelated registration happens between
native export calls within a save, a repeated expression can retain an earlier
self-contained export.
Byte identity with the uncached writer is therefore checked under stable
registration state; semantic roundtrips are checked after registry growth.
The uncached writer itself does not promise identical bytes across registry
changes. This optimization does not add synchronization inside Symbolica's
exporter: official `c3408e4` reads polynomial/finite-field registry counts and
iterators separately, so concurrent mutation of those registries during one
native export is outside this safety claim. The concurrent regression registers
only unrelated plain symbols. No native locking fix, new serialization format,
or custom symbolic decoder is introduced.

This source change still changes the implementation fingerprint. Unchanged wire
format does not permit importing cache files produced by a different source
fingerprint; normal incompatibility rejection continues to apply.

## Reproducing the checks

In the ordinary native development environment, with the license supplied
privately, use a separate target directory:

```sh
export CARGO_TARGET_DIR="$PWD/target"
export CARGO_BUILD_JOBS=8
cargo fmt --all -- --check
cargo test --locked --release --lib transport_cache::serialization_tests -- --nocapture
cargo test --locked --release --test transport_cache --test canonical_cache --test algebraic_cache --test prescribed_cache --test complex_cache
cargo clippy --locked --release --all-targets -- -D warnings
cargo test --locked --release --lib transport_cache::serialization_tests::profile_canonical_snapshot_exports -- --ignored --exact --nocapture
```

The focused unit checks run in subprocesses so unrelated parallel tests cannot
change Symbolica's registration state during byte comparisons. They cover dense
and canonical identities, physical prescriptions, all coordinate origins, root
germs, complex MPFR values, uncertainty/provenance metadata, a fresh-process
import with different symbol allocation, mixed old/new polynomial and field
blobs imported after pre-registering different lists/fields in a fresh process,
concurrent unrelated symbol registration, and failed-save cleanup/preservation.
The existing integration checks retain corruption and incompatible-source
rejection tests.

## Timing scope

The ignored profile builds the actual 48/61 canonical systems and exact physical
source coordinates. Its 48 cache entries contain explicitly synthetic numerical
values and synthetic evidence. It reads no numerical reference or acceptance
bank. Every timed snapshot performs the ordinary boundary validation and
bincode payload encoding. Five uncached/memoized pairs alternate execution order
in one process after warmup, with byte equality checked outside the timed region.
A second profile registers 1,000 unrelated symbols to measure sensitivity to
Symbolica's process state. Atomic production saves and a checked production
reload are measured separately.

These measurements describe serialization, not transport, amplitude computation,
or the full notebook's process state. They do not establish notebook speedups or
performance parity with upstream solvers. The frozen notebook environment and
its acceptance banks were not modified for this measurement.

The 2026-10-05 release profile produced byte-identical 22,472,150-byte payloads
for all ten paired comparisons. Of 362,019 matrix entries, 357,430 were zero;
605 distinct Atoms were exported across the entire snapshot.

| Registered symbols in isolated process | Uncached payload, median | Memoized payload, median | Ratio |
| --- | ---: | ---: | ---: |
| 332, fixture setup | 1.283 s | 0.206 s | 6.23× |
| 1,332, including 1,000 extra unrelated symbols | 4.323 s | 0.211 s | 20.47× |

Three complete memoized atomic saves after registry growth took 0.289–0.315 s
(median 0.290 s). Reload through the production compatibility/validation path
preserved the payload. The raw timings, source digests, toolchain, and validation
record are in [the profile report](../reports/validation/2026-10-05-cache-snapshot-export.json).
