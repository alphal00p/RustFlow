# Symbolica native export registry race

This standalone Rust host has **only Symbolica** as a direct dependency, pinned
to official community commit `c3408e4ba1d3bdd4ea55678fad50e27009be13d4`.
It does not depend on RustFlow and changes no runtime pins.

The `Write` callback deterministically schedules unrelated registry registration
after its count has been serialized, before Symbolica constructs the iterator
for that registry. This is a valid reentrant writer and also demonstrates the
interleaving possible with registration on another thread. The callback writes
all requested bytes unchanged; it does not produce a custom serialization format.

Run each mode in a separate process in a licensed native Rust environment:

```sh
cargo run --locked --release -- control
cargo run --locked --release -- finite
cargo run --locked --release -- polynomial
```

Expected official-c340 results:

| Mode | Native export | Native import |
|---|---|---|
| control | `Ok`, 51 bytes | exact `1`, zero trailing bytes |
| finite | `Ok`, 59 bytes | `UnexpectedEof: failed to fill whole buffer` |
| polynomial | `Ok`, 64 bytes | `InvalidData: Cannot load legacy atom format...` |

Each mode additionally exports/imports after registration has finished; all
three post-registration controls succeed. Exit zero means the mode produced its
expected baseline outcome, including the deliberately asserted errors. To test
a fix, the race modes must instead require successful import and exact equality.

`baseline-verification.json` preserves the measured source/binary/stream hashes,
commands and logs. This session reused an already compiled official-c340 rlib
rather than rebuilding dependencies; the exact direct-rustc argument vector is
`build-command.json`, with compile-time `CARGO_CRATE_NAME` set to
`symbolica_export_registry_race`. No secret environment values are recorded.
`runs.json` records separate bounded processes (20 seconds, 2 GiB address space,
core dumps disabled); these are seconds-scale serialization checks, not physics
measurements.

The inherited defect is in native `State::export_partial`: a registry length and
its append-only iterator are acquired separately. Polynomial-list symbol discovery
also precedes a later independent traversal of the list registry. `State::export`
has the same count/iterator pattern. A correct owner fix snapshots coherent
registry prefixes under the existing State read lock before any caller I/O,
uses those prefixes for dependency discovery/counts/payload, and releases the
lock before invoking `Write` so reentrant registration cannot deadlock. There is
no proposed RustFlow serializer or new wire format.
