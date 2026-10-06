# Isolated owner fix verification

The fixed-owner report records commit
`e771a7b92c3f63c3bb49089f03ed4da108598a99`, based on official Symbolica community
`c3408e4ba1d3bdd4ea55678fad50e27009be13d4`. It was recorded before publication.
Only `src/state.rs` and `tests/import_export.rs` change in that owner commit.

The fix captures six registry offsets/counts twice and retries until both
collections agree. The append-only lengths are monotonic; agreement establishes
an overlapping stable interval. Full and partial export use those same bounded
prefixes for dependency discovery, serialized counts and entries. Symbol user
data and polynomial lists can depend on each other, so merely reading their
lengths in a fixed order would not suffice.

No State lock is acquired during snapshot capture or held across writer calls.
The original read-lock proposal was rejected: public symbol-generator callbacks
already hold the State write lock and may export an existing expression.
`generator-repro.rs` and its logs preserve the official control, rejected
candidate timeout and corrected successful run.

The three corrected standalone modes all export the same 51 bytes as the
callback-free baseline and import exact `1` without trailing bytes. The owner
integration test also covers six full/partial writer callbacks and export from
a symbol generator, alongside its existing import/export checks. Strict Clippy
checks passed for the selected owner library and the isolated test host.

[The validation report](patched-verification.json) records commands, hashes,
results and limits. `patched-host/` is the measured host, with explicit local
owner paths; it is evidence, not a second installed runtime. No notebook runtime,
source pin or acceptance bank was changed. This is focused serialization
validation, not a full Symbolica feature-suite or scientific benchmark run.

The snapshot covers published prefixes, not atomic publication of a whole
multi-symbol group. Concurrent unsafe State reset remains excluded by the
existing contract. Continuous registration can cause snapshot retries; no
bounded-wait performance guarantee is claimed.
