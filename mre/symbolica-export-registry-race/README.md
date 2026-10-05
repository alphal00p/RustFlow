# Symbolica export registry race

This standalone MRE depends only on official Symbolica community commit
`c3408e4ba1d3bdd4ea55678fad50e27009be13d4`. Numerica and Graphica are pinned to
that same commit. It does not depend on RustFlow or change the notebook runtime.

Run each mode in a separate process with a licensed native Rust toolchain:

```sh
cargo run --locked --release -- control
cargo run --locked --release -- finite
cargo run --locked --release -- polynomial
```

A custom `Write` implementation writes all bytes unchanged, then registers an
unrelated coefficient domain at a deterministic point in native serialization.
This legal reentrant callback reproduces an interleaving also possible from a
second thread. Export succeeds in all three modes, but the two mutation cases
produce streams which native import cannot read:

| Mode | Export length | Import result |
| --- | ---: | --- |
| control | 51 bytes | exact `1`, no trailing bytes |
| finite | 59 bytes | `UnexpectedEof` |
| polynomial | 64 bytes | `InvalidData` |

Every process also checks that another export/import succeeds after registration
has completed. The program deliberately asserts the baseline failure: exit zero
means the expected control success or mutation failure occurred. A fixed-owner
regression must instead assert successful import in all modes.

## Evidence and reproduction scope

[The measured baseline](baseline/baseline-verification.json) used the unchanged
`main.rs` and an already compiled official-c340 library. The exact direct-rustc
command, bounded-process runs, logs and exported streams are preserved under
`baseline/`, with file hashes in `baseline/baseline-bundle.json`. No executable,
build output directory or license is committed.

The portable manifest/lock at this directory were generated separately to pin
all three owner crates coherently; see [provenance](provenance.json). That fresh
Cargo build was not the route used for the recorded baseline. The original
manifest/lock are retained under `baseline/` for accurate provenance.

## Cause

`State::export_partial` reads polynomial/finite-field counts and constructs their
iterators separately. Its polynomial symbol-dependency collection also precedes
an independent list traversal. `State::export` has the same count/iterator issue.
Registry growth during the writer callback can therefore add payload entries
without updating the corresponding counts or dependencies.

The owner fix needs coherent bounded prefixes of all three append-only registries
without holding a lock across caller I/O. Holding the State lock across `Write`
would deadlock valid registration callbacks; taking that lock at export entry
also requires care because existing symbol-generator callbacks already hold it.
Symbols and polynomial lists can refer to each other, so capturing their lengths
in one fixed order alone does not establish a coherent dependency closure.
The binary layout need not change. The recorded baseline's suggested locking
approach is an initial proposal, not a validated fix. This inherited defect is
separate from RustFlow's reuse of native export blobs between calls.
