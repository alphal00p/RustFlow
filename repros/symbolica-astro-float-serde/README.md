# Astro Float serde roundtrip changes its value

This reproducer uses only the official Symbolica community revision
`6defcca968ca8411977fb1f641a9dee49ee7b7a7` and bincode. No RustFlow or HEPKit code,
scientific fixture, Python runtime, or local dependency patch is required.

```sh
cargo run --release --manifest-path repros/symbolica-astro-float-serde/Cargo.toml
```

A Float rounded from 1/3 with declared precision 415 roundtrips exactly through
Symbolica's `bincode::Encode`/`Decode`. Encoding that same object through its serde
implementation and decoding preserves the reported precision but changes its
exact rational value. The final assertion is expected to fail on the pinned owner.

The Astro backend's serde implementation writes `(prec, BigFloat.to_string())`
and reparses the decimal string. Its direct binary implementation uses hexadecimal
text. A separate scan of 192 nonzero scientific input components found 72 changed
serde values and zero changed direct-binary values. RustFlow cache schema 8 uses
the existing direct binary owner throughout its numerical payload; its outer
version, source compatibility, checksum and validation gates are retained.
