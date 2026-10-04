# Symbolica root-refinement performance reproducer

Symbolica commit `75f8350094b90254ee71dc2a391fde0d14b0204a` isolates the two
roots of `(x - 1/2)^2 + 10^-100` in about 2 ms, but the first request to refine
an enclosure to radius `2^-232` exceeds a 20-second bound. The exact roots are
`1/2 ± i*10^-50`. This program uses only Symbolica, with exact rational inputs.

With a Rust toolchain and a configured Symbolica license:

```sh
cd repros/symbolica-close-root-refinement
cargo build --release
timeout 30s target/release/symbolica-close-root-refinement-mre 100 shifted 232
# Controls that completed quickly in the recorded experiment:
target/release/symbolica-close-root-refinement-mre 80 shifted 232
target/release/symbolica-close-root-refinement-mre 100 shifted 128
```

`minimal.rs` constructs the polynomial directly, calls `isolate_roots`, then
`IsolatedRoot::refined`, then classifies the roots. The first refinement is the
slow operation. Exit status 124 from `timeout` means the run reached its time
bound; it is not an error returned by Symbolica. No license material is included.

The recorded executable was built with `rustc` against RustFlow's validated
release Symbolica library. The portable manifest pins the same revision but has
not been independently rebuilt. `final-results.json` and the three `*-final.log`
files preserve the observations. `provenance.json` also records earlier
diagnostics and hashes of local profiling artifacts under `target/`; those
artifacts and binaries are not distributed with this reproducer.

Earlier bounded diagnostics also found slow refinement for `x^2+10^-100` at
radius `2^-448`. A profile reaches `refine_root_to_tolerance` and
`refine_root_disk_with_newton`, with substantial integer multiplication/GCD work
in the exact rational fallback. This suggests intermediate rational growth;
the reason the earlier bounded-center attempts fail remains undiagnosed.

These observations establish a performance problem, not nontermination or an
incorrect root certificate. A fix should preserve distinct conjugate roots,
exact certification, and the requested enclosure radius. The underlying native
code is in `src/poly/univariate/roots.rs` at the pinned Symbolica revision.
