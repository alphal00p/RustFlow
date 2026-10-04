# Symbolica exact-root isolation performance reproducer

This second Symbolica-only reproducer is distinct from the `IsolatedRoot::refined` case. On Symbolica `75f8350094b90254ee71dc2a391fde0d14b0204a`, the exact real quadratic

```
delta = 16/10^43 + 1/10^68
P(x) = (x - 1/2)^2 + delta^2
```

has exact roots `1/2 ± i*delta`, but native `P.isolate_roots()` does not return within the recorded 20-second process bound. The complete contour regression that exposed it was terminated after 60 seconds. There is no subsequent refinement or classification in `minimal.rs`, and no RustFlow/RustRed/HEPKit dependency. Inputs use exact `Rational` coefficients; no numerical input conversion is involved.

With a configured Symbolica license and Rust toolchain:

```sh
cd repros/symbolica-close-root-isolation
cargo build --release
timeout 30s target/release/symbolica-close-root-isolation-mre
```

The Cargo manifest pins Symbolica and its companion crates. The actual recorded minimal program was built with a local standalone rustc script linked to the previously validated Symbolica rlib; that script and binary remain under `target/`, and the portable Cargo build was not separately rerun. No license is stored here.

The 20-second bound ended with signal 15 after the harness requested termination; the process had not returned an error or root result. This is evidence of severe performance, not proof of nontermination. Sampled RSS stayed near 11 MiB. The 5-second profile from the full contour case reached `isolate_square_free_roots_from` → `roots_hot_start`, with MPFR/GMP work. This differs from the prior refinement case's exact rational Newton fallback. Inspect the adaptive precision/hot-start path at the pinned `src/poly/univariate/roots.rs` around 2662–2755; no particular cause is asserted without further upstream instrumentation.

The physical geometry regression is preserved explicitly but ignored in the default contour suite because this synchronous native call is not cancellable. Its rounding-collapse assertion is also tested directly against exact analytic radius-zero disks for these known roots, so the geometry guard remains covered without relaxing or merging roots.

`result.json` records the bounded run, and `provenance.json` records hashes and profiling scope. No production or dependency source was changed.
