# Symbolica cross-factor root-isolation performance reproducer

At Symbolica commit `75f8350094b90254ee71dc2a391fde0d14b0204a`, calling
`isolate_roots()` on the exact degree-eight polynomial in `polynomial.txt`
exceeds 30 seconds. Using Symbolica's exact `factor()` first, then calling the same
`isolate_roots()` on each of its four factors (degrees 1, 2, 2, 3), returns all eight
roots in about 0.013 seconds in the recorded process. No RustFlow code is linked
into this reproducer.

With a Rust toolchain and a configured Symbolica license:

```sh
cargo build --release
timeout 30s target/release/symbolica-cross-factor-isolation-mre
target/release/symbolica-cross-factor-isolation-mre --factors
```

The two commands run in separate processes so cached factor certificates cannot mask the
first observation. Exit 124 means the external time bound was reached; Symbolica
does not itself return a failure here. This is a bounded performance observation,
not a proof of nontermination or an incorrect root certificate.

The polynomial is an exact denominator-domain condition pulled back along the
NP Higgs-jet Euclidean-to-W physical line. Its constant rational denominator is
removed because it is nonzero and has no roots. `minimal.rs` uses only native
Symbolica parsing, exact rational polynomial conversion, factorization and root
isolation. `--factors` preserves and prints factor multiplicities and all native
root counts; it does not approximate or merge roots.

The recorded binary was compiled by rustc against the already validated release
Symbolica library. The portable Cargo manifest pins the same revision but has
not been independently rebuilt. Source/input/library/binary/log hashes and
commands are in `provenance.json` and `results.json`. No license material is
included.

A separate 180-second full-planner replay also stopped within native root
isolation on this polynomial. Its short profile reached `isolate_roots`,
`RootCache::build_root_multiset`, `separate_isolated_roots` and
`refine_root_disk_with_newton`, with substantial GMP multiplication/GCD work.
The source's global separation loop refines roots from different defining
factors. Intermediate exact rational growth is a plausible explanation; this
reproducer does not establish the precise cause. No physical continuation or
numerical integral result is claimed by either mode.
