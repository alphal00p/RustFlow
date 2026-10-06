# Unchanged-engine probe of epsilon-singular physical transport

This experiment changes no production source or dependency. It tests the existing
ODE and physical boundary-error owners at exact nonzero epsilon values, before
implementing a sampled physical-transport fallback.

The supplied connection is `I' = N I / epsilon`, with
`N = [[1,1],[-1,-1]]`, `N²=0` and `I(0)=(1,0)`. Its exact solution is
`I(x)=(1+x/epsilon,-x/epsilon)`. Diagonal epsilon shearing rejects the negative
constraint cycle. An exact Symbolica basis transformation
`T=[[1/epsilon,1],[-1/epsilon,0]]` gives the regular matrix `[[0,1],[0,0]]`,
proving this fixture has a meromorphic solution despite the diagonal failure.
The probe checks this identity but does not apply that gauge during integration.

Every fixed-epsilon system is passed to the unchanged ordinary one-channel engine.
The full physical route separately uses a supplied boundary with **30 verified
digits**, explicit component errors **1e-35**, and **90-digit storage**. The default
options request 20 digits. The actual source-error floor and all existing physical
uncertainty checks remain active; an exact-boundary shortcut is not used.

| Epsilon | Distance | log10 of existing amplification estimate | Exact transfer infinity norm | Full physical result |
|---|---:|---:|---:|---|
| 1/100 | 1 | 86.859 | 201 | Accuracy rejection |
| 1/2500 | 1 | 2171.472 | 5001 | Accuracy rejection |
| 1/5800 | 1 | 5037.816 | 11601 | Accuracy rejection |
| 1/100 | 1/1000 | 0.0869 | 6/5 | Accepted, 27 verified digits |
| 1/2500 | 1/1000 | 2.1715 | 6 | Accepted, 26 verified digits |
| 1/5800 | 1/1000 | 5.0378 | 63/5 | Accepted, 23 verified digits |

These epsilon values are representative of the existing fitter: through the finite
term with leading power -1 and a 20-digit request, its first two grids have
25 samples up to 1/100 and 29 samples up to 1/200; the smallest values are
1/2500 and 1/5800. This probe is not a completed Laurent fit.

All six central integrations pass independent precision/order refinement in one
accepted step with zero rejections, using 282 bits and order 112. Their endpoint
values agree with the exact solution at 50 digits. That lower-level result checks
integration of the stated central boundary only. It is not a replacement for
physical input-uncertainty propagation.

All three long physical routes fail specifically because the supplied uncertainty
is amplified beyond the requested accuracy, and their caches remain unchanged.
All nearby routes pass while retaining that uncertainty. For this known exact
solution, the true transfer is `I+x*N/epsilon` and its infinity norm grows only as
`1+2*x/epsilon`. The report records the corresponding propagated source-error
comparison, but it is **not supplied to the numerical controller**. The inherited
exponential matrix-norm estimate is the cause of the long-route rejections; this
is useful conservative behavior until a stronger general owner is validated.

A generic finite-epsilon fallback therefore needs both valid finite-epsilon source
values and improved conditioning evidence for difficult bases. A finite cached
Laurent prefix cannot initialize it without a remainder bound, and stable fits
cannot by themselves establish meromorphy (compare `I'=I/epsilon`). The independent
design audit is included in the archive. No generic fallback, tighter bound,
performance parity or complete physics calculation is claimed here.

The archive retains unchanged production source, exact harness, raw JSON, source
and dependency fingerprints, build/Clippy/format logs and the initial design audit.
Runtime fields are observational only; this is a conditioning probe, not a timing
comparison. Reproduce in a private native target with the official dependency graph:

```bash
cargo run --release --locked --example epsilon_singular_probe
cargo clippy --release --locked --example epsilon_singular_probe -- -D warnings
cargo fmt --all -- --check
```
