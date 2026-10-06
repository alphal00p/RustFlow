# Supplied-boundary notebook controller acceptance

The native notebook controller completed the calculation from empty caches using
the shipped sixteen starting configurations. Automatic boundary generation was
replaced with a failing test hook: this run cannot silently regenerate missing
or insufficient seeds. Numerical comparison files were read only after native
transport and amplitude assembly completed.

| Controller stage | Seconds |
|---|---:|
| Import supplied boundaries | 0.662 |
| Initial physical transport, all sixteen configurations | 101.775 |
| Initial amplitude preparation, projection and evaluation | 38.376 |
| Binary restart and exact repeated transport | 2.424 |
| Warm transport, all sixteen exact hits | 0.520 |
| Warm amplitude projection and evaluation | 2.681 |

These are one-core native measurements on CPU 44 of the shared AMD EPYC 9754
host, with other compilation active. Stage times include controller checkpoint
I/O but exclude Python/module/model startup and the notebook UI. The runtime is
a release cdylib from the shared community host, built with `python_stubgen`,
installed only in an isolated development environment. It is not a newly built
distribution wheel or a browser runtime.

The extension uses RustFlow `7096ba8`, Symbolica `6defcca`, HEPKit `b96600b` and
RustRed `7c1ed037`. In particular, it predates the separately validated
[staged amplitude optimization](../../performance/2026-10-06-amplitude-staged-contraction/README.md).
The 38-second amplitude stage above includes the original constructor; it does
not measure that newer optimization. The controller uses physical guard 30 and
order 32 with the adaptive accuracy checks intact.

All 4,360 transported coefficients passed the 20-digit mixed comparison against
the archived reference. The maximum scaled difference was approximately
`1.70e-30`. All eight W/Z form factors and three observables agreed within their
combined uncertainty allowances. The native EW square, interference and HEFT
square retain 35, 36 and 47 estimated relative digits. The independent archived
EW reference still supports only 19 comparison digits.

The resulting cache has 48 entries. Binary restart preserves every value,
comparison error, root sheet and stored provenance. Warm transport takes zero
ODE steps for all configurations, and warm amplitude values and errors equal
the initial outputs. Runtime/source attestation remained unchanged throughout
the run.

Focused checks on this same host passed 55 controller tests, 14 scientific-bundle
tests, six exporter tests, seven notebook smoke/static checks and Marimo lint.
An added provenance assertion initially expected an evaluated cache hit to omit
its deliberate history annotation. The corrected test checks raw stored
provenance exactly and requires returned history to retain the original
provenance; the failed log is retained. No numerical implementation changed.

`report.json` records the runtime and source hashes, comparisons, timings and
scope. `evidence.tar.gz` preserves the runner, exact steering sources, manifest,
logs and both source-sensitive native cache files. The portable input bundle is
the community data at commit `0744f89496cfa3d4e46e1149a45ee0b6bf8112cd`.
The archive contains no extension executable; its exact SHA256 is recorded in
the runtime evidence. Its source-sensitive binary caches cannot be imported
into a different build by ignoring compatibility checks.
