# Smaller initial series for the live Higgs-jet demonstration

Initial guard 30/order 32 reduces all sixteen checked native physical transports
to **91.512 seconds**, compared with **270.213 seconds** for the original guard
60/order 96: **66.13% less time, or 2.953 times faster**. Requested accuracy stays
at 20 digits. Adaptive precision/order refinement, residual checks, conditioning,
supplied-error propagation and cache admission are unchanged.

This follows the [guard 40/order 64 comparison](../2026-10-06-initial-physical-transport/README.md),
which took 163.129 seconds. The same frozen executable, exact supplied inputs,
two systems, shared growing cache, projector and native HEPKit amplitude kernel
are used. All three profiles retain independently generated 40-digit starting
boundaries; no destination or amplitude values enter the calculation.

| Component | Original 60/96 | Compact 30/32 |
|---|---:|---:|
| Sixteen physical transports | 270.213 s | 91.512 s |
| W/Z form-factor evaluation | 1.688 s | 1.677 s |
| Three observable evaluations | 0.873 s | 0.864 s |
| Initial systems and supplied seeds | 0.694 s | 0.689 s |
| Initial projector construction | 1.312 s | 1.288 s |
| Initial amplitude kernel construction | 33.635 s | 34.338 s |

These are native release component measurements, excluding notebook startup,
Python/UI overhead, downloads and checkpoint I/O. They are not WASM or upstream
code timings. Each profile was run once on the same shared AMD EPYC 9754 host,
with the original on logical CPU 43 and this follow-up on CPU 44. Other builds
continued. The substantial first-use amplitude construction is separate from
the transport timing.

The unchanged exact-rational comparison checks all **4,360** final complex
coefficients, every progressively inserted cache boundary, eight form factors
and three coherent observables. Every difference lies within combined admitted
errors, and both profiles satisfy the requested 20-digit tolerances. There is
no binary64 conversion in these scientific inequalities. Starting coordinates,
root germs, source caps and provenance match. The cache grows from 16 to 48
entries. All sixteen transports finish in one accepted step, without rejections,
at 249 working bits/order 64, with 37 verified digits and input cap 40.

A second, independently written rational checker also passes using stricter
L1 differences and conservative maximum-component scales. It checks the 4,360
final coefficients, 8,720 inserted cache coefficients, eight factors and three
observables; its script and results are included in the archive.

The EW square, interference and HEFT square retain **35, 36 and 47** estimated
relative digits. The separate archived EW reference supports 19 comparison
digits; this run does not upgrade its precision. The lower budget applies to
physical transport; automatic auxiliary-mass boundary generation is unchanged.

`report.json` retains the complete comparison and diagnostics.
`evidence.tar.gz` contains the new exact outputs, run log, unchanged comparison
script and its output, plus a manifest. The baseline outputs, harness and source
provenance are archived with the previous report. Reproduce this follow-up with
that frozen harness and the same input:

```bash
INITIAL_TRANSPORT_GUARD=30 INITIAL_TRANSPORT_ORDER=32 taskset -c 44 \
  ./full16-native boundaries.json compact.json
python3 compare_full_profiles.py full16-high.json compact.json comparison.json
```

Source digest:
`a8e17584df583d823059db180f72a08fd18082d057578b16a3b3e4348bad76a2`.
Dependency digest:
`9ee6eced2bc39b05bb333b1b021ab556d2a1bb790a255205d3f141fc64cc96ff`.
