# Dispersive mode and the fixed-precision patch kernel

The dispersive mode (`src/dispersive`, feature `native`) computes the value at `z = 1` of the master integrals of a fixed-dimension
auxiliary-mass differential equation from the **discontinuity across the tachyonic cut** of one line. It was built to lift
five-loop vacuum masters far beyond the accuracy that ordinary eta-flows reach in reasonable time, and it is the production path of
the grid-B samples of the five-loop vacuum campaign (class 16359).

## What it computes

Let `z = 1 + eta` be the squared mass of one line (the "e slot") and `D = d0` a fixed rational dimension `4 - 2 eps`. RustRed supplies
`dI/deta = A(eta) I` at that `d0`; the dispersive mode keeps the **e-block** (masters with a positive power of the e line; the e-free
masters drop out of the discontinuity, which obeys the homogeneous e-block system).

1. **Structure (exact, mod p).** Singular points (factor base of all denominators, Aberth roots), local Frobenius families at `z = 0`
   (the soft region, exponent `d0/2`) and at infinity (`-k d0/2`, `k = 1..L`), their ranks, pivot coordinates and stable orders.
   The linear-algebra Frobenius needs no normal form, so it works on non-Fuchsian Laporta bases of Poincaré rank 5.
2. **Soft boundary.** The soft family is fixed by the soft region's coefficients: Gamma factors times `(L-1)`-loop all-massive vacuum
   integrals, reduced by RustRed onto lower-loop masters whose values come from RustFlow samples (`--v native:...`). An optional
   **lane stage** (`--lane SPEC`) first regenerates the four-loop prism-lane values natively from four-loop P14 values.
3. **Transport** of the real discontinuity `d(s)` along `z = -s` from `s0` to `S` with the **patch kernel** (below): Taylor patches
   with polygon detours above every singular point on the axis; on real patches the state is projected to the real axis and the
   dropped imaginary part is recorded; the Cauchy weight `1/(1 - z)` accumulates `int d(s)/(s + 1) ds` patch by patch.
4. **Matching at `S`** to the infinity families (row-scaled least squares) and the **dispersion integral**
   `I_j = head [0, s0] + mid (patches) + tail [S, inf)`, `J_j(z = 1) = -(sin(pi d0/2)/pi) I_j`.

All polynomial algebra is exact (rug rationals and integers); values are MPFR at explicit working precisions: `W` for the transport and
`WF = W + guard` for families, fit, matching, head and tail. Pivoting of the MPFR families is the reference tool's threshold pivoting
(`threshold:30`, default); `complete` and `structural` are available.

## The patch kernel (`src/patch_kernel`)

A planned-path, fixed-precision transport of `dy/dx = A(x) y` with exact rational-function coefficients, separate from the adaptive
certified transport `CompiledSystem::transport`:

- each patch shifts every distinct polynomial **exactly** to its centre (integer synthetic division) and rounds once;
- the Taylor recurrence `D_g u_g = sum_t N_t y_j` per (row, denominator) group is scheduled over worker threads so that every value
  receives the same operations in the same order: **output is bit-identical for any thread count**;
- the state after every patch is written as an **exact checkpoint** (hexadecimal MPFR, fingerprint of system, path, precision, mode,
  Cauchy pole, projection, initial vector and kernel algorithm); a run never overwrites a checkpoint of another run or a more
  advanced one, and holds an advisory lock;
- non-finite values or a series tail above the bound stop the run with the state at the last good patch.

`TransportConfig::dispersive(bits)` is the dispersive configuration (real projection, Cauchy pole at `z = 1`);
`TransportConfig::generic(bits)` keeps complex states and no weight. Systems enter as `RationalSystem` (from entries, from Symbolica
atoms, from `DifferentialSystem`, or `CompiledSystem::patch_kernel_system`).

## Command line

`cargo run --release --example dispersive -- ...` (all options in `src/dispersive/cli.rs`):

```
dispersive --de DE.json --class-config CONFIG.json --m0 M0 --minf MINF --w W --wf WF --workdir DIR --out OUT.json
           --v native:SOFTPLAN.json,SOFTSOL.json,SAMPLE4L.json | native:SOFTPLAN.json,SOFTSOL.json,lane --lane LANESPEC.json
           [--two-precision dW,dWF] [--threads N] [--pivoting threshold:30] [--boundary-checks-min N]
           [--ref ORACLE_OUT.json]... [--code-version STR]
dispersive --progress DIR                       # one JSON progress line read from the checkpoint files
dispersive --system-info DE.json CONFIG.json    # e-block size and kernel system fingerprint
```

Everything class-specific comes from the DE file, the class config (e slot, loop number, named masters) and the command line.
`--two-precision dW,dWF` runs the sample a second time at `(W + dW, WF + dWF)` (`dW` may be negative) and writes per-entry digits.
Hidden diagnostic flags: `--stop-at-patch K`, `--lane-only 1`, `--family-check zero|inf-K`.

## Checkpoint and resume

Every stage checkpoints in the work directory: the zero and infinity families every 15 minutes between two Frobenius orders (exact
binary state), the transport after **every patch**, the infinity families and the final sample as stage files. Each file carries a
content fingerprint (DE terms, plan, boundary rows, configuration, kernel, code version); a work directory of another input or
configuration is refused, never reused. **Rerunning the same command resumes**, with any thread count. SIGTERM/SIGINT exit with code
75 (resumable): inside the transport after the kernel's next check, elsewhere at once (stage files are written atomically), so a shed
that allows only seconds of grace loses at most one patch or one family interval. A work-directory lock excludes a second process.

## Accuracy reporting

`out.json` (the campaign Store schema: `I` per e-block master, `J`, named masters) carries:

- `digits` per entry from the two-precision comparison (the measurement; agreement of two runs at different `W`, `WF`, different
  paths and orders), and `digits_summary`;
- `oracle_comparison` per `--ref` (per entry and named);
- the in-run gauges in `residuals`: soft-boundary fit, matching at `S`, Frobenius recursion residual of every family **relative to
  the family column's scale** and where the family is used (`s0` for the zero family, `w = 1/S` at infinity), transport tail and
  overlap, the imaginary part dropped after the detours; a fit or matching residual above `1e-10` stops the run;
- the boundary gate: full rank of the lower-loop endpoint system and its redundant equations at a numerical-zero bar (a consistency
  check of the inputs, not an accuracy statement); the gluing quality record of lower-loop samples is re-checked;
- `cost` per stage (wall, CPU seconds, peak RSS, summed over all resumed process segments).

The in-run gauges are not accuracy claims: final digits come from the two-precision comparison and from independent references.

## Validation (class 16359, node eps = 1/150000, d0 = 299999/75000; 14 Nov 2026)

569 e-block masters, 58 563 rational entries, 707 patches (288 complex), W = 165, WF = 195:

- Main run against the in-house reference sample at the same W: the five class masters agree to **59.2-63.8 digits**
  (P400 62.5, P441 63.8, P445 59.2, P551 63.4, P632 63.3); all 569 entries min 53.0, median 76.6. Against the reference's own W150 and
  W135 samples the agreement equals those samples' accuracy (46.9-55.6 and 32.6-47.3 digits).
- Residuals: soft fit 1e-111.6 (reference 1e-90.1), matching 1e-123.0 (reference 1e-101.4), recursion residuals of all six families
  1e-180..1e-191 on the column scale; detour imaginary part dropped 1e-82.3.
- Cost: 67.7 CPU-h, 4.6 GB peak (reference tool: 118.6 CPU-h, 13.5 GB); transport 59 CPU-h, 5.8 h wall at 16 threads.
- Controls: identical output for 1/2/4/8/16 threads (prism lane, full pipeline; five-loop transport segment); random SIGKILL/SIGTERM in
  every stage and resume with other thread counts; three memory-guard sheds through the launcher; all identical.
- **Pending:** the two-precision self-consistency run (W150/WF180) and a 1e-110 perturbation of all lower-loop inputs; this section
  will carry their numbers when they finish.
