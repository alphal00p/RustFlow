# pflow: massless propagator integrals at fixed d0 from a one-mass flow (route B)

Module `src/pflow` (feature `native`). Front ends: `examples/pflow.rs` (flags), and `examples/dispersive.rs` with a class config
whose `"kind"` is `"propagator"`. Both call `pflow::cli::run`.

## Method
- Give one line e of a massless two-point family a mass: 1/(k_e^2 + z)^(a_e) with z = m^2/p^2. Then J(z) obeys the one-mass DE of
  RustRed (e-block: the masters with a_e >= 1).
- The discontinuity rho(s) = -(1/pi) Im J(-s + i0) obeys the homogeneous e-block. The z-free masters are z-independent sources and
  drop out.
- Where rho comes from:
  - on (0, 1): the non-integer Frobenius families at z = 0 (classes k h, h = d0/2, k = 1..L);
  - on (1, inf): the families at infinity (classes -k h);
  - their coefficients are fitted to the leading region rows (soft sets at 0, hard subgraphs at infinity; graph.rs).
- Two half-line patch-kernel transports take rho to the threshold neighbourhoods. There, local families at z = -1 (classes
  frac(j h), j = 0..2L) give the termwise integrals to s = 1.
- Boundary directions that the rows leave free are fixed by the continuity of the integer-class projection at z = -1.

## What the outputs are (E-normalisation)
E is Euclidean, measure d^dk / pi^(d/2), p^2 = 1. RustRed's I = (-1)^(sum a) E, with D_e = k^2 - z (Minkowski). For a master a with
a_e = index of the massive line, `value(m)` = FP int rho s^-m ds:
- m = 1: E[a] at z = 0, the massless value.
- m = n + 1 >= 2: (a_e)_n / n! x E[a + n 1_e]. The CLI also writes `integral` and `integral_value` = E[a + n 1_e].
- m = -k <= 0: (-1)^(a_e - 1) (a_e)_n / n! x E[a_e -> -n] with n = k + 1 - a_e. This is the e-pinched integral with line e as a
  numerator of power n; it vanishes when n < 0. `integral` / `integral_value` as above.
- `--extra a;b` gives non-master targets through their RustRed rule and the moments. The z-free masters in a rule need
  `--efree FILE`.

## Acceptance (every run)
All checks are gates with threshold 10^-T, T = `--target-digits` or by default min(w, provided digits) - gate_margin
(gate_margin = 20, `--gate-margin`):
- boundary fits: residual of every over-determined fit;
- order budget: wf - w - 10 must cover the digits lost per order at a spurious point, and the budget must not cost rank (counted
  from the first nonzero row);
- every nonzero row representable (n >= n_lo; a leading row within the family's orders). n_lo = min(--nlo, lowest nonzero
  row - 2) per class; the threshold families use min(--nlo, -10);
- threshold families span the e-block; the threshold decompositions; the integer-class continuity at z = -1. This is an
  accuracy estimate of the threshold join: measured at four loops it tracks the precision of the outputs that depend on the
  join (no1 at w45 vs w55). At small eps the near-degenerate classes frac(j h) cost digits (~25 at |eps| = 1.3e-3 at four
  loops with a census-top basis); the decomposition's condition number is reported but is a far too pessimistic bound;
- the continuity-fit residual. Directions that continuity cannot see are UNDETERMINED: their effect on every requested output
  must be below the gate, else the run is not accepted;
- every master has boundary rows (`--allow-skipped` to override). Every irreducible lower-loop graph is provided;
- entries of the provided table agree with the exact evaluator wherever it can reduce them;
- independent checks: exactly vanishing masters are numerical zeros; masters with a closed form (graph evaluator) agree with it.

A run that fails a gate still writes its output (`"accepted": false`, failures in `diag.acceptance`), and the front ends exit
with code 3. `--accept-diagnostics` makes the exit code 0. `digits_cap` = min(w, provided digits) bounds the digits a run may
claim; `digits_estimate` = min(digits_cap, -log10 continuity); `--two-precision K` estimates them independently.

## Lower-loop boundary values (`--provided FILE`)
```
{"schema": "pflow-provided.v1", "normalisation": "E, P^2=1", "d0": "p/q", "digits": D,
 "values": [{"edges": [[u, v, a], ...], "ext": [s, t], "value": "..."}]}
```
- Indices are integers or "p/q" strings; d0 must be the run's.
- Graphs are keyed by `graph::canonical`: index-0 lines and p-carrying bridges contracted, series lines merged, the minimum over
  relabellings.
- The exact evaluator runs first. An entry is used only for a graph it cannot reduce, and entries it can reduce are cross-checked.

## Checkpoints (`--workdir DIR`)
- Frobenius families write their exact state every 900 s. Every transport writes after every patch.
- Each precision pass uses its own subdirectory `p<w>-<wf>-<m>`. The family fingerprint is the DE content (d0, masters, every
  term).
- Rerunning the same command resumes.

## DE blocks
`--de` takes one RustRed harness solution, or several per-cut BLOCKS on one common master basis: a comma list, or a manifest
`{"blocks": [...]}`.
- The blocks are united.
- RustRed's identity rules (master = 1 x master) are ignored.
- The same target must have the same rule in every block.
- A master of one block must not be reduced by another block.

## Tests
`cargo test --release --lib pflow:: -- --include-ignored` (about 10 minutes at 2 threads without the four-loop regression; that
one takes 30-60 minutes).
