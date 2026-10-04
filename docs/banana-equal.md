# Equal-mass banana: analytic infinity boundary

The full four-master equal-mass three-loop banana system from the pinned
DiffExp example has been evaluated through ε⁴ at `t=-1` and `t=32`, starting
from the analytic gamma-function boundary at infinity. All 20 coefficients
agree with a live original-DiffExp calculation to 20 absolute decimal digits.
This validates the differential equation and boundary workflow; it is not a
new IBP reduction or full amplitude calculation.

The exact system and comparison data are in
[`banana-equal.json`](../fixtures/diffexp/banana-equal.json). The native
prototype and original solver measurements are recorded in
[`banana-equal-validation.json`](../reports/diffexp/banana-equal-validation.json)
and [`banana-equal-original.json`](../reports/diffexp/banana-equal-original.json).
Both use DiffExp commit `784c8229bf92369a03f011a48e161522c8c54bbd` as the
scientific-data reference. Source hashes and original precision metadata are
retained. DiffExp's implementation is loaded unchanged; its computational
source is not incorporated into this MIT library.

## Boundary and continuation

The dimension is `D=2-2ε`; the variable is `t=p²/m²`, with finite differential
singularities at 0, 4 and 16. Infinity is parameterized by `t=-1/x`, approaching
`x=0` on the positive real axis. Only master B3's leading `x` coefficient,
including its logarithms, and the constant master B4 are supplied. Masters B1
and B2 remain unknown at matching. The fixture also records their leading
coefficients deduced independently from the exact DE, but the Rust evaluation
does not use those coefficients as input.

The boundary follows from the four analytic gamma-function regions recorded
in the fixture. Expanding `log Γ(1+z)` cancels Euler's constant exactly. Expanding
to ε⁶ before multiplying the ε⁻² prefactor gives all required coefficients
through ε⁴, with powers of `log(x)` up to six. The four-master epsilon hierarchy
is represented as a 20-component exact rational system. The existing
generalized Frobenius recurrence constructs its complete basis; partial
coefficient constraints fix the integration constants. No numerical reference
value seeds the integration.

The native calculation evaluates the matched series at `x=1/128`, continues
to `x=1` (`t=-1`), then uses the existing pole-avoiding contour in the upper
half of the `t` plane to reach 32. This implements the example's `+i0`
prescription. Working precision and series order are increased together:
60 digits/order 80 and 80 digits/order 112. A preliminary order-20 diagnostic
failed the Taylor-tail accuracy check; accepted runs retain that check.

| Endpoint | Native refinement: maximum absolute change | Refined native vs original |
|---|---:|---:|
| `t=-1` | 3.46 × 10⁻⁵⁸ | 2.14 × 10⁻⁵⁰ |
| `t=32` | 1.38 × 10⁻⁵⁸ | 4.52 × 10⁻⁴⁶ |

The original calculation uses its notebook settings (500-digit working
precision/order 50, Möbius maps and Padé), then a 600-digit/order-75 refinement.
Its corresponding maximum changes are 6.73 × 10⁻³¹ and 3.59 × 10⁻²⁷.
The comparison claim is capped at 20 digits; long stored mantissas and working
precision are not treated as measured accuracy.

## Observed phase times

These are single observed runs, with other coordinated work on the host.
Original and native solvers use different precision, expansions and paths, so
this table does not establish a matched speedup or general performance parity.
The prototype links the immutable native library at source commit
`15471e4bbbc7736f1cce62bfb8ff1bad50821dd1`; its driver and executable hashes
are in the validation report.

| Solver/settings | Initialization | Infinity to −1 | −1 to 32 |
|---|---:|---:|---:|
| Native 60 digits/order 80 | 2.486 s | 0.335 s / 32 steps | 0.583 s / 42 steps |
| Native 80 digits/order 112 | 4.454 s | 0.477 s / 32 steps | 0.868 s / 42 steps |
| Original 500 digits/order 50 | 43.309 s | 20.805 s / 3 segments | 33.859 s / 7 segments |
| Original 600 digits/order 75 | 0.281 s (warm) | 36.295 s / 3 segments | 77.923 s / 7 segments |

Native initialization includes Frobenius construction and constraint solving;
original initialization is `PrepareBoundaryConditions`. Original refinement
runs in the same kernel, so its warmed initialization time cannot be compared
with the first initialization. Saved Padé representations at `t=20` are also
retained in the raw original output, with separate construction timings. Native
saved-series evaluation at that point is not part of this acceptance check.

The offline regression is `tests/banana_equal.rs`; it independently repeats
both native settings and compares every coefficient at both endpoints.
`tests/asymptotic.rs` covers logarithmic coupling, distinct integer resonance,
formal regulator sectors, fractional exponents, missing constants,
inconsistency, nonfinite inputs and unavailable series orders.

## Reproduction

Run the native regression without Mathematica:

```sh
cargo test --release --locked --test banana_equal
```

With the pinned checkout and a licensed kernel, set `DIFFEXP_ROOT` and a new
`DIFFEXP_OUTPUT` directory, then run the unchanged original package through the
independently written driver:

```sh
DIFFEXP_ROOT=/common/dev/diffexp \
DIFFEXP_OUTPUT=/common/dev/amflow/target/banana-oracle-repeat \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
timeout --signal=TERM --kill-after=10s 600s \
bern-wolfram -noinit -noprompt -script scripts/diffexp_banana_oracle.wl
```

`scripts/extract_diffexp_banana.wl` uses the same environment variables to
export the exact matrices and independently expand the analytic boundary,
without running DiffExp's numerical solver. It verifies all five matrix-file
hashes before loading them. The reusable extractor was run successfully and
its output exactly matched the inputs used by the numerical prototype.
