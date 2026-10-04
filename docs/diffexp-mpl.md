# Remaining DiffExp multiple-polylogarithm profiles

The pinned DiffExp notebook's trailing-zero, complex-letter and weight-20 examples
pass independent native transport and original Wolfram comparisons at 20 decimal
digits. Together with `mpl-101-4.json`, these cover every numerical MPL profile
in `MultiplePolylogarithms.nb` at commit
`784c8229bf92369a03f011a48e161522c8c54bbd`.

The new input sequences are `G(1,-10,0;4)`,
`G(10,-10+i,-1/2,-50;1)`, and `G(1,2,...,20;21)`.
Native code builds the triangular differential equation directly. The first
case obtains its logarithmic boundary from exact regularized constants and
native generalized Frobenius recurrence; the other cases start from exact
finite constants at zero. No original numerical value enters a native boundary.
The frozen native ordinary solver accepts the complex letter using an exact
symbolic parameter mapped to MPFR `i`. A literal complex polynomial coefficient
was rejected by that frozen version; its failed probe remains recorded. The
production compiler now encodes exact complex literals before rational-polynomial
conversion. Its integrated MPL regression supplies the literal complex letter
directly; the historical timings below remain those of the frozen symbolic-input
implementation.

The original harness executes the GPL notebook's own helper definitions from
the pinned checkout, including shuffle/logarithm extraction. This is oracle
execution, not source copied into the MIT library. It uses the notebook's
250-digit working precision, chop 225, Mobius transformation, division 3, and
orders 50/75 or 100, with parallel execution disabled. Native settings and
continuation contours differ, and the trailing-zero computation solves a
different number of systems, so the following single-core measurements do not
establish a matched-algorithm speedup. They were collected on CPU 28 of a shared
host while other scientific jobs ran on different cores.

| Case | Original order | Original helper seconds | Native digits/order | Native boundary+compile+transport seconds | Native/original absolute difference |
|---|---:|---:|---:|---:|---:|
| trailing | 50 | 0.768252 | 80/50 | 0.512075 | 3.998e-25 |
| trailing | 75 | 0.926068 | 100/75 | 0.355881 | 6.140e-38 |
| complex | 50 | 0.443258 | 80/50 | 0.085049 | 8.781e-28 |
| complex | 75 | 0.510206 | 100/75 | 0.060672 | 2.368e-40 |
| weight20 | 100 | 33.669108 | 100/100 | 5.974607 | 5.504e-57 |

The weight-20 independent refinement used 140 digits/order 140 and took
10.181645 seconds internally. All native basis components were checked against
increased working precision and order; maximum differences were 2.96e-76,
1.13e-79 and 8.27e-96 respectively. These small differences are evidence of
stability, not a claim that working precision certifies the same accuracy.
Both absolute and nonzero-MPL relative 20-digit checks pass.

Every original component retains its recorded `Precision` and absolute
`Accuracy`; comparison uncertainty estimates round absolute accuracy downward.
The weight-20 imaginary output is a finite-accuracy zero (absolute accuracy
56.5038), and is checked with an absolute bound only. The original shuffle
helper's signed `pm` error coefficient may cancel and is retained as an estimate,
not a rigorous bound. Two early harness-only decoding failures introduced a
machine zero via `pm -> 0`; the accepted runs extract the exact coefficient of
`pm^0` instead. They do not alter DiffExp's algorithm.

`reports/diffexp/mpl-additional-validation.json` retains all output arrays,
component metadata, source/binary hashes, internal phase times and step counts.
Two-second process polling can overestimate short-run wall time and miss the
RSS peak, so internal phase timings are reported independently. The offline
regressions need no Wolfram installation:

```sh
source .dev-env
cargo test --locked --release -p symbolica-amflow --test mpl_additional
```

To regenerate one original profile, use a fresh existing output directory and
the pinned DiffExp checkout. `MPL_CASE` is `trailing`, `complex`, or `weight20`;
`MPL_ORDER` is 50, 75, or 100 respectively. The harness checks package/notebook
hashes and preserves native Wolfram output alongside JSON:

```sh
mkdir -p target/mpl-original-new
DIFFEXP_ROOT=/common/dev/diffexp MPL_OUTPUT="$PWD/target/mpl-original-new" \
MPL_CASE=trailing MPL_ORDER=50 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5s 600s \
  bern-wolfram -noinit -noprompt -script scripts/diffexp_mpl_additional_oracle.wl
```

The larger weight-20 oracle used a 1800-second cap and 4 GiB process-tree RSS
limit in the recorded run. It completed in 33.669108 seconds in the original
helper (36.012178 seconds including launch/polling).
