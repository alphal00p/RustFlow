# Unequal-mass banana: all 15 masters

The complete 15-master three-loop banana system from DiffExp is evaluated
through ε⁴ at `p²=50`, with squared masses `(2, 3/2, 4/3, 1)`, in `D=2−2ε`.
All 75 coefficients agree with a live calculation using the unchanged original
DiffExp package to **20 absolute decimal digits**. The native calculation starts
from analytic gamma factors at equal-mass infinity; numerical reference values
are used only for comparison. This is a supplied-DE and boundary acceptance
test, without new IBP reduction or amplitude evaluation.

The [fixture](../fixtures/diffexp/banana-unequal.json) contains the exact rational
connections and all endpoint references. The [original report](../reports/diffexp/banana-unequal-original.json)
and [native report](../reports/diffexp/banana-unequal-validation.json) retain
source hashes, complete values, precision metadata, commands and timings.
Scientific input comes from DiffExp commit
`784c8229bf92369a03f011a48e161522c8c54bbd`. Its GPL implementation remains in
the separate, unchanged checkout; no solver source is copied into this library.

## Analytic boundary and independent paths

The [equal-mass analytic boundary](banana-equal.md) determines four masters.
At equal masses they populate the 15-master basis through the exact map
`[0,0,0,0,0,0,1,1,1,1,2,3,3,3,3]`. The first six entries are double-raised
integrals, the next four single-raised integrals, then the scalar banana and
four factorized tadpole products. No extra numerical boundary is required.

The full native route continues from infinity to equal-mass `t=1/2`, deforms
the squared masses along `(1+x, 1+x/2, 1+x/3, 1)` for `x=0…1`, then continues
`p²=1/2…50`. A second route continues the equal-mass system to `t=50` before
the same mass deformation at fixed `p²=50`. Upper-half-plane continuation
implements the `+i0` prescription. The coordinate `mm_i` in the upstream
matrices denotes a squared mass.

The native implementation uses the existing generalized Frobenius expansion
and partial asymptotic constraints at infinity, followed by `EpsilonSystem`
transport. Only the infinity initialization uses the exact augmented
20-component hierarchy; regular transport retains the four- or 15-master
epsilon structure. The fixture contains the exact chain-rule pullback for
each of the three lines. Only the ε⁰ and ε¹ connection matrices are nonzero.

| Check, over every applicable coefficient | Maximum absolute difference |
|---|---:|
| Full route: 60 digits/order 80 → 80 digits/order 112 | 5.97 × 10⁻⁵⁸ |
| Full route versus the independent native route | 3.79 × 10⁻⁷³ |
| Native full route versus original shortcut, all 75 coefficients | 1.42 × 10⁻³⁴ |
| B11 versus the paper's printed 55-place table | 8.13 × 10⁻⁵⁶ |
| Four analytic tadpole products, 20 coefficients | 1.54 × 10⁻⁸³ |

The last two figures were measured with the earlier augmented native solver;
the production epsilon solver differs from its refined endpoint by at most
10⁻⁸¹. The offline production regression independently checks these identities,
both routes, both precision/order settings, and every original endpoint
coefficient. The original comparison claim remains capped at 20 digits;
85-digit stored references, 1000-digit working precision and reported error
estimates are not treated as measured accuracy.

## Original runs and observed times

The original full notebook route was stopped by its configured 16-GiB
process-tree RSS limit after **1599.51 s**. It had completed 10 of 26 physical
segments; its last printed point was approximately `p²=5.61546`. That display
is not a saved reusable boundary. The run requested saved physical series,
as in the notebook, and produced no final `p²=50` result or completed runtime.
Its completed mass-deformation stage at `p²=1/2` agrees with all 75 native
coefficients to 1.52 × 10⁻³⁴. The memory-censored result is preserved separately.

A distinct original run followed the equal-`t=50` route and completed in
**350.10 s**, with a peak polled process-tree RSS of **3.68 GB**. This provides
the complete original endpoint oracle. Its settings were:

- Equal-mass initialization/transport: 500-digit working precision, chop
  precision 250, order 80, division order 2, radius parameter 1.
- Unequal-mass transport: the notebook's 1000-digit working precision, chop
  precision 500, order 70, division order 4, radius parameter 10.
- Möbius maps and Padé enabled; one kernel and one CPU. The full-route seed
  instead used order 70 and division order 4, as recorded in its report.

The original shortcut spent 42.358 s preparing the analytic boundary,
35.880 s reaching `t=−1`, 65.492 s reaching `t=50`, and 193.315 s on the mass
deformation. Its reported error maximum at the final point was approximately
2.10 × 10⁻³². The original notebook duplicates values into the unequal basis
without propagating the equal-stage error array; both arrays are retained in
the report. Independent comparisons therefore determine the accuracy claim.

The production native standalone driver, linked against the immutable release
library at `c74011e`, took **16.78 s** for the full 60-digit/order-80 route,
**25.63 s** for the full 80-digit/order-112 route, and **8.04 s** for the
80-digit/order-112 alternate route. Each includes analytic initialization.
These are single observations on a shared host: original CPU 28, native CPU 29.
The solvers use different precision, orders, paths and saved-series work.
These measurements do not establish a matched speedup or general performance
parity.

The report also preserves the paper/notebook distinctions: notebook chop
precision 500 versus the paper's baseline 250, and a conflicting extra
`(1+3ε)` factor in one plot caption. The test follows the exact supplied
matrices and notebook basis. Historical paper timings are not current-host
measurements.

## Reproduction

The native test requires no Mathematica:

```sh
cargo test --release --locked --test banana_unequal
```

To reproduce the original shortcut with an already licensed kernel and the
pinned checkout, choose an available CPU and a new output directory:

```sh
python3 scripts/run_diffexp_banana_unequal.py \
  --upstream /common/dev/diffexp \
  --kernel /home/ben/.local/bin/bern-wolfram \
  --output target/banana-unequal-oracle-repeat \
  --route shortcut --cpu 28 --timeout 1800 --rss-bytes 17179869184
```

`--route full` reproduces the full notebook route, including saved physical
series and the planned Padé evaluation at `p²=10`; that original run did not
finish under the recorded memory cap. The runner polls RSS every two seconds,
so its measured peak may slightly exceed the requested cap. It records the
exit reason and terminates the process group on resource exhaustion.

The reusable driver checks the package and all equal/unequal matrix hashes,
accepts only a fresh output directory, and preserves original precision and
error metadata. The numerical reports retain the hashes of the actual target
drivers used for the measured runs; the reusable driver generalizes their
paths and configuration without an additional full numerical rerun.

`scripts/extract_diffexp_banana_unequal.wl` uses `DIFFEXP_ROOT` and a new
`DIFFEXP_OUTPUT` directory to regenerate the three rational line connections
without running a numerical solver. Its exported matrices, coordinate maps,
endpoints and pole polynomials were checked exactly against the accepted fixture.
