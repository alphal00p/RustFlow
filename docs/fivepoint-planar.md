# Full planar one-loop five-point transport

The benchmark transports all 13 master integrals through epsilon 4 from PH1 to PH6 in the ancillary data of [arXiv:2005.04195v2](https://arxiv.org/abs/2005.04195v2). The 58-letter canonical alphabet has 30 nonzero letters in this family and uses two square-root generators, `tr5` and `sqrtG3`. Rust constructs the six-variable differential system, differentiates the roots in each dlog, and pulls it back to the exact affine physical path.

The boundary is the supplied ancillary numerical vector. The fixture preserves each input's Wolfram precision and absolute accuracy, including finite-accuracy zero entries. Inexact inputs have a conservative absolute uncertainty floor of 1e-130. The low-level algebraic transport does not yet propagate that uncertainty into a certified output bound; independent precision/order refinement and comparison of all 65 coefficients provide the reported numerical validation. This is supplied-boundary differential transport, not automatic boundary generation or a full scattering amplitude.

The local prescriptions follow the pinned original DiffExp run: every listed physical invariant and the relevant square-root discriminants have `+i0`. Their simple zeros require both upper and lower detours in the affine parameter. The fixture records 80-digit root locations and derivative signs. The explicit 65-waypoint contour has a radius bounded by the nearest other pole. Additional apparent poles use upper detours; the complete endpoint comparison validates this particular route. This fixture-specific construction is not a general prescription planner.

The original DiffExp notebook settings are working precision 150, initial order 40, accuracy goal 15, division order 3, epsilon order 4, no Mobius transformation, no Pade, and one process. Native validation uses working precision/order pairs 40/48, 60/80 and 80/112 with 20 requested digits. All 65 coefficients in the 40/48 run agree with 80/112 to at least 30 observed absolute digits; the 60/80-to-80/112 difference is 2.29e-51 or smaller. These settings and algorithms do different work: original DiffExp crossed 97 real segments, while the native lower-precision contour accepted 575 steps and rejected 306. The timings are separate single-run phase observations, not a matched-settings or statistical speedup claim. Both were single-core runs on the same shared host (original CPU 28, native CPU 29).

| Solver / working digits / initial order | Matrix preparation | Compilation / configuration | Transport | Accepted / rejected |
| --- | ---: | ---: | ---: | ---: |
| Original DiffExp / 150 / 40 / goal15 |0.035 s (log matrix) |0.743 s |89.440 s |97 original segments |
| Rust / 40 / 48 / goal20 |20.019 s |0.071 s |59.553 s |575 /306 |
| Rust / 60 / 80 / goal20 |20.011 s |0.075 s |175.430 s |575 /306 |
| Rust / 80 / 112 / goal20 |20.027 s |0.073 s |343.004 s |574 /305 |

Original matrix preparation constructs the logarithmic canonical matrix, while native preparation constructs and pulls back all six physical derivative matrices. These preprocessing phases do different work; the reported transport times are the respective public transport calls.

The cheaper native run completed in 79.646 seconds including its exact physical-system construction; the original process took 92.030 seconds including configuration and boundary preparation. Native peak polled RSS was 26.5 MB for that run, compared with 1.09 GB for the original process tree. The original endpoint differs from the ancillary endpoint by 6.17e-23, whereas the cheaper native result differs by 1.88e-31. Both meet the 20-digit comparison used here. The original run requested only 15 digits and used a different segment algorithm, so these observations do not establish a general performance ratio.

The earlier exact-root-series native implementation was stopped at 600 seconds without an endpoint. A 2376-sample profile identified exact symbolic square-root series at numerical chart centers as the bottleneck. The optimized version reuses the library's rational Taylor solver for `r'=(R'/2R)r`, preserving branch and residual checks. No endpoint equality or timing ratio is claimed against the incomplete baseline.

Run the dedicated production benchmark explicitly with `cargo test --release --test fivepoint_planar -- --ignored --nocapture`. The explicit test takes about 11 minutes on the recorded single-core setup and runs all three native settings. It compares every coefficient against both the ancillary endpoint and the live original endpoint, checks precision/order refinement, and checks the endpoint signs of both square roots. No Mathematica runtime is needed for the Rust regression.

The original package is pinned to `784c8229bf92369a03f011a48e161522c8c54bbd`, and all three scientific input files are hash-pinned. The reusable independent Wolfram harness uses `DIFFEXP_ROOT`, `DIFFEXP_ANCILLARY` (the v2 `anc` directory), and `DIFFEXP_OUTPUT` (a fresh output directory). Raw original output, source metadata, complete native output, source/binary hashes, timing settings and the censored slow-chart baseline are retained in the accompanying reports.

The complete settings, outputs and hashes are in [the validation report](../reports/diffexp/fivepoint-planar-validation.json); the exact scientific data and finite-accuracy input metadata are in [the fixture](../fixtures/diffexp/fivepoint-planar-1loop.json). Regenerate the original with [the independent harness](../scripts/diffexp_fivepoint_oracle.wl) after setting its three environment variables.
