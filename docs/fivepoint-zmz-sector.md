# Connected two-loop ZMZ sector

The original ZMZ masters 14 and 15 form a coupled two-loop sector with six propagators, including the mixed line `(l1-l2)^2`. Their exact differential-equation dependency closure contains 13 original masters: `[14,15,49,50,51,52,55,56,59,60,65,68,69]`. The fixture retains this original ordering and the scientific topology/numerator definitions from `pureBasis-zmz.m`; it does not relabel the integrals as one-loop masters.

The benchmark uses the pinned 2005.04195v2 PH1 and PH6 points and the supplied PH1 numerical boundary, through epsilon 4. The canonical subsystem has 152 alphabet coefficients, 53 matrix entries and the two registered roots `tr5` and `sqrtG3`. For every retained row, all nonzero columns of every exact alphabet matrix belong to the retained set. The original oracle independently asserts this closure before restricting its full 75-master tensor and boundary vector.

All 65 coefficients passed three comparisons: native 40 digits/order 48 against 60 digits/order 80 (maximum difference 3.97e-31), refined native against ancillary (3.37e-51), and native against the live original oracle (8.51e-23). The endpoint sheets also match the expected radicands. The comparison claim is 20 absolute digits. Source accuracy metadata and finite-accuracy zeros remain in the fixture with a conservative 1e-130 absolute allowance; the low-level solver does not return a certified propagated boundary-error bound.

| Run | Exact physical matrix preparation | Compilation | Transport | Accepted / rejected |
| --- | ---: | ---: | ---: | ---: |
| Native 40 digits/order 48 |5.341 s |0.067 s |66.157 s |575 /306 |
| Native 60 digits/order 80 |5.232 s |0.064 s |195.977 s |575 /306 |
| Original WP150/order40/goal15 | See raw phases | See raw phases |96.989 s |97 original segments |

These are single-run observations on a shared host, with native CPU 29 and original CPU 28. The original and native algorithms do different preprocessing and contour work, and the original requested 15 digits versus the native 20. The raw report preserves those differences rather than claiming a matched-settings speedup.

Run the ignored scientific regression with `cargo test --release --test fivepoint_zmz_sector -- --ignored --nocapture`. It takes roughly five minutes on the recorded single-core setup. The helper is shared with the full one-loop five-point benchmark. The reusable original harness requires `DIFFEXP_ROOT`, `DIFFEXP_ANCILLARY` and a fresh `DIFFEXP_OUTPUT`; no runtime is needed for the Rust regression.

This milestone is a genuine connected two-loop subsystem with supplied boundaries. It does not establish automatic boundary construction, an amplitude, or acceptance of the complete 75-master ZMZ benchmark. Full source hashes, phase timings and coefficient comparisons are retained in `reports/diffexp/fivepoint-zmz-sector-validation.json`.
