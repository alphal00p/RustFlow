The connected equal-mass two-loop ordinary boundary coefficient passed all three native profiles against an independent Schwinger-parameter reference at effective dimension `12/5`.

The native family is `K²−1`, `L²−1`, `(K−L)²−1`, with powers `(1,1,1)` and measure `d^D K/(i π^(D/2)) d^D L/(i π^(D/2))`. The harness uses the existing recursive boundary owner with `TadpolesOnly` and `bubble_subloops=false`. Each profile performed three native preparations. No reference value or special formula was supplied to that owner.

| Digits / series order | Relative real discrepancy | Absolute imaginary remainder |
|---|---:|---:|
| 18 / 60 | 1.18e−54 | 6.54e−53 |
| 28 / 60 | 3.13e−56 | 2.16e−55 |
| 28 / 80 | 2.02e−64 | 4.39e−64 |

All profiles used 40 guard digits and fresh recursive contexts. The adjacent native precision/order changes were below `1.31e−53` and `4.32e−56`. The entire native run completed in 2.17 seconds with peak child RSS 18 MiB, within the predeclared 300-second bound.

The independent calculation used standard-library decimal Gauss–Legendre quadrature on a sector-resolved simplex, plus a separately bounded positive-real Gamma evaluation. It refined precision and quadrature order and repeated the calculation with a different endpoint map. Its total measured runtime was 4.70 seconds. The native value is negative because three Minkowski denominators contribute three Wick signs. The empirical quadrature checks support the predeclared 10-digit comparison; the much smaller observed discrepancies are not a rigorous 54-digit error certificate.

The native predictions were saved before the comparison process opened the reference. `comparison.json` binds the complete predictions, reference, scripts, and frozen coherent v2 libraries. `artifact-manifest.json` records the final report files. No Cargo invocation, production change, supplied oracle read, or tolerance adjustment was made for this validation.

This checks the connected ordinary hard coefficient that can occur in the three-loop occupied boundary. It does not establish full three-loop finite-density closure, boundary assembly, or final target accuracy; those require their own complete-cut comparisons.
