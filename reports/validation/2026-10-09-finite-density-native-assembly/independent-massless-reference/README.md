# Independent massless two-loop sunset reference

The independent reference generator passes, with four ordinary mathematical
checks and one explicit generation check. It reads only the exact massless
sunset definition and uses native arbitrary-precision Gamma/Beta arithmetic
and direct convergent quadrature. It imports no finite-density evaluation owner,
reads no supplied oracle and reads no RustFlow prediction. Native numerical
acceptance is a separate comparison of predictions saved beforehand.

The reference covers the complete raw Euclidean scalar and raised original
`g1_2+u1*u2` targets at D=7, 13/2 and 15/4, with chemical potentials 1 and 3/2,
and Laurent orders [-2,0] at mu=1. The double-cut mass derivative includes all
three bulk terms and the nonzero moving upper surface. Vacuum and single-cut
zeros use the stated common regulated thermal dimensional prescription; they
are not separately assigned propagator prescriptions. The full derivation is
in [finite-density-massless-reference.md](../../../../docs/finite-density-massless-reference.md).

Direct D7 quadrature evaluates the original differentiated kernels before their
radial/angular analytic integration. The independent angle map removes the
square-root endpoints. The configurations `(digits, radial nodes, angular nodes)`
are `(50,6,32)`, `(80,6,32)` and `(80,8,48)`. Each of the six recorded scalar,
bulk, surface and total quantities passes a 17-digit relative check against
its analytic Beta value. The largest initial relative discrepancy is
1.4017296143e-19; the finest discrepancy is 3.9948091328e-31. Changing only
precision shifts the values by at most 1.5889428067e-54. These are empirical
comparisons, not rigorous quadrature error intervals. The Beta references also
pass separate 60/90-digit refinement checks with a 50-digit relative criterion.

At D7 the raised bulk sum cancels exactly while the upper surface remains
nonzero, making this a direct check that the support derivative is retained.
The tests additionally check D13/2, the meromorphic continuation to D15/4,
homogeneity under mu=3/2, the rational-pi special values, and Laurent pole/finite
coefficients by independent small-epsilon limits.

The four ordinary checks took 0.091402 seconds, peak child RSS 15412 KiB. The
explicit reference generation took 0.073674 seconds, peak child RSS 15380 KiB,
excluding compilation and Nix startup. Binary, reference-source and fixture
SHA256 digests are in `ordinary-resources.json` and `generation-resources.json`.
An earlier generation resource record is preserved: its wrapper read the large
linked binary into memory before forking, contaminating the inherited RSS
measurement. The corrected resource run uses a fresh streaming-hash wrapper;
both generation runs passed, and no reference values were changed to match
native results.

To reproduce generation after building the test target:

```sh
RUSTFLOW_DENSITY_MASSLESS_REFERENCE_REPORT=/tmp/independent-massless-reference.json \
  nix develop --command cargo test --locked --release \
  --test finite_density_massless_reference \
  generate_independent_massless_sunset_reference -- --ignored --exact --nocapture
```

The separate `compare_massless_reference.py` requires four completed sample or
five completed Laurent native profiles. It checks all sample cut sectors and
the total, precision/order/start refinements, and the separate Laurent epsilon
grid. It uses relative 1e-12 above magnitude 1e-20 and absolute 1e-25 below that
threshold, including imaginary-zero checks. Subsequent saved native comparisons now pass at D=13/2, D=15/4 and through
Laurent order zero; they are recorded separately in
[the native numerical aggregate](../massless-native-validation.json).
The independent reference generation itself reads no native predictions.
