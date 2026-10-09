# Independent massive Laurent reference

Status: independent reference generation passed. No full finite-density AMF
prediction or supplied-oracle answer was used in this calculation.

The exact input is recorded in `independent-laurent-reference.json`: the massive
two-loop sunset with squared masses `(1/4,1/4,1)`, chemical potential 1, the
scalar target and the original raised numerator `g1_2+u1*u2`. Values use the
unscaled Euclidean measure. Every vacuum, single-cut, double-cut and explicit
Fermi-surface contribution is retained, together with both complete sums.
The reference contains Laurent powers −2 through 0.

The derivation is in [finite-density-reference.md](../../../../docs/finite-density-reference.md).
Six Schwinger sectors have their local ultraviolet boundary subtracted
analytically. The compact and Feynman-parameter integrals use independent
multiprecision Gauss–Legendre quadrature. Ten exact epsilon samples per run
are interpolated with exact rational Lagrange weights after multiplying by
epsilon². No production reduction, AMF boundary, transport or fitting owner
is called. Samples are saved before interpolation and comparison.

The four runs independently vary nodes (48→64), epsilon spacing
(1/10000→1/20000), and working precision (60→80 decimal digits). All 108
coefficient comparisons pass relative 1e−12 for resolved nonzero values and
absolute 1e−20 for zeros in the raw Euclidean normalization. The largest
observed relative changes are:

| Refinement | Largest relative change |
|---|---:|
| Quadrature nodes | 3.481e−17 |
| Epsilon grid | 1.005e−22 |
| Working precision | 8.982e−53 |

Six independently derived analytic ultraviolet residues also pass after
interpolation; they are checks, not fit constraints. Three preliminary tests
verify quadrature moments, rational Laurent interpolation, and agreement of
the new subtraction with the previous convergent D=12/5 reference.

The reference process took 86.70 seconds and peaked at 12,372 KiB child RSS.
Compilation and Nix startup are excluded; see `resources.json`. Refinement
changes provide empirical error evidence, not rigorous interval bounds.

Reproduce the generation from the repository root:

```sh
cargo test --release --test finite_density_reference --no-run
RUSTFLOW_DENSITY_LAURENT_REFERENCE_REPORT=reports/validation/2026-10-09-finite-density-native-assembly/independent-laurent-reference \
  cargo test --release --test finite_density_reference \
  generate_independent_massive_sunset_laurent_reference -- --ignored --exact --nocapture
```

The main feature still needs a source-closed double-cut system and complete
AMF predictions before these coefficients can establish its numerical acceptance.
