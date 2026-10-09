# Independent four-loop E7 reference

The native validation-only reference generator passed both tests, including its
explicitly invoked generation test. Its complete derivation is
[`docs/finite-density-e7-reference.md`](../../../../docs/finite-density-e7-reference.md).
It integrates neutral tensor bubbles and compact radial/angular beta functions,
with the independent charged mass derivative and nonzero upper surface retained.
Two agents independently checked its algebra, dimensions and normalization.

The reference file supplies both original targets in the seven-edge fixture:
the mandatory I37 definition and the supplemental raised occupied-line numerator
target. It binds the complete incidence, routings, masses, charges, polynomial
numerators and physical powers, and contains exact analytic Laurent coefficients
through order zero in the prescribed unexpanded reference normalization. Orders
below -2 vanish by the derived Gamma expression. The generator read neither
numerical oracle values nor native AMF predictions.

| Independent check | Largest observed relative difference |
| --- | --- |
| 50-to-80-digit native arithmetic at epsilon=1/5,1/10,1/20 | 1.131983000135411e-54 |
| Original beta/tensor versus duplicated Gamma expression at 80 digits | 2.291164926332169e-84 |
| Finite-coefficient extrapolants at epsilon=1e-20,1e-30 versus exact expression | 5.095944890816612e-18 |

All six finite-dimension precision comparisons, six independent-expression
comparisons and four Laurent extrapolants passed their recorded assertions.
The exact coefficient expressions have an analytic provenance; the extrapolant
changes are empirical checks, not rigorous interval error bounds. Actual native
AMF prediction stability and acceptance remain separate, unperformed gates.

Prebuilt execution took 0.03672115 seconds with peak resident memory 12332 KiB;
compilation and Nix startup are excluded. `reference.json` retains all native
decimal strings, exact formulas, source and lockfile digests. `resources.json`
and `resources.log` preserve the invocation and complete actual test output.

```sh
RUSTFLOW_DENSITY_E7_REFERENCE_REPORT=reports/validation/2026-10-09-finite-density-native-assembly/independent-e7-reference/reference.json \
  nix develop --command cargo test --locked --release --test finite_density_e7_reference \
  -- --include-ignored --nocapture --test-threads=1
```

This supplies one previously missing raised-line reference. The other two
four-loop supplemental references and all three families' native AMF evaluations
remain required. Production code does not import this test or its results.

After the independent reference was saved, a separate validation-only comparison
with supplied I37 passed all three Laurent coefficient identities exactly in
Q[pi], resolving the supplied `a1=243/4-3*pi^2/4`. The saved independent native
100-digit coefficient values agreed with a Decimal Machin-pi evaluation to a
maximum relative difference of 2.522e-100. The full evidence is in
[`supplied-I37-reference-comparison.json`](supplied-I37-reference-comparison.json),
including immutable source hashes and the original expressions. The original
reference and oracle bytes were preserved. This is **reference versus reference**:
one supplied record and three coefficients compared, zero AMF predictions read,
zero native AMF-to-oracle comparisons or acceptance claims.

```sh
python3 tools/finite_density/compare_e7_references.py
```
