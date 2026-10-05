# Exact coefficient collection loses large coefficients

This reproducer depends only on official Symbolica
`c3408e4ba1d3bdd4ea55678fad50e27009be13d4`, with no RustFlow or local Symbolica
patch. All inputs are exact integers and symbols.

For `(10^309 + 10^309*epsilon)*x`, `coefficient_list(&[x])` returns an empty
list. Reconstructing the expression from that list gives zero. The same example
with `10^308` returns the expected single coefficient.

```sh
cargo run --release --locked --manifest-path tools/mre/symbolica-coefficient-list/Cargo.toml
cargo run --release --locked --manifest-path tools/mre/symbolica-coefficient-list/Cargo.toml -- collect
cargo run --release --locked --manifest-path tools/mre/symbolica-coefficient-list/Cargo.toml -- exact
```

The first command reports every case. On the pinned revision, `collect` fails
its exact reconstruction assertion at power 309; `exact` passes all cases.
The passing control supplies `AtomField { statistical_zero_test: false, .. }`
to `try_to_polynomial`, then uses Symbolica's
`to_univariate_polynomial_list`. It does not implement another collection
algorithm.

`coefficient_list` calls `to_polynomial_in_vars`, whose internal fields use
statistical zero tests. A coefficient sum beyond the range of its internal
`f64` evaluator produces `Inconclusive`; `AtomField` treats that result as zero
and removes the coefficient. The public collection documentation has no
coefficient-size restriction. An exact collection API should preserve such
coefficients and reserve statistical identity testing for an explicit request.

This error was exposed by physical nonplanar Higgs–jet transport. Exact row
clearing and polynomial cross products passed, but coefficient extraction
discarded nonzero cleared numerators. The solver's independent differential
residual correctly rejected every step. The captured native sample and the
minimal reproducer are summarized in
[the validation report](../../../reports/validation/2026-10-05-symbolica-coefficient-list.json).
