# Canonical physical differential systems

`algebraic::CanonicalAlgebraicSystem` stores an exact connection
`dY = epsilon * sum_a M_a dlog(L_a) Y` as its ordered letters, constant matrices,
physical variables and registered square roots. The representation is immutable
and shared through `Arc`. It postpones physical matrix assembly until an exact
path is supplied.

```rust,ignore
let connection = CanonicalAlgebraicSystem::new(
    epsilon, &variables, &letters, &constant_matrices, roots,
)?;
let system = connection.pullback(&path, last_epsilon_order)?;
let compiled = system.compile(precision)?;
```

The pullback substitutes the path in each letter and radicand before taking
its total logarithmic derivative. Native Symbolica differentiation includes
root derivatives. This avoids constructing and cancelling a separate dense
multivariate matrix for each physical variable. The resulting `AlgebraicSystem`
uses the existing Taylor engine, root charts, branch continuation and supplied
boundary interface. Dense `AlgebraicKinematicSystem` inputs remain supported.

Source restrictions are computed from the original expressions before
cancellation. Each letter must be nonzero and defined, each radicand must be
nonzero and defined, and each rational coordinate map must be defined. Original
negative-power bases survive even when a stationary Jacobian, a matrix sum, or
native rational cancellation removes their poles from the final connection.
These restrictions remain available through `nonzero_conditions()` and are
carried into the compiled system. Formal quotient norms conservatively require
invertibility on all registered sheets; they may exclude a removable pole on
one selected sheet.

Root-power grouping uses native exact rational polynomials and native
`to_polynomial_in`. It does not use statistical general-expression zero tests.
Signed powers and exact complex coefficients retain their algebraic values;
numerical evaluation still uses Symbolica `Complex<Float>` at explicit precision.

The cache adapter is `RustFlow::new_canonical` or
`RustFlow::with_canonical_conditions`. It uses the same progressive
`RustFlowCache`, source selection, independent refinement, inherited-error
propagation, and transactional insertion as dense connections. Its `evaluate_to`
takes `Some(&root_germ)` when roots are registered, or `None` for a canonical
system without roots. The regular real-path restrictions in
[algebraic caching](algebraic-cache.md) still apply.

A canonical cache identity retains the ordered variables, letters and constant
matrices, root registry, source guards, basis and normalization. Its representation
is distinct from a dense connection identity: algebraically equivalent inputs
are not silently assumed to define interchangeable cached boundaries. Binary
schema 5 stores this representation and optional prescribed continuation directly,
rejecting older snapshots.

Regressions compare canonical and dense pullbacks exactly, preserve removable
source and coordinate holes, and check a complex root letter against the
analytic epsilon series of a logarithm at increased precision and order. The
cache tests cover restart, progressive source reuse, roots-empty data, identity
separation and inherited accuracy. Large-family preparation and transport timings
are separate experiments; these tests alone establish no performance parity.
