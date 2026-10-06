# Reusable exact algebraic preparation

`AlgebraicSystem::prepare()` creates a `PreparedAlgebraicSystem` containing an
immutable snapshot of one exact, pulled-back differential system. It validates
and normalizes registered-root coefficients with the existing Symbolica-owned
quotient algebra, retains the root equations and source domain conditions, and
records the sparse coupling/epsilon indices needed for compilation.

Preparation also retains native exact rational row numerators/denominators,
denominator factors with their multiplicities and constant content, and cleared
common-denominator rows. These are produced by the same Symbolica conversion,
GCD, division, multiplication and factorization operations as direct compilation.
Factor scanning retains its original first-occurrence order.

`prepared.compile(precision)` creates a fresh numerical system at that precision.
It still performs coefficient extraction, numerical coefficient specialization,
exact stored-dyadic source specialization, source residual setup and pole solving.
No rounded coefficients, pole approximations, root values,
boundary values or accuracy evidence are reused across precision profiles.

For example:

```rust,ignore
let prepared = system.prepare()?;
let first = prepared.compile(Precision::decimal(60)?)?;
let refined = prepared.compile(Precision::decimal(100)?)?;
```

The source is borrowed read-only through `prepared.source()`. Later modifications
to the caller's original `AlgebraicSystem` do not affect a preparation; prepare
a new snapshot to incorporate them. A preparation retains its exact source,
normalized entries, cleared rows, factors and index tables for its lifetime.
It is neither a global cache nor a serialized boundary record. Native algebra retains its existing
resource and interruption limits.

`prepare_with_context` and `compile_with_context` check cancellation around
native operations. A single native algebra/compiler call is not preempted.
An interrupted operation returns the typed cancellation error and does not
publish a partly compiled object or insert a physical cache entry.

Regular and prescribed physical routes prepare once after pulling back the
source, then share that exact object across route admission, precision/order
refinement and boundary uncertainty propagation. Every numerical compilation
remains fresh. Exact refinement, root-sheet comparison, conditioning, propagated
input uncertainty and source accuracy caps retain their existing roles.

Singular-endpoint setup keeps its raw pullback and existing Frobenius preparation
path. Coordinate transformations that change an exact local system still compile
that transformed source separately. Rational-only routes use the compatibility
compiler, which prepares exact rows afresh on each call; they do not yet retain
prepared rows across precision profiles. Numerical parameter values are passed
fresh to every row compilation, and exact source checks use their full stored
precision independently of the requested numerical working precision.

This does not change the physical boundary cache format, compatibility contract
or admission requirements. Its normal implementation source fingerprint changes
with the library source, as for other implementation updates. Preparation alone
is not an accuracy certificate or a dimensional-sector endpoint selector.
