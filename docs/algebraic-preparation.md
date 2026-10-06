# Reusable exact algebraic preparation

`AlgebraicSystem::prepare()` creates a `PreparedAlgebraicSystem` containing an
immutable snapshot of one exact, pulled-back differential system. It validates
and normalizes registered-root coefficients with the existing Symbolica-owned
quotient algebra, retains the root equations and source domain conditions, and
records the sparse coupling/epsilon indices needed for compilation.

`prepared.compile(precision)` creates a fresh numerical system at that precision.
It still performs native rational row compilation, exact denominator clearing
and factorization, numerical coefficient specialization, source residual setup,
and pole solving. No rounded coefficients, pole approximations, root values,
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
normalized entries and index tables for its lifetime. It is neither a global
cache nor a serialized boundary record. Native algebra retains its existing
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
that transformed source separately. Rational-only routes are unchanged apart
from cancellation checks around their existing compilation calls.

This does not change the physical boundary cache format, compatibility contract
or admission requirements. Its normal implementation source fingerprint changes
with the library source, as for other implementation updates. Preparation alone
is not an accuracy certificate or a dimensional-sector endpoint selector.
