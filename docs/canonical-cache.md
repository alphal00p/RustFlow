# Canonical forms in the physical point cache

`RustFlow<CanonicalAlgebraicSystem>` uses the same mutable `RustFlowCache` as
rational and dense algebraic systems. Its constructor keeps exact logarithmic
letters and their constant matrices separate until a starting boundary and a
physical path have been selected. A path is substituted into the letters before
the differential matrix is assembled. The existing registered-root compiler,
shared adaptive series engine, independent refinement and weighted propagation
of supplied boundary uncertainty then evaluate that matrix.

```rust,ignore
let flow = RustFlow::with_canonical_conditions(
    canonical, &basis, &normalization, Prescription::PlusI0,
    "regular physical sheet", &reduction_conditions,
)?;
let result = flow.evaluate_to(
    &mut cache, &destination, Some(&root_germ), epsilon_range,
    &options, &context, &cost_policy,
)?;
```

The constructor reuses `CanonicalAlgebraicSystem::nonzero_conditions()`, which
retains the original letter and radicand domains before differentiation and
cancellation. Additional reduction conditions remain part of identity and are
checked before source ranking. Registered-root systems require a complete
explicit germ. A canonical form without roots is supported by passing `None`;
a germ is rejected in that case.

The cache identity records the canonical representation explicitly, including
its dimension, ordered physical variables, letters, constant matrices and root
registry, together with the basis, normalization, prescription, branch domain
and exact conditions. It is deliberately distinct from a dense system, even
when their differential equations are mathematically equal. Automatic
cross-representation boundary equivalence is not inferred. One bank may contain
both representations.

Binary schema 4 stores that representation and reconstructs/validates its exact
identity on load. Versions 1 through 3 and incompatible source/dependency
fingerprints are rejected. No numerical value is migrated or assigned stronger
accuracy during loading. Save/restart, compatible exact hits and accepted
intermediate insertion use the existing bank implementation.

The physical scope is unchanged: regular straight paths with real nonzero
radicands and explicit local root germs. The caller's cost policy and branch
domain must also establish any integral/logarithmic monodromy. Source errors
and independently compared transport changes remain numerical estimates, not
interval certificates. A large working precision does not improve the recorded
source evidence. Explicit complex contours and threshold crossings use the
separate native transport interface.

Validation covers progressive reuse and binary restart against an analytic
registered-root dlog solution, independent transport through its exact dense
connection at higher precision/order, distinct ordered identities, rootless
canonical transport, canceled-letter source holes and insufficient supplied
accuracy. These tests are distinct from the full five-point benchmark.
