# RustFlow physical transport and boundary cache

`RustFlow` is the native physical-transport object; `RustFlowCache` is a separate,
mutable bank shared across requests. Load a nonempty binary snapshot, pass a
mutable reference for each evaluation, and save the growing bank afterwards.
Neither object requires a Python runtime. The future Python API belongs in
HEPKit; this repository supplies the native core and the CLI.

The [physical-cache example](../examples/physical_cache.rs) seeds an exact
analytic boundary, evaluates one destination, evaluates a nearby destination
from the newly populated bank, and repeats the second request as a cache hit.

```rust,ignore
let mut cache = RustFlowCache::load(cache_directory)?;
// Insert an AMF, analytic, or supplied boundary with accuracy evidence if needed.
let result = flow.evaluate_to(
    &mut cache, &destination, epsilon_range, &options, &context, &cost_policy,
)?;
cache.save(cache_directory)?;
```

`KinematicSystem` holds exact partial-derivative matrices in a common basis.
`KinematicPath` performs simultaneous coordinate substitution with the full
chain-rule Jacobian. `EpsilonSystem` expands rational epsilon dependence exactly
and integrates its coefficient hierarchy directly. The recurrence clears
denominators once per physical row and shares those finite polynomials between
epsilon orders; it does not construct a dense augmented matrix.

`KinematicDerivative` derives invariant derivatives from an integral family,
including the change of the external Gram matrix. The older
`parameter_derivative` deliberately keeps scalar products fixed for auxiliary
and Feynman parameters. They describe different derivatives.

## Reusing computed physical points

Each cache record includes the full system/basis/normalization/branch identity,
physical coordinates, epsilon range, coefficient values, and accuracy evidence.
Exact input coordinates, numerical coordinates, and exact coordinate images of
numerically selected path parameters remain distinguishable. Retaining the exact
image of a path parameter prevents separately rounded invariants from leaving
the original kinematic line.

`RustFlow::evaluate_to` first selects a compatible starting boundary. An exact
coordinate hit returns the stored value; a nearby point supplies the starting
condition for a new transport to the actual requested destination. The default
provided distance policy uses arbitrary-precision scaled coordinate distances
and requires a caller-supplied path-admissibility check. Applications may provide
their own `TransportCost` policy, including singularities and expected step count.

Accepted intermediate points are compared at the same coordinates between
independent precision/order runs. Endpoint accuracy is not automatically assigned
to every intermediate segment. Boundary uncertainty is propagated with a
matrix-norm estimate, using lower bounds for rational denominators on each local
disk and subdivisions when those bounds are inconclusive. A cache record's
working precision cannot replace its verified accuracy. Stronger cached accuracy
requires explicitly passing tighter independent comparisons and input-error
checks, with a conservative arithmetic ceiling and two-digit margin.

These are consistency/error estimates, not interval-arithmetic proofs. Repeated
transport from the same uncertain input cannot improve the input's accuracy.
If propagated uncertainty exceeds the request, use a closer or more accurate
boundary. Input reference precision and inherited errors must be supplied
honestly, even when a reference contains hundreds of displayed digits.

Insertions are indexed by content keys and committed in validated batches;
evaluation does not clone the entire existing bank. Binary persistence retains
Symbolica expressions and precision-bearing numbers, checks schema and dependency
compatibility plus an authenticated payload, and replaces snapshots atomically.
Concurrent processes should coordinate snapshot writes or explicitly merge
their banks; independent writers can otherwise replace one another's snapshots.

## Current scope

The cache-first orchestration currently handles supplied epsilon-regular rational
systems and regular straight paths. The standalone coefficient solver supports
explicit complex polygonal contours and retains logarithm monodromy. Both use
the shared AMF Taylor engine and its tail/defect checks. Matrix epsilon poles
require a basis change; common Laurent leading powers of the solution are allowed.

Square-root branch-aware transport, automatic threshold prescriptions in the
physical orchestrator, general partial singular boundaries, and the larger
DiffExp examples remain in the [parity plan](parity-plan.md). The live
polylogarithm oracle is checked through independent boundary construction, not
by using its target value as a boundary.

The HEPKit integration will use its native DOT/graph conventions, Linnet for
graph manipulation, and Spenso/HEPKit tensor facilities. The ecosystem audit is
part of this work: native graph scalar-product ordering must be mapped explicitly
to the RustRed coordinate order, and graph symmetry factors must not be applied
twice.
