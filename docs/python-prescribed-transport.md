# Prescribed physical transport from Python

`ContinuationPrescription` binds exact physical threshold declarations to a
`KinematicTransport`. It lives in `symbolica.community.hep.integration` alongside
the existing transport and boundary-cache classes. Pass an ordered sequence of
`(Symbolica expression, "+i0" or "-i0")` pairs, a homotopy-domain description,
and an optional `unprescribed_side`. Exact expressions never pass through Python
floats. The native `PhysicalContinuation` owner checks allowed kinematic
variables, polynomial form, and original denominator domains when the declaration
is bound. Epsilon-dependent, root-dependent and nonpolynomial declarations are
rejected. Declaration order and the domain are preserved in the cache identity.

Both the ordinary constructor and `KinematicTransport.canonical` accept
`continuation=...`. For example, given native expressions `s`, `epsilon`, `Y`
and `one` (the exact constant 1):

```python
from symbolica.community.hep.integration import (
    ContinuationPrescription, KinematicTransport,
)

continuation = ContinuationPrescription(
    [(s, "+i0")], domain="upper s-plane approach", unprescribed_side="+i0"
)
transport = KinematicTransport(
    epsilon, {s: [[epsilon / s]]}, [Y], one,
    branch_domain="logarithm normalized at s=1", continuation=continuation,
)
# `boundaries` contains values and evidence inserted for this transport identity.
result = transport.evaluate(
    boundaries, {s: -one}, 0, 2, admit_prescribed_path=True
)
```

The existing native planner generates and checks the detours. Polynomial sides
do not establish global monodromy: `admit_prescribed_path=True` is the caller's
assertion that this declared homotopy applies to the integral. It is separate
from `admit_straight_path`, which cannot authorize a prescribed detour. Exact
compatible cache hits need no new route admission. Opposite prescriptions or
different homotopy domains have different identities and cannot share numerical
evidence merely because coordinates coincide. Registered roots retain explicit
source and destination sheets through the same native continuation owner.

For [singular endpoints](supplied-singular-endpoints.md), a bound continuation
also selects the native prescribed matching route in `evaluate_endpoint`.
`admit_matching_path` admits that route; `admit_endpoint` separately admits the
final local chart. A cached terminal still requires endpoint admission. Terminal
limits remain separate from regular ODE starting values, and binary save/load
and cache merging preserve both the terminal evidence and its regular anchor.

Computation releases the GIL and accepts the existing `ComputationControl`.
Typed native failures preserve unsupported inputs, missing evidence, inconsistent
sheets, exhausted limits and cancellation. The Python layer adds no contour,
Frobenius, graph or tensor algorithm, and no cache format. Local and accumulated
accuracy retain the native solver's documented scope; a declared contour is not
an accuracy certificate.
