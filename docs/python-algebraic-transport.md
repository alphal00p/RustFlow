# General registered-root matrices from Python

`KinematicTransport` accepts exact general differential matrices with named
square-root generators through its optional `roots` argument. This extends the
ordinary dense constructor; canonical dlog construction remains available
through `KinematicTransport.canonical`. Both use the existing native connection,
continuation, cache and precision machinery in
`symbolica.community.hep.integration`.

For native expressions `s`, `epsilon`, `r`, `I1`, `I2`, `one` and `zero`:

```python
transport = KinematicTransport(
    epsilon,
    {s: [[epsilon / r, one / r], [zero, epsilon / r]]},
    [I1, I2],
    one,
    roots={r: s},
    branch_domain="chosen square-root sheet and admitted logarithm branch",
)
assert transport.roots == {r: s}
result = transport.evaluate(
    boundaries,
    {s: 4 * one},
    0,
    2,
    root_sheets={r: 1},
    admit_straight_path=True,
)
```

Here `r` denotes a registered root satisfying `r**2 = s`; it is not a substituted
numerical value. Source boundaries must belong to this transport identity and
carry their own `root_sheets`. Signs `+1` and `-1` are relative to the principal
root at each exact point. Missing or inconsistent signs produce typed errors.
`roots=None` and `roots={}` retain the rational constructor and its identity.
The `roots` property also exposes declarations on canonical transports.

Each physical derivative is supplied in the same basis. The native owner admits
epsilon-regular matrices, including epsilon-independent couplings and higher
epsilon orders. It retains original denominator domains and handles sums of
registered roots by exact quotient arithmetic. Conservative norm-domain checks
can reject a point regular on one sheet; Python does not relax that admission.
Radicands depend on the declared physical variables, not epsilon or other root
generators. Unsupported inputs fail through the existing typed exceptions.

For threshold crossings, bind a
[`ContinuationPrescription`](python-prescribed-transport.md) and explicitly
admit the prescribed homotopy. Ordinary and prescribed transport both support
binary cache restart, nearby-source selection, inherited uncertainty and exact
hits without new route admission. Finite singular limits use
[`EndpointRoute`](supplied-singular-endpoints.md), with matching root signs and
separate matching/final-approach admission. Terminal limits remain separate
from regular starting values.

All expressions and numerical values remain native Symbolica objects; the
binding introduces no root, ODE, graph or tensor algorithm. Computation releases
the GIL and uses the existing cancellation/progress control. Reported digits
remain limited by input evidence and independent refinement, not merely working
precision. General epsilon-singular basis transformations and constrained
registered-root endpoint sectors retain their native implementation limits.
