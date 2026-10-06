# Exact endpoint relations in Python

`AsymptoticCoefficient`, `AsymptoticRelation` and `EndpointConstraints` live in
`symbolica.community.hep.integration`. They declare exact physical information
for the existing [singular endpoint solver](supplied-singular-endpoints.md).
Expressions remain native Symbolica objects. Construction and evaluation reuse
the Rust validation, exact recurrence, matrix and cache owners.

For the system `y1' = y2/s`, `y2' = 0`, assert the absence of the logarithmic mode:

```python
from symbolica import E
from symbolica.community.hep.integration import (
    AsymptoticCoefficient, AsymptoticRelation, EndpointConstraints,
)

log_mode = AsymptoticCoefficient(component=0, power=E("0"), log_power=1)
constraints = EndpointConstraints(
    [AsymptoticRelation([(log_mode, E("1"))], E("0"))],
    provenance="Regularity of the physical first component at s=0",
)
result = transport.evaluate_endpoint(
    cache, route, leading=0, last=0, constraints=constraints,
    admit_matching_path=True, admit_endpoint=True,
)
assert result.constraints.provenance == constraints.provenance
```

The transport, its numerical boundary and the `EndpointRoute` must already be
defined. The selected power is measured in that route's local parameter. A
relation asserts `sum(weight * coefficient) = value`; weights and values are
exact Gaussian-rational constants and powers are exact real rationals. A tiny
nonzero divergent coefficient remains nonzero regardless of requested accuracy.
Numerical boundary values alone never establish these relations.

For a common epsilon range, a selector's component is
`epsilon_offset * basis_size + component`, with zero offset at the requested
leading epsilon power. All supplied regular components must still match the
constrained solution within their propagated input/numerical error estimates.
Assertions do not turn working precision into measured accuracy.

The optional `constraints` argument works with ordinary and prescribed endpoint
matching. Omitting it preserves the existing unconstrained interface. This
milestone supports rational connections, finite epsilon hierarchies and
registered-root connections whose local leading constants are Gaussian rational.
The native solver derives exact sheet actions on the complete local solution
space and enforces their consistency together with the supplied relations.
Leading constants such as `sqrt(2)` remain unsupported. It does not perform
symbolic dimensional-regulator sector selection.

Relations and their provenance are retained in `EndpointResult.constraints`,
including results returned by `BoundaryCache.endpoint_entries()`. Binary
restart preserves them, and constrained terminal reuse requires the same exact
epsilon range and mathematical declarations. Different provenance is a distinct
assertion. Work limits (`max_dimension`, `max_order`, `max_coefficient_bits`,
`max_scalar_cells`) are validated but do not change a completed record's
mathematical identity. Source-sensitive schema-7 compatibility remains native
owned. Computation releases the GIL, and cancellation or failed admission does
not publish partial endpoint records.

These bindings are in this crate's optional `python`/`python_stubgen` features;
they are newer than the frozen community notebook extension.
