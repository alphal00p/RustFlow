# Supplied-boundary endpoint transport

`RustFlow::evaluate_endpoint` evaluates a finite singular limit from an existing
verified regular boundary bank. It supports rational, registered-root and
canonical physical connections. `evaluate_prescribed_endpoint` first reaches
the regular matching point through the existing prescribed-contour engine and
requires its usual `HomotopyAdmission` as well as an `EndpointAdmission` for the
final approach.

An `EndpointRequest` contains the epsilon range, resource limits and an exact
`EndpointChart`: a `KinematicPath` with endpoint at parameter zero, a nonzero
matching parameter, matching root signs, logarithm winding and a declared
terminal homotopy. A chart such as `s = 1/z` describes infinity directly. The
same winding enters every fractional power and explicit logarithm. The caller
must admit the physical approach; a local convergence disk does not establish
global monodromy.

The chart is pulled back exactly, with its Jacobian. The existing registered-root
restriction of scalars and rational `PreparedFrobenius` produce the local basis.
For a rational connection the root list is empty. Native Symbolica matrices
invert the numerical matching basis once per profile and apply that inverse to
all physical unit directions together. Native complex balls check the inverse
residual and propagate arithmetic errors for that stored numerical matrix.
Exact polynomial coefficient extraction certifies a punctured disk free of
other source, root, normalization and chart poles. Only powers of the endpoint
parameter are removed from those guards. No graph, tensor, polynomial
convolution or independent Frobenius recurrence is introduced.

The result contains epsilon coefficient limits, their numerical evidence and
the regular matching boundary. Independent increased precision and order
profiles must agree. The endpoint transfer map propagates every recorded
matching error and an additional floor from its verified digit count. Its
achieved count cannot exceed the matching boundary's count or inherited input
evidence. These are numerical consistency and propagated-error estimates, not
a rigorous global forward-error certificate. A tiny central result does not
make a large absolute input uncertainty disappear.

Every independent physical input direction must have a finite limit. A central
solution obtained by tuning away a divergent sector is insufficient for this
independent-component error model. The endpoint recurrence preserves small
nonzero logarithmic rows and treats exact resonances without magnitude-based
zero classification. Exact support admission additionally refuses unresolved
nonpositive sectors rather than trusting rounded cancellation. Correlated or
constrained physical boundary sectors remain required future work for broader
parity; they are not silently discarded.

These are limits **after** a finite epsilon expansion. For example, the
logarithmic coefficient in `Y'=epsilon/s Y` diverges as `s` tends to zero and is
rejected. Selecting a symbolic dimensional sector such as `s^epsilon` is a
different operation requiring symbolic exponent provenance. No small finite
mass stands in for the endpoint.

`RustFlowCache::endpoint_entries()` contains terminal records separately from
`entries()`. Terminal limits never enter ordinary nearest-source selection:
finite values at a singular point do not specify initial data for its ODE.
`len()` counts regular records only; use `endpoint_entries().len()` for terminal
records (or `len(bank.endpoint_entries())` in Python).
Each endpoint retains a regular matching anchor, which is also inserted into
the ordinary bank for future continuation. Compatible terminal requests can
reuse a stored endpoint directly. The identity includes the exact physical
system, basis, normalization, guards and branch convention; terminal matching
also requires the complete chart, root signs, winding, homotopy and compatible
epsilon range.

All pending regular checkpoints and the endpoint are committed together only
after final validation and cancellation checks. Admission failure, divergent
input directions, insufficient evidence and cancellation leave the bank
unchanged. Schema 6 serializes both record types with native Atom payloads,
source/dependency fingerprints and integrity checks. Older schemas are rejected
explicitly; existing publication banks are not migrated or overwritten.

The API retains the rational solver's restrictions on irregular sectors,
indicial roots and Fuchsian normalization. Lift dimensions, series orders and
refinement attempts are bounded explicitly. Conservative support and disk
admission may reject valid problems requiring a more expressive sector or
convergence analysis. These limits, and large-system performance, remain open
scope rather than an assertion of full upstream parity.

The Python binding exposes `EndpointRoute`, `EndpointResult`,
`KinematicTransport.evaluate_endpoint` and `BoundaryCache.endpoint_entries`.
Construct the route from Symbolica expressions, including the exact matching
parameter. `admit_matching_path` admits the regular affine route;
`admit_endpoint` separately admits the final chart approach. Both default to
false, and terminal cache hits still require endpoint admission. The result
exposes its route, finite coefficients, comparison errors, digit evidence and
regular `matching_boundary`; ordinary transport diagnostics are attached to
that matching result. Cache merging preserves both terminal and regular records.
The prescribed-contour endpoint interface is currently exposed in Rust.
