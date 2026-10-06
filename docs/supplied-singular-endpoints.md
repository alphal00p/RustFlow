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
nonpositive sectors rather than trusting rounded cancellation. The separately requested exact-constraint interface below admits correlated
physical boundary sectors without inferring exact relations from numerical
values. Other constrained sectors remain open scope.

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
unchanged. Schema 7 serializes both record types with native Atom payloads,
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


## Explicit exact asymptotic constraints

The Rust methods `evaluate_constrained_endpoint` and
`evaluate_prescribed_constrained_endpoint` take an additional
`EndpointConstraints`. Each `ExactAsymptoticRelation` asserts an equality
between exact Gaussian-rational weighted coefficients selected by component,
exact rational power and logarithm power. A nonempty provenance string states
where this physical information comes from. Components in an epsilon hierarchy
use `epsilon_offset * physical_dimension + component`; the leading epsilon
power remains part of the request. These declarations are physical assumptions,
not facts recovered from a rounded boundary or a digit label.

For example, in `y1' = y2/z`, `y2' = 0`, declaring the coefficient of
`z^0 log(z)` in `y1` to be exactly zero removes the divergent integration
constant while leaving the finite amplitude free. A numerical `y2 = 0` alone
cannot establish this relation. An explicitly declared nonzero logarithmic
coefficient, including `10^-300`, fails the finite-limit proof.

Numerical and exact coefficient work use one shared Frobenius recurrence and
normalization. Symbolica owns the exact series, Gaussian-rational matrix row
reduction, inverse and products. Native exact row reduction returns an affine
space `c = c0 + N a` without fixing the remaining free amplitudes. Every
nonpositive power and logarithm in the complete known projected prefix is
checked on both `c0` and every column of `N`; padding from normalization does
not extend the known prefix. Inconsistent equations, unresolved coefficient
fields and a remaining divergent direction produce typed failures.

At each independently increased precision/order profile, the existing full
conditioned matching inverse recovers the free coordinates. The corresponding
constrained solution must reconstruct **every** supplied regular component
within the recorded input errors, their evidence floor and estimated numerical
errors. Both the finite transfer map and reconstruction propagate input
sensitivity and the native inverse residual bound. No numerical rank decision,
least-squares fit or tolerance-based exact zero enters the support proof.
These component checks may conservatively reject a feasible uncertainty box;
they are not a rigorous joint or global uncertainty certificate.

Terminal cache entries retain the exact relations, their provenance and the
original regular matching anchor. Constraints enter terminal identity;
unconstrained requests cannot reuse constrained endpoints. Constrained reuse
currently requires the exact epsilon range, because discarding higher-order
relations needs a separate projection proof. The limits controlling exact
work do not change the meaning of an already validated cache record. Schema 7
stores the declarations using the same native Atom codec and explicitly rejects
older schemas. Existing publication banks and runtimes are left untouched.

`ExactFrobeniusLimits` bounds dimension, order, coefficient height and scalar
storage/work estimates. Checks precede dense recurrence allocations, exact
series and native matrix operations, with cancellation checks between bounded
operations. These conservative admission estimates are not allocator or wall
clock quotas; symbolic preparation still has the existing Frobenius owner's
restrictions. The standalone `ExactFrobeniusBasis::constrain` API also exposes
the reusable affine space independently of endpoint transport.

This constrained interface currently admits rational connections, including
coupled finite epsilon hierarchies. A registered-root restriction of scalars
requires an additional exact proof that the constrained space stays on the
declared root sheet; the implementation refuses it until that proof is owned.
The ordinary registered-root endpoint interface remains available. Exact
non-Gaussian coefficient fields, symbolic dimensional-regulator sector
selection, larger systems and Python exposure of the new constraint declarations
remain open work. None of these restrictions redefines broader endpoint parity
as complete.
