# Regular square-root transport

`algebraic::AlgebraicSystem` extends supplied differential systems with explicitly
registered generators `r` satisfying `r² = R(x)`. Radicands are exact rational
functions of the path parameter. Exact rational-complex constants are allowed;
undeclared variables, nested roots and floating input coefficients are rejected.
Specialize external parameters before compilation.

An ordinary system wraps one coefficient matrix with `AlgebraicSystem::ordinary`.
An epsilon system supplies `EpsilonSystem` matrices directly, using its existing
common Laurent leading power and triangular epsilon recurrence. The compiled
object is immutable and may be shared among concurrent transports. Every call
provides its own boundary, contour and initial root seeds.

Integer powers of every registered root are reduced with `r² = R`. This applies
to numerator and denominator before the denominator is checked. Thus `r⁻¹`
becomes `r/R`, and `1/(r²+1)` becomes `1/(R+1)`. After normalization the denominator
may have only one root monomial; sums such as `1/(r+1)` require algebraic
rationalization and are currently rejected. A denominator which vanishes under
the root relations is invalid.

`RootSeed::Principal` and `Opposite` choose the initial discrete sign.
`RootSeed::I0(Prescription::PlusI0)` chooses the upper side of a **real radicand**;
`MinusI0` chooses its lower side. This is a radicand prescription, not a global
choice of contour side in `x`. For example, the sign of `R'(x0)` determines how
an `x`-contour near a simple real zero approaches the radicand cut. Callers supply
the complex contour explicitly and are responsible for its physical homotopy.
`RootSeed::Value` uses a finite nonzero supplied value solely to select one of
the two signs; its magnitude is recomputed from `R` at working precision. Zero
or ambiguous hints are errors. This does not restore precision to a boundary.

At each accepted point, the next chart uses exact binary rational coordinates
for the rounded MPFR center and Symbolica's exact normalized expansion
`(R(x)/R(center))^(1/2)`. Native `Series` multiplication uses a fixed-precision
complex coefficient ring; all arithmetic remains MPFR. Root series predict the
analytic sheet, and an independently evaluated radicand determines its magnitude.
Both midpoint and endpoint root defects, root truncation, and squared-root
residuals must pass in addition to the shared solution checks. Rejected trials
leave branch state unchanged. The returned `BranchState` records the endpoint,
root identities in declaration order, and values.

Root zeros and poles join rational matrix poles in the shared step-radius
constraints. This entrypoint transports between regular points. Algebraic
Frobenius endpoints, arbitrary algebraic extensions, automatic physical path
planning, and algebraic `RustFlowCache` reuse are not implemented here. The
returned `EpsilonSolution` has no `verified_digits` or independently verified
checkpoints: compile and transport again at higher precision/order to establish
accuracy. Algebraic cache reuse additionally needs root-sheet identities and
root-aware uncertainty bounds; rational cache APIs must not be used to bypass
those requirements.

Regression tests cover one and two windings of `sqrt(x)` with
`y'=y/sqrt(x)`, giving `exp(-4)` and `1`, independent precision/order refinement,
direct epsilon coefficients, simultaneous transports on opposite sheets,
rejected-step nonmutation, sparse root-series defects, exact power normalization,
and invalid seed/registry inputs. The implementation is independently derived
from the chain rule and Taylor recurrences; it does not translate DiffExp source.
