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

`AlgebraicKinematicSystem` stores exact derivative matrices in several physical
invariants together with the same root registry. Its `canonical_dlog` constructor
uses the exact total derivative, including `dr/dv = (dR/dv)/(2r)`, for each letter.
The constant letter matrices remain exact. `pullback(&KinematicPath, order)`
substitutes all physical coordinates simultaneously in the matrices and radicands,
applies the path Jacobian, and returns an `AlgebraicSystem` through the requested
epsilon order. The path must specify every physical variable, with a parameter
distinct from those variables, epsilon, and all root symbols. Source entries
and root radicands are checked before multiplying by the path Jacobian. A
stationary coordinate therefore cannot hide a source pole or a path that lies
identically on a root branch locus; such paths require a separate limiting
prescription and are rejected.

Integer powers of every registered root are reduced with `r² = R`. This applies
to numerator and denominator before the denominator is checked. Thus `r⁻¹`
becomes `r/R`, and `1/(r²+1)` becomes `1/(R+1)`. Root sums in denominators are
inverted with native Symbolica `AlgebraicQuotient` arithmetic over rational
functions, eliminating one generator at a time while preserving its name and
sheet. For example, `1/(r+1)` becomes `(r-1)/(R-1)`. Fallible inversion rejects
nonunits, including dependent-root cases requiring a sheet-specific domain.
A denominator which vanishes under the root relations is invalid.

`AlgebraicKinematicSystem::nonzero_conditions()` exposes the exact source domain,
including root norms and radicand numerators/denominators. Pullback retains these
in `AlgebraicSystem.nonzero_conditions` before Jacobian or matrix cancellation;
epsilon-dependent conditions retain their first nonzero epsilon coefficient to
preserve the regular expansion's valuation. Compilation keeps the resulting
poles in its exact singularity metadata and checks them at starting points and
subsequent charts. A stationary coordinate and a canceled radicand denominator
therefore cannot erase an isolated source hole.

These formal norm conditions are conservative: `1/(r+1)` with `r²=x` excludes
`x=1`, even on the positive sheet where that particular expression is regular.
Continuation through such sheet-specific removable norm poles is not supported
by this interface. Supplied rational nonzero conditions receive the same checks.

`RootSeed::Principal` and `Opposite` choose the initial discrete sign.
`RootSeed::I0(Prescription::PlusI0)` chooses the upper side of a **real radicand**;
`MinusI0` chooses its lower side. This is a radicand prescription, not a global
choice of contour side in `x`. For example, the sign of `R'(x0)` determines how
an `x`-contour near a simple real zero approaches the radicand cut. Callers supply
the complex contour explicitly and are responsible for its physical homotopy.
`RootSeed::Value` uses a finite nonzero supplied value solely to select one of
the two signs; its magnitude is ignored after a scale-invariant angular sign check and recomputed
from `R` at working precision. Zero
or ambiguous hints are errors. This does not restore precision to a boundary.

Initial principal roots and endpoint germ references substitute the exact binary
rational coordinate in the exact radicand before converting it to MPFR. This
preserves exact cancellation and the side of the principal cut even when
separately rounded polynomial coefficients would choose another side.
At each accepted point, the next chart checks the radicands at exact binary
rational coordinates corresponding to the rounded MPFR center. Root coefficients
then use the existing rational differential-system recurrence for
`r' = (R'/(2R)) r`, starting from the accepted root values. This avoids constructing
large exact rational Taylor coefficients at each numerical center. Native
`Series` multiplication assembles the kernels with fixed-precision complex
arithmetic. Independent tests compare these root coefficients with Symbolica's
exact normalized `(R(x)/R(center))^(1/2)` series, including difficult rational
and complex charts.

Root series predict the analytic sheet, and an independently evaluated radicand
determines its magnitude. Both midpoint and endpoint root defects, root
truncation, and squared-root residuals must pass in addition to the shared
solution checks. Rejected trials leave branch state unchanged. The returned
`BranchState` records the endpoint, root identities in declaration order, and
values.

Root zeros and poles join rational matrix poles in the shared step-radius
constraints. This entrypoint transports between regular points. Algebraic
Frobenius endpoints and arbitrary algebraic extensions are not implemented here.
The separate [physical algebraic cache](algebraic-cache.md) supports progressive
`RustFlowCache` reuse on regular affine complex paths with certified root germs. The separate
[prescribed-contour planner](prescribed-contours.md) constructs detours from
explicit polynomial prescriptions; it does not infer the physical prescription
from a scattering process. The
returned `EpsilonSolution` has no `verified_digits` or independently verified
checkpoints: compile and transport again at higher precision/order to establish
accuracy. The physical cache wrapper performs this refinement and carries
inherited source uncertainty through root-aware weighted bounds. Explicit contour
transport must not be relabeled as verified cache evidence without those checks.

Regression tests cover one and two windings of `sqrt(x)` with
`y'=y/sqrt(x)`, giving `exp(-4)` and `1`, independent precision/order refinement,
direct epsilon coefficients, simultaneous transports on opposite sheets,
rejected-step nonmutation, sparse root-series defects, exact power normalization,
and invalid seed/registry inputs. Multivariate tests compare the exact pullback
with an independent chain-rule expression and verify its direct epsilon solution
against the analytic logarithm of a letter ratio. The implementation is independently derived
from the chain rule and Taylor recurrences; it does not translate DiffExp source.
