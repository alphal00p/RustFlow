# Endpoint arithmetic conditioning

Taylor transport now checks the conditioning of each candidate endpoint with
Symbolica's native `ComplexBall` polynomial evaluator. It assigns an explicit
relative one-unit working-precision perturbation to each stored real and
imaginary component of the retained coefficients and local argument. Exact zero
components have zero perturbation radius; a low-storage-bit zero imaginary part
does not degrade the real part. This avoids artificial uncertainty from large
coordinate scales and padded zero Taylor coefficients. The check compares the resulting radius on the same mixed scale as the
requested accuracy: `10^(-d) * max(1, |value|)`.

This catches severe cancellation that a zero numerical differential defect and
two agreeing rounded answers can miss. For example, the exact system
`y0' = (1 - 10^100) y1`, `y1' = 0`, starting at `[10^100, 1]`, has endpoint `[1, 1]`
after a unit interval. At 60 and 80 decimal working digits the unguarded port
returned the same incorrect zero for the first component. The conditioning
check requires more precision before this result can acquire accuracy evidence.

`Error::InsufficientPrecision { minimum_bits, context }` distinguishes this
condition from ordinary truncation, singular-domain, cancellation and resource
errors. The automatic epsilon transport and AMF/FT Laurent reconstruction wrappers
share an integer conversion of its conservative bit hint. They rebuild the system
and request fresh boundary data or epsilon samples within the existing attempt
limit. They still require two successful independent profiles.
Fixed-precision compiled transport returns the typed error to its caller. A
proposed physical step that rounds back to its center also requests fresh
precision before committing any values or branch state.

The automatic wrapper also retains exact waypoint expressions until checking
whether their numerical separation is resolved. It skips a provably identical
exact waypoint and preserves the identity behavior of an explicitly numeric
zero-length leg. The boundary API supplies only a floating-point coordinate;
its first point is compared as that exact dyadic value, without inventing an
original exact coordinate. A fresh provider can improve that coordinate during
a retry, but stale rounded inputs cannot gain missing digits. The low-level
compiled solver evaluates values, physical residuals, branch state and saved
segments at the same declared endpoint, including when a supplied coordinate
retains more bits than its working arithmetic. Midpoint sampling preserves the
input coordinate precision before averaging; it is not an exact-dyadic claim
for endpoints with arbitrarily separated exponents.

`FlowDiagnostics::conditioning_digits` records the minimum checked tolerance
along the accepted path and across both successful profiles. Saved-point checks
contribute their own ceilings. Stronger physical-cache evidence is capped by
this value as well as the existing comparison, source-accuracy and working-
precision constraints. Cached input uncertainty is still propagated separately;
raising arithmetic precision does not increase the source's evidence.

This is a conditioning diagnostic, **not** a proof that one unit bounds every
error incurred while constructing the Taylor coefficients. Native directed ball
arithmetic encloses the polynomial for the stated coefficient/argument balls;
the choice of those perturbations is heuristic. Ordinary rational charts also
check their residual against exact Gaussian-rational source coefficients with
native directed balls; registered-root residual assembly retains its separate
MPFR-estimate scope. Neither check is global forward-error certification. Very severe
cancellation in repeated nonlinear chart expansions can still exhaust the
bounded refinement budget; such runs fail explicitly rather than acquiring
verified digits from agreement alone.

The exact-source check preserves the existing local truncation tolerance, which
tightens with working precision. Enclosure widths for severely cancelling
non-dyadic coefficients can therefore remain too large through every permitted
retry. This is an explicit limitation of the current error budget; there is no
silent relaxation to the requested output tolerance. A residual whose central
value dominates the enclosure width still requests shorter steps. Supplied
numerical parameters are specialized as their stored dyadics, with external
uncertainty kept in the existing boundary/input contract.
