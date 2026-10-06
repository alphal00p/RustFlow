# Optional rational trial transport

`FlowOptions::pade = Some(PadeOptions::default())` enables native Padé trials
in ordinary rational and epsilon-expanded differential systems. The CLI option
is `"pade_degree": 16`; Python uses `EvaluationOptions(pade_degree=16)`. Omitting
this option retains Taylor transport. Registered-root systems explicitly keep
Taylor trials, including their branch tracking.

The shared controller proposes physical steps using its existing source-pole,
coordinate, and step-budget rules. It first tries a rational approximation at
the same proposed point. A construction, domain, defect, or conditioning failure
selects Taylor at that point; a cancellation aborts the operation. Diagnostics
record rational trials, committed rational steps, fallbacks, and the last reason.
This option does not extend the source convergence radius or bypass a source hole.

Symbolica's native rational approximant constructs each component from exact
stored MPFR dyadics. Native univariate GCD removes common polynomial factors;
we verify degree bounds, a nonzero denominator at the chart center, and exact
Taylor matching. A common denominator is formed with native GCD and products.
For a source row `D Y' = A Y` and candidate `P/Q`, the retained exact numerator
of the differential defect is `D(P'Q-PQ')-Q A P`. All source degrees and epsilon
couplings are retained, including in a Möbius chart. Directed enclosures bound
this defect on the whole trial disk and exclude zeros of `D Q²`. Endpoint and
midpoint differential checks remain independent admission checks.

The arithmetic-conditioning diagnostic evaluates the rational function with
native balls around its exact coefficients and argument. The existing Taylor
seed conditioning check is retained separately. Saved segments retain their
accepted rational representation; `EpsilonSolution::evaluate_segment` and
intermediate cache insertion therefore evaluate the accepted function. The
public `TaylorSegment::coefficients` are seed coefficients when `pade` is present;
call `evaluate_local` to evaluate that segment's actual candidate.

These are local defect and arithmetic-conditioning checks, **not a rigorous
global forward-error proof**. Boundary uncertainty, transport amplification, and
independent order/precision comparisons retain their existing roles. Every fresh
comparison profile raises the requested rational degree as well as Taylor seed
order; if the hard degree cap would be exceeded, that reference uses Taylor. Fresh
precision retries must regenerate source and boundary data. No fitted coefficient
is certified merely by the working precision or by agreement of two Padé degrees.

`PadeOptions` bounds component degree (at most 32), common denominator degree,
input representation height, estimated intermediate coefficient height, and
estimated native polynomial work. The binary exponent and stored precision are
checked before converting a Float to a rational. Exact source terms are checked
before shifts. Cancellation is checked between bounded native operations; one
native polynomial operation is not asynchronously interruptible. The work limit
is a conservative construction estimate, not a wall-clock or allocator quota.
Large or ill-conditioned common denominators can trigger Taylor fallback. This
first rational/epsilon slice is not a performance-parity claim for large systems.

Point-bank files persist values and accuracy evidence, not approximants, so their
binary codec is unchanged. Source-sensitive identities still include the full
native source and dependency graph; completed finite-epsilon sample keys also
include these options. Existing caches do not become independent validation.
