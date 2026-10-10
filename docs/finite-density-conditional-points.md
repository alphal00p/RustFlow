# Native reduction on exceptional index loci

A valid recurrence can fail at particular integer indices. For example, the
massless two-loop double cut exposed a recurrence with denominator `a3-a4`, while
the requested numerator has `a3=a4=-1`. Removing that condition would invalidate
the reduction. The full original source corpus instead supplies an independently
replayed recurrence at that exact point, whose remaining denominator is `eta`.

Residual-focused discovery now distinguishes an ordinary `NoApplicableRule`
leaf from a native `ConditionVanished` failure. Ordinary residuals retain their
existing ray searches. A vanishing guard requests a singleton integer box from
the same native source owner. This also handles coupled exceptional loci that
cannot be represented by fixing just one coordinate of a ray.

Both searches share the existing total and per-domain discovery allocations.
Newly exposed guards require an available refinement pass. An exact point is
attempted at most once per discovery call. Diagnostics record its indices,
pass and conservative allocation separately from finite-axis subdivisions.

Searching a point establishes no coverage by itself. Successful older rules
retain precedence, and a new point rule must pass native source replay before
it can supply a replacement. Target and derivative application is repeated
against the resulting program. An unresolved guard remains a failure; it is
never cancelled between targets, installed as a terminal, or declared zero.
Other native failures retain their existing failure paths. The new rule's
physical parameter conditions remain part of the eventual flow contract.

The retained failure evidence and native point diagnosis are in
[`2026-10-10-native-zero-domain-projection`](../reports/validation/2026-10-10-native-zero-domain-projection).
Its `pre-refinement-*` artifacts describe the earlier failing build. Subsequent
`guard-point-*` artifacts use a separately frozen build and explicitly record
their refinement budget; they do not rewrite the earlier experiment.

The revised build passes 179 native/RustFlow regression tests (one additional
native test remains ignored). Full massless two-loop predictions pass 40
fixed-dimension reference checks and 30 refinement checks, followed by 30
Laurent reference checks and 24 refinements. The massive two-loop single-profile
regression also agrees with its independent reference and previous predictions.
Both newly exposed double-cut diagonal points are resolved by native point
proofs. The three-loop input still exceeds the 1,024-label frontier budget; no
three-loop numerical acceptance follows from this fix.
