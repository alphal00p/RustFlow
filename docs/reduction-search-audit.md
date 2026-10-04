# Paper-example reduction search audit

The required four-target, 20-digit numerical gate is still incomplete. The
historical full run stopped during IBP reduction, before numerical propagation.
This document separates that search bottleneck from the boundary-series work.

## Historical checkpoint census

The retained first-epsilon (`epsilon = 1/2700`) depth-three checkpoint
`native-stage-7efec0e2d85461458d8106175294642b05551a48bbd504825fedbabbaef752cf.bin`
has 717,122 completed concrete searches and 105,582 pending integrals. Its size
is 18,559,733,623 bytes. These are observations of an old, interrupted run, not
timings or coverage claims for the current Symbolica revision.

| Active propagators | Completed searches |
| --- | ---: |
| 7 | 76 |
| 6 | 3,323 |
| 5 | 56,409 |
| 4 | 375,656 |
| 3 | 281,658 |

Five three-line sectors dominate the pending work (line numbers are one-based):

| Sector | Pending integrals |
| --- | ---: |
| 2, 5, 7 | 32,584 |
| 3, 5, 7 | 23,170 |
| 3, 4, 7 | 18,687 |
| 1, 5, 7 | 17,954 |
| 1, 6, 7 | 13,077 |

The requested numerator rank is three. The reachable concrete searches have
already reached numerator rank eight; pending integrals have up to 15 extra
denominator powers. The seven-line sector itself is a small part of this work.

## Why the frontier grows

The adapter calls RustRed's `solve_numeric_cases` on small groups of concrete
integrals. Each call shares a modular search and exact replay within that group,
then returns local exact rules. The adapter recursively searches their RHS
integrals. Its support-first order makes pinches simpler even when the numerator
rank or extra denominator powers increase. Changing only the final coordinate
tie order does not change that property.

The historical run used `max_exact_frontier = 0`, which deliberately postpones
exact coefficient cancellation until the structural frontier has been searched.
This avoids enormous intermediate rational expressions but explores every raw
RHS dependency, including contributions that later cancel. More workers or a
larger target budget cannot by themselves resolve that algorithmic tradeoff.

The analytic bubble shortcut cannot eliminate the dominant `(2,5,7)` sector
inside the same integer-power family. Integrating its massless loop-one pair
leaves transfer `(l2-p1)^2`, which is neither external-only nor the remaining
active denominator `(l2+p3)^2-1-eta`. The current shortcut explicitly requires
one of those forms. This rejection prevents introducing undeclared or
noninteger-powered propagators.

## Available RustRed routes

`solve_numeric_cases` currently offers sparse and factorized sparse exact
replay. The adapter already selects the latter. `solve_domains` can instead
search formulas with symbolic positive powers while fixing numerator indices,
including exceptional branches under explicit limits. Such formulas could
amortize many concrete searches, but their guards, source conditions, descent
and recursive RHS coverage must all be preserved when applying them.

RustRed's `CandidateReducer` provides guarded formula application but its current
family admission requires a unit-mass vacuum shape. It cannot be substituted
directly for the nonvacuum paper family. Its finite residuals are also explicitly
not certified independent masters.

The isolated `paper_reduction_probe` example exercises the native numeric and
directed-domain search APIs on representative three-line targets. It reports
preparation, search and exact replay times, raw RHS breadth and index growth.
The `terms` source-order option changes the order in which existing IBP rows are
visited; it does not change the integral order or the source identities. This
probe does not recursively close the RHS and always reports
`complete_reduction: false`.

## Comparison scope

Pinned upstream AMFlow delegates IBPs to external reducers; its documented
default is FiniteFlow plus LiteRed, with initial seed rank three and zero dots
(raised when targets demand more). That workflow constructs and reduces a
shared IBP system, unlike the adapter's repeated local searches. Consequently,
native local-search timings must not be presented as end-to-end AMFlow speedups.
A valid full comparison must include reduction, boundary construction and
numerical solving at the same target accuracy, and separately identify cached
and cold runs.

## Guarded-ray prototype measurements

An isolated diagnostic on the updated Symbolica dependency generated a ray with
fixed numerator rank three for each of the five sectors above. Each ray left all
three positive powers symbolic. Generation took 3–22 ms and produced four to
seven guarded rules. These are successful searches for the nominated rays, not
proofs that their recursive descendants are covered.

A matched-work prototype in sector `(2,5,7)` used all 125 triples of positive
powers from one through five, with the same fixed rank-three numerator:

| Operation | Observed time | Other observations |
| --- | ---: | --- |
| Generate one guarded ray and specialize all 125 targets | 29 ms | 125 applicable specializations |
| Shared native concrete search for the same 125 targets | 299 ms | 176 ms exact replay; 927 replay rows; 306 distinct raw RHS integrals |

These single-run measurements used an `opt-level=1` standalone diagnostic linked
to optimized dependencies. They motivate the integration; they are not a
production benchmark, a full reduction comparison, or an AMFlow speedup claim.
Both outputs still require recursive RHS reduction.

The library now exposes `RustRedBackend::parametric_rules`, initially opt-in. Its
rule bank checks fixed coordinates, coupled affine equalities, complete
exception conjunctions, original denominator guards, source conditions, compact
index bounds and exact descent. It retains specialized nonzero conditions in
`Reduction`. Missing or exhausted rays fall back to the concrete search, and
finite search residuals never become rules. `ParametricReduction` progress events
report generated rays, applied formulas and uncovered targets. The bank is
scoped to a reduction's source family, sector and source ordering; it is currently
an in-memory optimization, while specialized concrete rules use the existing
restart checkpoints.
