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

## Bounded production pilots with guarded rays

These pilots use the archived optimized library at commit
`ce21f7693ae7d642006cb0896efda25eb2d80f20`, pinned to Symbolica main
`75f8350094b90254ee71dc2a391fde0d14b0204a`. They attempt the first paper target
`[1,1,1,1,1,1,1,-3,0]` at `epsilon = 1/2700`, or construct that family's
recursive boundary regions. Each invocation runs on one pinned CPU with a
300-second wall limit. These are diagnostic runs, not numerical acceptance or
complete AMFlow timing comparisons.

The first-target runs use the archived `paper_subsector` executable, depth
three, and `--parametric-rules`. The regions runs use the same validated library.
The forced-pruning regions executable adds only command-line limit controls to
`paper_regions`; that wrapper was compiled separately against the archived
library, rather than being part of the archived source commit.

| Pilot | Exact-frontier threshold | Target budget | Outcome |
| --- | ---: | ---: | --- |
| First target, structural expansion | 0 | 8,192 | Limit: 4,612 searched + 8,558 pending |
| First target, small exact threshold | 512 | 8,192 | Limit: 4,362 searched + 8,220 pending |
| First target, forced exact cancellation | 100,000 | 8,192 | Limit: 8,135 searched + 93 pending |
| Resume forced cancellation | 100,000 | 16,384 | Five-line work closed; limit: 8,550 searched + 22,943 pending |
| Resume into four-line work | 100,000 | 65,536 | Wall timeout during substitution; checkpoint: 35,271 searched + 3,407 pending four-line targets |
| Regions, structural expansion | 0 | 32,768 | Limit: 11,683 searched + 21,156 pending |
| Regions, forced exact cancellation | 100,000 | 32,768 | Wall timeout during substitution; checkpoint: 7,245 searched + 2,679 pending four-line targets |

The small threshold skips exact pruning once a frontier exceeds 512, so it
barely affects growth. Forced cancellation closes the active five-line layer
but then exposes 20,621 four-line, 2,321 three-line, and one two-line target.
It postpones lower layers and cancels some coefficients; it does not by itself
solve the breadth of the next layer.

During the 65,536-budget continuation, guarded rays supplied 28,066 reductions
from 1,704 generated rays, with 13 concrete fallbacks. Ray generation and
specialization consumed 57.4 seconds in the reported events; concrete fallback
search took 1.0 second. The timeout occurred while substituting 38,649 accumulated
rules, after the checkpoint above. The checkpoint is 69,783,278 bytes; a live
memory observation was about 585 MiB. Thus this run shifted the main cost from
repeated concrete searches to exact coefficient propagation. The terminal
checkpoint excludes work completed after the last atomic save.

In the forced-cancellation regions pilot, 653 generated rays supplied 9,897
reductions, with 27 concrete fallbacks. Reported ray work took 89.0 seconds and
concrete searches 4.5 seconds. Its final event was substitution of 9,897 rules;
the last checkpoint is 27,022,535 bytes. A live memory observation near the
wall limit was about 654 MiB. Neither production pilot completed differential
preparation, boundary initialization, or a numerical sample.

Logs and checkpoints are retained under the ignored task directories
`target/parametric-paper-{pilot,exact-pilot,pruned-pilot}` and
`target/parametric-regions-{pilot,pruned-pilot}`. Restart counts are cumulative;
per-invocation ray/search times are not summed across different experiments.
