# Implementation coverage

This is an implementation status document, not a claim that the requested port
is complete. The original paper's four two-loop targets at s=30, t=-10/3, m²=1
remain the mandatory completion gate.

| Upstream operation | Rust implementation | Present limitation |
|---|---|---|
| Family setup, exact IBP | `family`, `reduction` | RustRed finite exact search; runtime arity 1–12 |
| Auxiliary differential equations | `differential_system` | Stable derivative closure required; residuals are not certified independent masters |
| SkipReduction | `differential_system_skip_initial` | Retains redundant target directions |
| RefineBasis | `refine` | Bounded D-factorizing swaps, reports if complete |
| Block structure | `DifferentialSystem::blocks` | SCC eigenvalue and recurrence LU blocks; full connection matrix |
| Ordinary propagation | `ode` | Pole-bounded Taylor steps with tail checks |
| Singular endpoints/infinity | `frobenius` | Rational indicial roots; integer shearing, Jordan-chain projector refinement, resonant logarithms |
| Auxiliary-mass placement | `MassMode`, `IntegralFamily::deform` | Intrinsic-mass groups, propagators, branches, loops, all lines, explicit positions; ordinary quadratic lines |
| Region enumeration/expansion | `regions` | Ordinary rank-one quadratic propagators; bounded routing budget |
| Tensor reduction | `tensor` | Repeated single hard vector: rank 32; other tensors: rank 14 and at most 128 pairing orbits |
| Massless bubble subloops | `bubble`, native adapter | Tensor rank at most eight, sum of bubble powers at most 32; transfer square must be external or an active denominator |
| Partial fractions/ISP completion | `integrand` | Affine scalar-product denominators, bounded term budget |
| Automatic boundary recursion | `recursive` | Depth 32; adaptive boundary order at most 16 in inverse mass; half-power numerator series |
| Feynman trick | `ft` | Recursive parameter DEs, regular rational reference points, Gaussian tensor terminals; Euclidean sectors only |
| Epsilon reconstruction | `epsilon`, `engine` | Conservative -2L pole bound; independent sample/precision refinement |
| Cache restart | `cache` | Reduction, native RHS-search checkpoints, and prepared AMF/FT systems; numerical memo is in-memory |

Missing for full parity: the full published four-target validation and
broader multiloop/complex-kinematics coverage. Pole-disk contour planning is
implemented and tested on both logarithm branches. Exact complex-mass tadpoles
are tested, but this is not general complex-kinematics validation.

A backend search budget must never be interpreted as a proof that a residual is
a master. The implementation checks repeated reduction and differential closure
and requires sufficient consistent boundary constraints. It returns typed
errors for incomplete closure, missing table rules, cycles, exhausted budgets,
unsupported normalization, nonfinite arithmetic, or failure of precision
refinement. Native RustRed calls are checked for cancellation before and after
the call; they cannot currently be interrupted from within this adapter.

## Numerical checks

The tests cover physical-region massless and massive bubbles, an on-shell
massless triangle, Euclidean massive bubbles, a non-vacuum single-mass sunrise compared against the
pinned upstream result with independent precision refinement, independent top
sectors, single-mass and two-mass two-loop sunsets with automatically generated region
boundaries, a three-loop single-mass banana vacuum with recursive FT boundaries, factorized
tadpoles, equivalent momentum routings and denominator permutations, linear
numerator boundaries, recursive FT sunsets, supplied-boundary analytic ODEs,
resonant logarithms, successive nilpotent balances, and Laurent fitting. Full
Laurent accuracy is tested with increased
precision/order and a second sample grid. Reference fixture digits are never
used as solver boundary data.

The test suite requires no Mathematica installation. Upstream regeneration is
optional and would require Mathematica and an upstream-supported reducer.

## Runtime limits and exceptional samples

The native runtime bridge accepts 1 through 12 scalar-product slots. Its default
adapter search uses depth 2 and at most 4096 concrete integrals; differential
closure uses at most 12 reduction rounds. RustRed's compact integral indices
range from -64 through 63, including generated seed indices. Exact rational epsilon specialization
is enabled by default for high-level multiloop AMF evaluation. It is covered by
a test comparing the resulting dimensional sectors with an analytic sunset.
Enabling basis refinement retains symbolic epsilon before factorization, even
when sampled reduction is otherwise enabled. Refined `new_at_epsilon` flows
therefore remain reusable at other epsilon values.
Recursive multiloop boundary families follow the sampled-reduction option too;
their memo keys include the exact epsilon sample. Boundary calls at exceptional
samples where dimensional exponent classes collide retain symbolic epsilon.

Nilpotent normalization first searches small subsets of Jordan-chain
projectors, then applies exact Moser projector refinement. Large chain counts
skip subset enumeration. At most 32 balances are performed, and a balance must
decrease the pole-order/rank pair. Unsupported indicial roots and an unsuccessful
normalization return explicit errors. Automatic AMF
rejects epsilon samples that create extra resonances between symbolic exponent
classes. This guards against an observed endpoint-sector ambiguity at D=1.

Native factorized replay is the default coefficient backend. It preserves source
and denominator nonzero conditions, performs exact substitution before searching
the next RHS frontier, and can checkpoint completed work within a search round. A larger
`max_targets` budget can resume that checkpoint. Finite-depth search can still
produce expensive intermediate rules; a checkpoint or a completed reduction
round does not imply that the differential basis has closed.

The opt-in `parametric_rules` backend option reuses guarded formulas with symbolic
positive powers and fixed numerator indices. Each specialization checks its
domain, preserves nonzero conditions, and must decrease the native integral
order. Uncovered cases fall back to concrete searches, and all RHS integrals
still require closure. See the [reduction audit](reduction-search-audit.md) for
the paper-example search bottleneck and the scope of prototype measurements.

For large native frontiers, coefficient expansion can cost more than the extra
IBP searches needed by a conservative dependency traversal. The
`max_exact_frontier` option controls this tradeoff (512 by default). All reported
final reductions still undergo exact substitution and a check that every
surviving residual was searched. Taylor transport checks the differential-equation
defect at step endpoints and midpoints as well as its final series terms, and
its truncation tolerance tightens when guard precision increases.

The native adapter also reduces eligible massless two-point subloops by exact
beta-function identities before requesting IBPs. These rational identities use
the same integral order as the sector solver, and each proposed rule is checked
to strictly decrease that order. Unsupported routings, complex coefficients,
exceptional dimensions, and nondecreasing rules fall back to IBPs. Set
`bubble_subloops: false` to disable this optimization. Checkpoint keys distinguish
the changed ordering; legacy partial checkpoints are reused only if no already
searched sector changes ordering.

`max_sector_batch` bounds the number of targets sharing one native exact
elimination (unbounded by default). Smaller batches can reduce intermediate
algebra at the cost of repeated seed searches. They retain the same strict
integral order and final exact substitution. `SectorReduced` progress events
report seeds, generated rows, exact-trace size, and total/exact milliseconds.
The acceptance example accepts `--case-batch N` for profiling this tradeoff.

Native search checkpoints use a versioned binary codec with authenticated
payloads and atomic replacement. Legacy JSON checkpoints remain readable;
new writes use `.bin` files. An incompatible or corrupt binary checkpoint is
rejected rather than silently replaced by an older JSON snapshot. Completed
reduction and differential-system caches retain their JSON containers.
With native checkpoints enabled, completed batches are saved every 60 seconds
by default, at search completion, and when cancellation is observed. Configure
`checkpoint_interval` to change this cadence. Only completed searches enter the
saved visited set; unprocessed targets and new right-hand-side candidates remain
pending. Resuming a partial round may search a conservative superset before the
next exact substitution.
Short intermediate rounds share that interval; rewriting a large checkpoint
after every small frontier change can otherwise dominate the actual search.

Increasing native search depth by one can also resume a shallower checkpoint.
Exact identities are retained, but every surviving residual is searched again
at the new depth. This matters even for a derivative-closed system: a redundant
basis can have spurious indicial modes. The physical paper subsector with lines
1, 3, 5, and 7 needs depth three in the current backend; depth two is rejected
by dimensional-sector validation.

`native_workers` enables bounded concurrent native batches (default one,
maximum 64 per reduction call). Results are merged in the same order as serial
search, and checkpoints retain only merged work. In-flight work is bounded by
one wave of batches. Account for both this setting and `FlowOptions::workers`
when budgeting total CPU and memory; the acceptance example exposes it as
`--native-workers N`.

Large reduction tables reuse Symbolica state headers during serialization and
reuse imported state maps during loading; entries remain readable in the
existing format. Conservative frontier traversal borrows graph keys instead
of allocating a key for every edge. Final exact substitution groups incoming
rational coefficients in balanced sums to reduce intermediate polynomial
growth; this changes evaluation order, not the reduction identities.
When at most 128 terminal integrals remain, substitution runs from the leaves
toward the targets so apparent poles can cancel within each identity. Child
expansions are released after their last parent uses them. Wider frontiers use
forward accumulation to avoid storing a large residual map at every graph node.
Small-basis substitution shares dependencies across a batch of targets, and
completed frontier expansions are reused when producing the final table.
`DifferentialClosure` events report basis size and remaining derivative targets.

Successive native searches at the same sector threshold reuse each target's
exact weighted frontier, including searched residuals. Appended identities are
applied to those coefficients; changed existing rules invalidate the cache.
Entering a lower sector also clears it, and the original targets are expanded
again to recover discarded lower-sector contributions. The complete dependency
graph is checked for cycles and sector increases before reuse, including paths
that previously cancelled. This cache is in memory only and leaves checkpoint
formats and final closure checks unchanged. Regressions compare incremental and
fresh expansions across partial rounds, exact cancellations, changed parameter
sets, rule replacements, wide frontiers, and cycles. See the
[controlled replay measurements](performance.md#incremental-frontier-replay)
for the measured scope of this optimization.

Before coefficient propagation, the native adapter installs zero identities for
all reachable leaves whose sectors RustRed has already certified as scaleless.
This includes lower sectors that have not reached the search queue yet. The
certificate domain conditions are retained, and full dependency-graph validation
still runs first. The same pruning applies when final substitution resumes from
a completed checkpoint; uncertified sectors remain subject to ordinary search.

Pole finding normalizes each factor and scales its variable using a power-of-two
Fujiwara bound before calling Symbolica's root finder. Newton corrections and
polynomial reconstruction check the returned root set. Regressions include a
degree-12 factor from the paper's boundary reduction and roots at very large and
very small scales; unresolved clusters return a numerical error.
Ordinary Taylor propagation clears denominators exactly per matrix row and
recurs over finite polynomial coefficients, omitting structural zeros. Tests
compare these coefficients against the independent rational-series recurrence
at complex centers. Tail and differential-equation defect checks are unchanged.
The pinned Symbolica `main` revision also checks root corrections and relative
coefficient backward error, fixing the coefficient-scale-dependent convergence
report reproduced in `repros/symbolica-root-scaling`. Scaling and root-set
validation remain in the AMFlow adapter. Cache keys include the manifest and
lockfile fingerprint, so dependency updates invalidate previous prepared systems,
reduction tables, and native search checkpoints.

Vacuum tensor projection groups pairings under permutations of repeated hard
vectors. At rank eight, four copies of each of two hard vectors need a
three-dimensional invariant system instead of the full 105-dimensional pairing
matrix. General distinct vectors retain the full projection. Analytic rank-eight
contractions and mixed-vector contractions test the grouping. Frobenius
recurrences traverse nonzero series-matrix entries without allocating a fresh
numerical zero for every comparison.
Single-vector projection uses the exact angular-moment formula with a memoized
Wick pairing sum, bounded to 100000 states. This extends the boundary tensors
beyond rank eight without constructing the full labelled-pairing matrix.

Boundary orders are selected by the rank of coefficient constraints in the
Frobenius solution. The planner requests all overlapping regions at each selected
power, retains uncomputed regions as unknown, and uses known-zero power and
logarithm constraints. It avoids expanding every component to the first nonzero
term of every fundamental solution. `BoundaryPlan` progress events report the
number of selected region series, coefficient count, and largest half-order.
Boundary matching still checks every computed coefficient, and numerical
refinement remains required to establish accuracy.

Factorized boundaries collect exact integral products before applying the term
budget. A regression from the paper's six-line boundary family at half-order
eight reduces over 10000 monomial contributions to 75 distinct products.
This removes a demonstrated boundary-budget failure; it does not establish the
full four-target numerical acceptance result.

`RecursiveBoundary::evaluate_many(family, targets, epsilon, precision)` evaluates
related boundary powers together. Region construction first collects the factors
needed by all selected coefficients, removes products with a scaleless factor,
and deduplicates requests within identical family definitions. Analytic and
user-supplied terminal values retain priority. Remaining requests share one
`PreparedFlow` per auxiliary-mass placement; vacuum targets choosing different
massive lines remain in separate groups. Prepared flows retain the complete
ordered target set in their cache key. Their target reduction coefficients are
multiplied into the endpoint Frobenius series before selecting the physical
limit, so poles in these coefficients are retained.

This follows the grouping in `BoundaryIntegrals` and `ReduceBoundary`, followed
by recursive evaluation of boundary master lists in `AMFSystemsSetup`, in the
[pinned upstream AMFlow.m](https://gitlab.com/multiloop-pku/amflow/-/blob/26005517a288086c4cb4d1b26d829691bc088485/AMFlow.m).
The batching regression checks preparation/backend-call savings for distinct
powers, differing vacuum mass choices, independent tadpole-product values,
and stability between 60 and 80 working decimal digits. Other checks cover
custom terminal providers, scaleless requests, duplicate targets, and reuse of
canonical numerical memo entries after denominator permutations. A connected
two-mass sunset batch with powers `111`, `211`, and `121` checks recursive
siblings against scalar solves and an independent beta-integral identity,
including a fresh evaluation at higher precision and expansion order.

Batching currently requires identical family definitions; it does not merge
families through a new momentum-routing or denominator-completion transformation.
The numerical memo remains per canonical integral, epsilon sample, and precision.
Recursion limits count nested batches rather than the number of integrals in a
batch, and encountering an active ancestor integral still returns a cycle error.
