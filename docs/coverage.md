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
| Tensor reduction | `tensor` | Exact vacuum projection, rank at most eight |
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
