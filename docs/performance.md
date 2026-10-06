# Comparison with original AMFlow 2.0

A same-host comparison against the original, unmodified C++ differential-equation
solver is available. It measures regular continuation and singular matching with
supplied boundaries. A first complete automatic finite-epsilon comparison also
covers a connected three-loop single-mass vacuum, with the settings differences
recorded below. Broad matched automatic-workflow performance remains unmeasured.
A licensed Wolfram runtime is now available through the Bern server, and
[live automatic oracle checks](upstream-oracle.md) exercise the original package
with Kira. Earlier small oracle checks had no corresponding Rust timing; the
new vacuum comparison includes both complete workflows without repetitions.
The Rust four-target two-loop
acceptance has passed; its numerical validation and execution conditions are
recorded separately below.

## Automatic single-mass vacuum

Both implementations evaluated `111111` and `211111` for the six denominators
`[q²−3, l², (l−q)², k², (k−q)², (l−k)²]` at exact epsilon `1/10`, using automatic
IBP reduction, boundary recursion and transport. Each process began with a fresh
symbolic cache, one worker, 40 requested digits and 80 configured working digits.
Both outputs pass a 35-digit relative comparison against the original's recorded
40-digit precision.

| Workload | RustFlow process | Original AMFlow + Kira process | Result comparison |
|---|---:|---:|---|
| Both three-loop targets, epsilon=1/10 | 2.814 s | 175.447 s | Both pass 35 relative digits |

This is one observation per implementation on the same shared host, pinned to
CPUs 32 and 33 respectively. Both used order 160; original AMFlow also used its
separate 50-term `ExtraXOrder`, for which RustFlow has no identical control.
The resulting top bases contained three native and two original integrals.
Native preparation/reduction took 0.143 s and automatic boundary/transport
2.627 s. OS page caches and other host activity were uncontrolled. These results
do not establish a general speedup or replace the larger benchmark timings.

At this exact epsilon, native sampled reduction could not uniquely recover an
indicial sector and returned `Unsupported`. The successful run retained symbolic
epsilon during preparation, then evaluated the requested finite sample through
the public API. The [full report](../reports/performance/single-mass-full-workflow.json)
preserves that failed trial, both successful outputs, configuration differences,
source/binary hashes and all timed phases. No reference value was used as a
native boundary.

## Latest ordinary-transport measurements

Commit `2a1ffac` was measured on 2026-10-06 with the same driver, exact inputs,
CPU 35, original executable and independent precision/order checks as the
preceding comparison. All three regular systems pass 20 scaled digits. Each
entry is a median of five fresh processes after an excluded warmup, alternating
implementations; base settings are 201 bits/order 80 and refinement uses
267 bits/order 112.

| Workload | Rust process | Original process | Rust transport | Original transport | Accepted steps, Rust / original |
|---|---:|---:|---:|---:|---:|
| Analytic 3-component chain | 0.09119 s | 1.75172 s | 0.082911 s | 0.001927 s | 9 / 4 |
| Upstream 12-master system | 3.76605 s | 4.88680 s | 3.649177 s | 0.102201 s | 11 / 5 |
| Paper 27-integral system | 3.96480 s | 14.67071 s | 3.745475 s | 0.061878 s | 4 / 1 |

The [complete report](../reports/performance/2026-10-06-ordinary-transport.json)
retains outputs, source hashes, settings, runtime provenance and all comparisons.
The largest scaled differences remain `8.06e-28`, `6.40e-43` and `1.45e-21`.
Native rejections are 8/8/3. The faster 12/27-master results accompany fewer
trial charts after changes to denominator enclosures; this comparison spans
multiple numerical changes and does not isolate a single optimization.

The native ordinary transport phase remains slower than upstream. Original
process totals include MPSolve setup, and neither these supplied-boundary jobs
nor their process ratios measure automatic AMFlow or the gg→hg application.
Concurrent tests/builds and uncontrolled OS caches remain timing limitations.
A separate current 27-integral profile places 77.98% of sampled time in exact
source residual assembly and 70.15% in native ball-polynomial multiplication.
Those checks remain required; the profile identifies optimization work, not a
reason to weaken numerical acceptance. Local checks still do not certify global
forward error. The combined full release suite was running when this benchmark
report was captured; its three benchmark validations and strict Python/stubgen
Clippy had independently completed. The subsequent
[full release result](../reports/validation/2026-10-06-cache-arithmetic-release/summary.json)
passes 547 unit/integration tests and one doctest, with 13 opt-in tests ignored,
and verifies that the measured numerical sources remained unchanged.

An [exact-rational residual experiment](https://github.com/alphal00p/RustFlow/blob/d8738fc0a9f0e08572fbae84cef58f5d3378ef24/reports/performance/2026-10-06-exact-source-residual/report.json)
was evaluated separately and was **not integrated into production**. It reuses
Symbolica polynomials to assemble an exact defect before outward conversion.
Two 14-case runs passed 31,716 exact-component containment checks, and 21 focused
tests, formatting and strict release Clippy passed. Initial-chart gains did not
persist across genuine complex endpoint charts; adverse high-degree/rational
inputs made construction 3.7–6.2 times slower. Exact cancellation also removed
an existing precision-retry signal. These are local constructor measurements,
not end-to-end transport timings. The archived branch retains the prototype,
counterexamples, failed preliminary checks and complete measurement provenance;
the validated production ball-arithmetic path remains in use.

## Preceding checked release

Library commit `3f1daea` was measured on 2026-10-05 against the original,
unmodified AMFlow 2.0 C++ executable. All four supplied-boundary comparisons
pass 20 scaled digits, including independent precision/order refinement.
Ordinary transport includes exact-source directed defect checks and endpoint
conditioning. These remain local diagnostics rather than a global accuracy proof.

| Workload | Rust process | Original process | Rust transport | Original transport | Accepted steps, Rust / original |
|---|---:|---:|---:|---:|---:|
| Analytic 3-component chain | 0.09007 s | 1.72620 s | 0.083151 s | 0.001903 s | 9 / 4 |
| Upstream 12-master system | 5.49649 s | 4.71083 s | 5.382974 s | 0.102368 s | 16 / 5 |
| Paper 27-integral system | 32.60942 s | 14.34587 s | 32.403464 s | 0.061628 s | 33 / 1 |

These are five-run medians after excluded warmups, with alternating fresh
processes pinned to CPU 35, one worker, requested accuracy 20 digits, working
precision 201 bits and order 80. Independent refinement uses 267 bits and order
112. The largest scaled validation errors are `8.06e-28`, `6.40e-43` and
`1.45e-21`. The metric divides each complex difference by
`max(1, |left|, |right|)`. Native halving uses 17, 34 and 145 trial predicates,
including 8, 18 and 112 rejections. The original uses its own fixed step
prescription. The 27-integral boundary is a deterministic test vector; this row
does not measure automatic physical boundary construction.

Two-component resonant Frobenius matching takes 0.007303 s in Rust and 1.137460 s
upstream, with process totals 0.014015 s and 1.142808 s. All analytic, original
and refinement checks pass, with maximum scaled error `2.30e-25`. This case
uses singular matching rather than Taylor transport. Original matching uses at
least 450 decimal digits internally and its separate `ExtraXOrder=20`, so the
work differs.

The [complete report](../reports/performance/2026-10-05-conditioned-transport-checked.json)
preserves outputs, repetitions, phases, counters, settings, hashes, source-to-
commit verification and rejection of an obsolete cache snapshot. Other host
jobs and OS page caches are uncontrolled. Original process totals include
MPSolve startup and pole finding; transport and complete-process timings
therefore answer different questions. Current ordinary native phases are slower
than the original. The [preceding release measurements](../reports/performance/2026-10-05-whole-segment-checked.json)
remain preserved, including their earlier Nix-runtime setup repair.

A fresh [unequal-banana comparison](../reports/performance/2026-10-05-conditioned-banana-low-profile.json)
measures the verification-cost regression on the same driver/helper/fixture:
497.728 seconds of current evaluation versus 55.685 seconds in the preceding
native release. Both pass all 75 original coefficients at 20 absolute digits.
Physical transport now uses 127 accepted/246 rejected steps versus 83/108.
The current three-route integration test also passes, in 1309.87 seconds versus
156.96 seconds previously. These observations identify additional arithmetic
and changed step selection as performance work; they are not original-DiffExp
speed ratios. Original and larger-system comparisons are tabulated in the
[benchmark summary](benchmark-comparison.md).

## Measured baseline

The baseline uses Rust library commit `de5981dd18a3cd5ba5d56aaf8a45e977462a19aa`,
including the fixed Symbolica `main` revision, and original AMFlow commit
`26005517a288086c4cb4d1b26d829691bc088485`. The Rust benchmark driver was linked
against that already validated release library; its source and executable hashes
are recorded separately. Later changes are not included in this baseline.

Measured on 2026-10-04, AMD EPYC 9754, one process restricted to logical CPU 24.
Each implementation used one thread, 201-bit arithmetic (the C++ solver's
60-decimal-digit setting), series order 80, and identical exact matrix, epsilon,
boundary literals, and endpoints. Numbers below are medians of five fresh
processes after one excluded warmup. The order alternated between implementations.
Compilation is excluded; process startup, parsing, pole finding, continuation,
and output are included in total time.

| Workload | Rust total | Original total | Rust transport | Original transport | Accepted steps, Rust / original |
|---|---:|---:|---:|---:|---:|
| Analytic 3-component logarithmic chain | 0.1116 s | 1.6887 s | 0.10493 s | 0.001908 s | 9 / 4 |
| Upstream 12-master system | 1.6947 s | 4.6665 s | 1.58679 s | 0.102466 s | 7 / 5 |

The original solver's preparation includes launching MPSolve for pole finding.
Its transport is substantially faster despite the larger total time in these
short jobs. These are **not full-AMFlow speedups**, and repeated continuation
inside a prepared process would amortize setup differently. Rust also performs
adaptive tail and differential-equation defect checks; the original `RegularRun`
uses its fixed step prescription. The logarithmic-chain run rejected eight Rust
trial steps, while the 12-master run rejected none.

The minimum required accuracy is 20 digits, checked with arbitrary-precision
decimal comparisons. Both implementations were independently repeated at
267 bits (80 decimal digits) and order 112. Every measured repetition is checked
against the corresponding refined result, and both refined results are compared.
The analytic case additionally checks the exact logarithmic solution. Baseline
cross-implementation relative errors were at most `8.1e-28` and `6.4e-43`,
respectively. Accuracy is normwise relative with an absolute tolerance below unit
magnitude. The implementations achieve different accuracy beyond the requirement:
the analytic original result has roughly 27 stable digits, whereas Rust has
roughly 60. Identical settings therefore do not imply identical numerical work.

The 12-master matrix and boundary come from the pinned original DE example; see
[fixture provenance](../fixtures/performance/README.md). The benchmark transports
from eta=1/2 to eta=1/10+i/5 at epsilon=101/9999900. The boundary retains its
180-digit upstream precision annotations, but these are not interpreted as a
verified accuracy estimate. This supplied-boundary test does not exercise
automatic boundary construction, IBP reduction, or Laurent fitting.

[Raw results](../reports/performance/2026-10-04-baseline.json) preserve every
measured timing, output component, validation error, input, compiler/version,
thread setting, executable hash, and source provenance. The machine had other
work running; its initial load averages are recorded. Peak RSS is retained only
as an OS diagnostic, since Python's forked child can inherit a peak-memory floor.

## Polynomial-recurrence result

The optimization in commit `9f23cfaa31cd81c9c2f0e2b31af9b90683a6a126` clears each
row's denominators exactly before numerical evaluation. Taylor propagation then
recurs over sparse polynomial coefficients and their finite degrees. It no
longer expands every rational matrix entry to order 80 or multiplies all zero
matrix entries. Pole finding, path limits, adaptive tail checks, and endpoint
and midpoint differential-equation checks are retained.

The repeated comparison uses validated library commit
`ce21f7693ae7d642006cb0896efda25eb2d80f20`, with the same original executable,
driver source, harness, inputs, CPU 24, precision, series order, and five-repeat
procedure. Original AMFlow was rerun alongside the optimized Rust executable.
Both comparisons and both independent precision/order refinements still pass
20 digits.

| Workload | Optimized Rust total | Original total | Optimized Rust transport | Original transport | Rust transport improvement over baseline |
|---|---:|---:|---:|---:|---:|
| Analytic logarithmic chain | 0.01783 s | 1.71098 s | 0.010742 s | 0.001903 s | 9.77× |
| Upstream 12-master system | 0.27302 s | 4.69397 s | 0.160158 s | 0.102376 s | 9.91× |

The accepted and rejected step counts are unchanged from the baseline. The
12-master before/after result changes are at most `4.72e-57` at 201 bits and
`5.92e-77` at 267 bits under the same scaled-error metric. The analytic outputs
are identical at the printed precision. The original solver still propagates
faster: it uses five steps for the 12-master case, compared with seven in Rust,
and its error-checking policy differs as described above. These improvements
refer to Rust's measured transport phase, not to full automatic integral
evaluation or a claimed general speedup over original AMFlow.

Both Rust library builds use the same Cargo release profile. The baseline
driver was linked separately using `rustc -O`; the optimized driver was built
with Cargo's release profile. The reported improvement ratios use timers around
the library transport call, excluding the driver's preparation and output.
The after-run data and every before/after comparison are preserved in
[`2026-10-04-polynomial-recurrence.json`](../reports/performance/2026-10-04-polynomial-recurrence.json),
which also records the baseline report's digest. This is a two-workload result;
denominator clearing can increase polynomial degrees in other systems.

## Build and reproduce

The original solver requires GMP, MPFR, MPC, Boost, yaml-cpp, and MPSolve. The
builder downloads the exact upstream archive, verifies its SHA-256, uses the
repository's pinned Nixpkgs revision for dependencies, and builds `desolver`
without WSTP or Mathematica. It changes only dependency paths in a separate
configuration file; the upstream solver source is unmodified.

```sh
nix develop
# Configure the Symbolica license in your shell as usual.
cargo build --locked --release --example benchmark_de
python3 scripts/build_upstream_solver.py --output target/original-amflow-bench
python3 scripts/benchmark_de.py \
  --rust target/release/examples/benchmark_de \
  --upstream target/original-amflow-bench/amflow-26005517a288086c4cb4d1b26d829691bc088485/diffeq_solver/desolver \
  --upstream-build-info target/original-amflow-bench/build.json \
  --rust-source-commit "$(git rev-parse HEAD)" \
  --output target/performance-comparison --repeats 5
```

Use a clean worktree, or preserve its patch with the report. Each output directory
must be new. Add `--cpu N` to restrict both executables to an available logical
CPU, as in the baseline. The comparison rejects a C++ executable whose hash does
not match its build record. It returns failure if any accuracy check fails;
working precision alone is never treated as verified accuracy.

The measured original build used GCC 15.3.0 and its default release flags
`-std=c++17 -O3 -ffast-math -march=native -w`, with OpenMP and quad arithmetic
disabled. Nix's compiler wrapper removed `-march=native`; the build record and
log preserve this fact. Rust used release optimization, debug information level
1, and no target-native override. The build record includes exact dependency
store paths. Floating-point values used for integral comparison remain decimal
strings and Python `Decimal`; only timing counters use ordinary OS number types.

## What the profile shows

A separate 499-Hz CPU-cycle profile of the baseline Rust 12-master run collected
878 samples with no lost samples. Approximately 97% fell under transport.
Building Taylor expansions of matrix coefficients accounted for about 45% of
the run: about 38% in rational quotient expansion and 7% in polynomial shifts.
Allocation and MPFR initialization are substantial within the arithmetic calls.

The baseline Rust recurrence expands every rational matrix entry to the full
series order and performs a dense convolution, including zero matrix entries.
For 12 components and order 80, this is 466560 complex product-adds per step;
only 61 of the 144 matrix entries are nonzero. In contrast, the original solver
clears denominators per block and recurs over the finite polynomial degrees
(`diffeq_solver/src/mpsolver.cpp`, `regular_expansion`). This profile motivated
the row-wise denominator clearing measured above. The different step counts and
error checks account for additional differences; no single cause is inferred
from total timings alone.

## Incremental frontier replay

A separate exact-algebra prototype measured reuse of filtered target coefficients
between native reduction rounds. It used a frozen paper-example checkpoint with
40840 rules and 1471 residuals, targeting `[1,1,1,1,1,1,1,-3,0]`. Twelve existing
rules at one sector threshold were withheld, then restored in three groups of
four. Each update was evaluated both from the original target and from the
preceding weighted frontier; every coefficient difference was exactly zero.
These controlled rounds use real identities but do not reconstruct the original
search chronology.

| Retained sector threshold | Required initial expansion | Fresh updates, total of three rounds | Cached updates, total of three rounds | Update-time ratio |
|---|---:|---:|---:|---:|
| Five active lines | 6.361 s | 19.111 s | 0.6755 s | 28.3× |
| Four active lines | 86.891 s | 268.651 s | 2.0117 s | 133.5× |

These are single controlled sequences, not repeated-run medians. The standalone
prototype used `rustc -O`, the pinned Symbolica dependency, exact polynomial
arithmetic, and one pinned CPU per case (28 and 27, respectively). Other work
ran concurrently. Timers include constructing the coefficient field and exact
substitution, but exclude checkpoint loading, searches, full-root structural
validation, boundary construction, propagation, and epsilon fitting. The initial
expansion is still required, including after a restart. The ratios therefore
measure repeated substitution only; they establish no end-to-end speedup over
Rust's previous workflow or original AMFlow.

The implementation retains searched residual coefficients, invalidates cached
frontiers when an existing rule changes or the sector threshold decreases, and
recovers lower-sector contributions from the original roots. Full dependency
validation also rejects cycles inside previously cancelled branches. Checkpoint
formats and final exact closure checks are unchanged.

[Replay metadata and individual timings](../reports/performance/2026-10-04-incremental-frontier-replay.json)
record the local prototype/checkpoint hashes and the controlled rule selection.
Those large local artifacts are not included in the report; this is diagnostic
evidence motivating the optimization, not a portable benchmark suite.

## Configurable backward substitution

A frozen paper-target checkpoint retained 544 structural terminal integrals
after certified scaleless leaves were removed. Raising `max_backward_frontier`
above its default of 128 lets this graph use backward substitution, simplifying
each child before substituting it into its parents. In one bounded run, backward
expansion finished in 24.247 s (27.009 s total, 507 MiB peak RSS); forward
substitution on the identical graph reached its 600 s limit without completing
(600.683 s total, 1,590 MiB peak). The forward time is censored, so exact output
agreement is unavailable and no speedup ratio is claimed. Backward substitution
produced 125 nonzero terms, including 87 three-line terms still needing searches;
this is not a completed reduction or numerical acceptance result. Other work
ran on the host. [The experiment report](../reports/performance/2026-10-04-backward-frontier.json)
records hashes, conditions, strategy settings and timing provenance. The result
motivates an opt-in threshold increase for this snapshot; wider backward maps
can cost more memory on other graphs.

## Full automatic workflow

The Rust calculation of all four paper targets through epsilon power zero passed
20-digit independent-grid stability and the upstream recorded-precision checks.
It used 27 samples at 60 working decimal digits/order 80, then 31 samples at
80 digits/order 96. The recorded elapsed time was 3553.007 seconds. This run
started with the exact caches from the earlier epsilon=1/2700 calculation, and
four additional processes prepared exact validation-grid systems alongside the
four numerical workers. It is therefore an assisted warm run, not a cold
four-worker benchmark. No numerical answers or reference coefficients entered
those caches. The [acceptance report](../reports/validation/2026-10-04-paper-two-loop-acceptance.json)
preserves the cache manifest, assistance details, all coefficients, observed
refinement errors, and executable/source hashes.

[`scripts/benchmark_full.wl`](../scripts/benchmark_full.wl) is an optional harness
for a configured Wolfram installation and upstream-supported IBP reducer. It
provides the physical massless bubble, two-mass sunset, and all four required
paper targets through epsilon power zero at 20 requested digits, one thread,
and caching disabled. It checks the pinned `AMFlow.m` digest and writes results
and timing metadata in the current directory. Use a fresh directory for each run.

This full timing harness has **not been executed** here; the separate
[oracle harness](upstream-oracle.md) has passed small automatic examples. The
full timing harness's output explicitly marks accuracy
unverified; results must be compared and independently refined before timing
claims are made. Rust automatic evaluation already performs additional sample
and precision refinement, while upstream `SolveIntegrals` uses its own numerical
configuration. A future full comparison must report this difference and include
equivalent correctness checks. Stored upstream answer files contain no timing
evidence and are not used to manufacture an original runtime.

## Singular matching with an integer resonance

A third matched workload constructs a Frobenius expansion at eta=0 for
`A = [[1/(1-eta), 0], [1, 1/(eta*(1-eta))]]`, using the exact boundary
`y(1/2) = (1, 0)`. Its solution is
`y = (1/[2(1-eta)], eta*log(2*eta)/[2(1-eta)])`. The indicial exponents 0 and 1
produce a resonant logarithm. Both implementations evaluate the matched series
at eta=1/10; Rust additionally extracts the physical endpoint `(1/2, 0)`.
The C++ input uses one rational fraction per matrix entry, as its parser requires.

The Rust driver was linked with `rustc --edition=2024 -O -C debuginfo=1` against
the archived release library from validated commit
`917886e777a31638504fd16a78eaf65081808626`. The original executable is unchanged
from the comparisons above. Both processes were restricted to CPU 24 with one
thread, requested 201-bit working precision, and nominal series order 80. These are medians of
five fresh processes after one excluded warmup, alternating implementation order.

| Phase | Rust | Original C++ |
|---|---:|---:|
| Total process | 0.014369 s | 1.162599 s |
| Singular expansion and boundary matching | 0.007174 s | 1.155660 s |

The Rust phase timer includes Frobenius construction and matching. The original
phase includes translation, indicial reconstruction, recurrence and matching.
Final series evaluation is excluded from both phase timers. Rust's additional
physical-limit check is included in its total process time.

All measured repetitions pass the analytic check and comparison against an
independent 267-bit, order-112 run. Maximum analytic errors at the base setting
are `2.30e-25` for Rust and `2.19e-31` for the original; at the refined setting they
are `5.35e-35` and `5.10e-41`. The base Rust endpoint error is `2.07e-25`.
Truncation at the matching point dominates Rust's rounding error here.

The internal work differs substantially. Original AMFlow's
[`find_pow_log`](https://gitlab.com/multiloop-pku/amflow/-/blob/26005517a288086c4cb4d1b26d829691bc088485/diffeq_solver/src/mpsolver.cpp#L499)
uses at least 450 decimal digits for indicial reconstruction under these settings.
It [rationalizes the indicial polynomial and obtains its roots](https://gitlab.com/multiloop-pku/amflow/-/blob/26005517a288086c4cb4d1b26d829691bc088485/diffeq_solver/src/frobenius.cpp#L24)
through [an external MPSolve process](https://gitlab.com/multiloop-pku/amflow/-/blob/26005517a288086c4cb4d1b26d829691bc088485/diffeq_solver/include/polyrat.hpp#L375).
Its singular matcher also [enforces at least 20 extra expansion orders](https://gitlab.com/multiloop-pku/amflow/-/blob/26005517a288086c4cb4d1b26d829691bc088485/diffeq_solver/src/mpsolver.cpp#L786).
Rust retains exact exponents 0 and 1 in this example and uses the requested
truncation order. Equal nominal settings therefore give different internal
precision, expansion work and achieved accuracy. These times characterize this
small analytic system; they establish no general speedup or full-workflow result.
The original endpoint is not evaluated in this comparison; Rust's endpoint is
checked against the analytic limit. Both solvers receive supplied boundary data.

[Raw results and provenance](../reports/performance/2026-10-04-singular-resonance.json)
retain each phase timer, numerical output, exact input, executable and library
hash, analytic errors, source references and precision/order checks. The source
inspection annotations were added after measurement; their digest identifies
the unannotated harness output. Reproduce using the original build described above:

```sh
cargo build --locked --release --example benchmark_singular
python3 scripts/benchmark_singular.py \
  --rust target/release/examples/benchmark_singular \
  --upstream target/original-amflow-bench/amflow-26005517a288086c4cb4d1b26d829691bc088485/diffeq_solver/desolver \
  --upstream-build-info target/original-amflow-bench/build.json \
  --rust-source-commit "$(git rev-parse HEAD)" \
  --rust-build-description 'Cargo release profile' \
  --output target/singular-comparison --cpu 24 --repeats 5
```

Choose an available CPU on another host. Each output directory must be new.

## Automatically derived 27-integral paper system

The matrix from the completed four-target sample at `s=30`, `t=-10/3`, `m^2=1`
and `epsilon=1/2700` provides a larger matched ordinary-transport workload.
It has 27 integrals and 190 nonzero entries. Export verifies the cache digest
and exact rational equality of all 729 matrix entries after reparsing; the
[fixture provenance](../fixtures/performance/README.md) preserves its basis
ordering and pole checks. Both implementations receive the same synthetic
integer complex boundary vector and continue from `eta=-512i` to `eta=-256i`.
The computed poles at both working precisions lie inside `|eta|<64`.

These measurements use immutable Rust library commit
`15471e4bbbc7736f1cce62bfb8ff1bad50821dd1`, the original executable above,
CPU 28, one thread, 201 working bits and order 80. The unchanged Rust driver
was linked with `rustc --edition=2024 -O -C debuginfo=1`. Medians include five
fresh processes after one excluded warmup, with alternating implementation order.
The separate Laurent calculation ran concurrently on CPUs 24–27; host load is
recorded. Compilation, matrix export and IBP reduction are excluded.

| Measurement | Rust | Original C++ |
|---|---:|---:|
| Total process | 0.469874 s | 14.481909 s |
| Ordinary transport | 0.262298 s | 0.061754 s |
| Accepted steps | 4 | 1 |
| Rejected trial steps | 3 | Not reported |

Rust's timed preparation takes 0.053344 s. The original total includes its
preparation and external MPSolve calls; its transport remains faster. The two
step-selection and error-checking policies perform different amounts of work.
No overall speedup ratio is inferred from these short supplied-boundary jobs.

All measured repetitions meet the 20-digit scaled-error requirement against
independent 267-bit, order-112 runs. Base cross-implementation error is
`1.45e-21` at most, improving to `3.57e-31` after refinement. Rust's own
precision/order change is at most `5.53e-59`; the original change is `1.45e-21`.
This verifies ordinary transport on the same exact system. The synthetic
boundary is not a physical integral boundary, and the benchmark does not
validate recursive boundary construction, Laurent coefficients or a full
automatic-workflow timing comparison.

[Raw results](../reports/performance/2026-10-04-paper27-de.json) include all
outputs, phase timings, errors, exact inputs and hashes. Reproduce with the
regular benchmark build above and add
`--case fixtures/performance/paper27-eta-eps-1-2700.json` to
`scripts/benchmark_de.py`, choosing a new output directory and an available CPU.
The custom-case loader accepts exact integer/fraction epsilon values and
decimal-string boundary/endpoints; it rejects malformed dimensions and
nonfinite or expression-valued boundary/endpoints without converting through
binary floating-point values.

### Matched nearby W planar/nonplanar characterization

The [recorded comparison](../reports/performance/gg-hg-matched-nearby.json) uses the full 48-master planar and 61-master nonplanar systems, permutation 2, at the same exact W-region source and two nearby destinations. Both implementations requested 20 digits and used the same hashed supplied boundaries and physical root germs. All 3,270 coefficient comparisons, including three repeats and both destinations through epsilon 4, agree within 20 absolute digits. The source evidence remains capped at 24 digits.

Median times in seconds, three sequential repeats on CPU 35 after model/package setup:

| Operation | Planar 48 | Nonplanar 61 |
|---|---:|---:|
| Native one fixed integration | 0.976| 2.046|
| Native checked public first query | 4.533| 11.543|
| Original public first query | 1.060| 2.376|
| Native checked nearby reuse | 4.530| 11.578|
| Original nearby reuse | 1.091| 2.351|
| Native populated exact hit | 0.00190| 0.01746|
| Original populated hit, full driver | 0.0602| 0.2018|
| Original populated hit, inner selection timer | 0.000174| 0.000186|
| Native binary save/load | 0.222/0.117| 0.784/0.541|

Native fixed integration used 60 decimal digits/order 64. The public native call additionally refined to 282 bits/order 96, propagated inherited source uncertainty and populated the cache; it reported 21 verified digits at the first destination and 20 at the second, using the scaled bound `error ≤ 10^(-d) max(1, |coefficient|)`. The independent native/original comparisons require 20 absolute digits. The original used 100 decimal-digit arithmetic and adaptive order. Its inner cold/reuse timer includes boundary selection, integration and saving. The original public driver reloads its matrices; the native model is already prepared. Thus these columns distinguish core integration cost from checked API cost and do not establish an overall speed ratio. Native exact construction/identity/pullback/compile cost 0.610 s planar and 2.225 s nonplanar before repeated measurements; boundary parsing cost under 0.001 s each. Separate source/endpoint insertion medians were 0.989/0.794 ms planar and 3.193/3.059 ms nonplanar.

Both populated nearby queries select the newly computed first endpoint. Both routes already require just one accepted chart, so this short workload cannot reduce the chart count further through nearby reuse. Exact hits perform no integration. The original's supported `$grid=True` setting is required to register saved endpoints for subsequent selection. A historical source-filename mismatch was rejected before integration, and the initial `$grid=False` cache-hit expectation also failed; both harness setup failures are retained in the report. The native measurements were unchanged.

This characterizes canonical master-integral transport, not complete matrix-element runtime. It uses the generic canonical representation, shared numerical engine and persistent RustFlow cache. Future improvements should reuse preparation, validation and persistence work within those shared APIs, preserving the same checks rather than introducing family-specific solvers. Other scientific jobs ran on separate CPUs, so these three observations are a small characterization rather than a controlled hardware benchmark.

### Why fewer rejected steps did not improve the measured segment

A shared step-growth proposal was tested on one saved, expensive segment of the 75-master ZMZ system, using the identical 375 starting coefficients, root sheets and endpoint at 60 working digits/order 56. Both modes used one frozen instrumented binary based on revision `4ff74ea` on CPU 35, with unchanged pole bounds and numerical acceptance checks. The candidate capped each proposal at twice the previous accepted step.

| Controller | Accepted charts | Rejections | Chart construction | Other transport work | Total transport |
| --- | ---: | ---: | ---: | ---: | ---: |
| Existing | 239 | 1261 | 661.84 s | 56.19 s | 718.02 s |
| Growth cap | 252 | 255 | 700.74 s | 30.11 s | 730.85 s |

All 375 endpoint coefficients agreed within 1.25e−48 absolute and the root values were identical. Chart construction accounted for about 92% of the baseline transport time: fewer rejected attempts could not offset the extra charts. This is one run per mode on a shared host, so it establishes neither a statistical slowdown nor a general speed ratio. The candidate was left out of production. The complete scientific state, every coefficient comparison, exact source/binary provenance and timings are preserved in [the controlled-segment report](../reports/performance/pole-step-growth-leg100.json). This supplied-checkpoint experiment is separate from full-path precision refinement.

### Bounded upward step search on the same Taylor chart

A follow-up to the unsuccessful growth cap keeps the same initial proposal and
tries at most two fully checked points between the last failed and first
successful step. One supplied 75-master leg 100 comparison on CPU 35 used a common
binary from validated source `fa1110f`, 60 working digits/order 56 and the same
375-component boundary. All 375 endpoint coefficients agree within 20 absolute
digits (maximum difference `6.718190e-49`), and all three root values are identical.
The prior cap experiment remains separate.

| Controller | Charts | Rejected trials | Superseded successes | Total trials | Chart time | Transport time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Existing proposal | 239 | 1261 | 0 | 1500 | 672.035 s | 729.713 s |
| Two upward trials | 187 | 1179 | 178 | 1544 | 516.157 s | 581.216 s |

Fewer expensive charts outweighed the extra acceptance checks in this one pair.
These are single observations, not a statistical or full-path speed claim.
Every trial uses the same tail, midpoint/endpoint defect and branch checks;
all 1544 candidate trials include the 178 successful results later superseded.
Production adoption remains subject to complete outer-budget accounting and
broader validation. See the [raw outputs, source and provenance](../reports/performance/step-bracket-leg100.json).
