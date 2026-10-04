# Comparison with original AMFlow 2.0

A same-host comparison against the original, unmodified C++ differential-equation
solver is available. It measures preparation and regular continuation with supplied
boundaries. **There is no measured full automatic-workflow comparison yet.** No
Wolfram runtime was found in this environment, and the mandatory four-target
two-loop acceptance calculation remains incomplete.

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
(`diffeq_solver/src/mpsolver.cpp`, `regular_expansion`). This identifies a concrete
optimization target. The different step counts and error checks account for
additional differences; no single cause is inferred from total timings alone.

## Full automatic workflow

[`scripts/benchmark_full.wl`](../scripts/benchmark_full.wl) is an optional harness
for a configured Wolfram installation and upstream-supported IBP reducer. It
provides the physical massless bubble, two-mass sunset, and all four required
paper targets through epsilon power zero at 20 requested digits, one thread,
and caching disabled. It checks the pinned `AMFlow.m` digest and writes results
and timing metadata in the current directory. Use a fresh directory for each run.

This harness has **not been executed** here. Its output explicitly marks accuracy
unverified; results must be compared and independently refined before timing
claims are made. Rust automatic evaluation already performs additional sample
and precision refinement, while upstream `SolveIntegrals` uses its own numerical
configuration. A future full comparison must report this difference and include
equivalent correctness checks. Stored upstream answer files contain no timing
evidence and are not used to manufacture an original runtime.
