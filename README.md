# symbolica-amflow

A native Rust implementation of auxiliary mass flow, using Symbolica for
exact algebra and MPFR arithmetic and the sibling RustRed checkout for IBP
reduction. Rust 2024; independent MIT-licensed repository.

**Work in progress:** the full AMFlow 2.0 acceptance gate has not passed. See
[coverage and limits](docs/coverage.md). Successful tadpole, bubble, and vacuum
calculations do not establish support for arbitrary multiloop families.

## Build

Requires Rust 1.96 or newer. Symbolica and its companion crates are pinned to
upstream `main` commit
[`75f8350094b90254ee71dc2a391fde0d14b0204a`](https://github.com/symbolica-dev/symbolica/commit/75f8350094b90254ee71dc2a391fde0d14b0204a)
(package version 3.0.1), including the scale-independent root convergence fix.
The root manifest's Cargo patches make RustRed use that same Symbolica revision.

Keep RustRed at `../rustred`, including `crates/rustred-core`. Its experimental
reconstruction feature is disabled. The high-level multiloop solver defaults to
exact rational epsilon specialization before IBP elimination;
`sampled_reduction: false` selects fully symbolic epsilon reduction. This project does not modify RustRed.
Basis refinement also retains symbolic epsilon, taking precedence over sampled
reduction so that dimension-dependent denominator factors can be identified.

```sh
nix develop
# Configure SYMBOLICA_LICENSE in your own shell; never commit a license key.
cargo test --release
cargo run --release --example massless_bubble
cargo run --release --example two_loop_acceptance
```

The last command runs the required four-target calculation and returns an error
if any stage is unsupported or a limit is exhausted. It does not substitute
stored reference values. The fixtures retain upstream precision metadata.

## Interface

- `IntegralFamily`, `Propagator`, `Integral`, and `KinematicPoint` specify exact
  quadratic denominators, numerator slots, signed powers, and substitutions.
- `solve_integrals` prepares the auxiliary differential system, generates
  boundaries, evaluates independent nonzero epsilon samples, and fits a Laurent
  expansion with a second sample grid and increased precision/order.
- `PreparedFlow` and `evaluate_samples` expose preparation and sampled evaluation.
  `BoundaryProvider` permits caller-provided asymptotic data.
- `DifferentialSystem::compile().transport()` propagates a system with supplied
  `BoundaryData`. `frobenius()` constructs generalized power/logarithm solutions.
- `fit_epsilon` fits caller-supplied samples. It reports **no verified digits**
  unless an independent refinement was performed by the orchestrator.
- `ReductionBackend` has native `RustRedBackend`, explicit `TableBackend`, and
  versioned `cache::CachedBackend` implementations. `FlowOptions::cache_directory`
  also persists prepared AMF and FT differential systems, including recursive children and exact basis-transformation histories.

See [library usage](docs/usage.md), the runnable examples, and integration tests
for complete typed inputs.

The integration measure is `d^D l / (i pi^(D/2))` per loop. Denominators are
`q^2 - m^2 + i0`; auxiliary mass is inserted as `D_i - eta`. Consequently
`d I(n)/d eta = sum_i n_i I(n + e_i)` over shifted lines. No Euler-gamma
normalization factor is included. Scalar-product coordinates are the upper
triangle of loop-loop products followed by loop-external products in loop-major
order. Numerator slots must have nonpositive powers.

Defaults are `D=4-2 epsilon`, 20 requested digits, 40 guard digits, AMF recursion,
and one worker. Automatic placement shifts all physical lines at one loop and
one selected physical line at higher loops, preferring a nonzero intrinsic mass. Placement modes also include equal-mass groups, a single propagator, branches,
loops, all lines, and explicit indices. Basis refinement, skipped target reduction,
cancellation, and progress callbacks are available. FT evaluation is restricted to sectors with nonnegative real Euclidean Symanzik coefficients; its recursive Gaussian terminals support tensor numerators up to rank eight.

Arbitrary-precision values remain Symbolica `Complex<Float>` values. MPFR
operations round at an explicit working precision; working precision is not an
accuracy guarantee. Verified coefficient accuracy uses
`|new-old| <= 10^-digits max(|new|, |old|)` for resolved nonzero coefficients.
Coefficients smaller than `10^-digits` use an absolute tolerance of
`10^-digits`; relative accuracy is not claimed for numerical zeros.

## References

Target: [AMFlow 2.0](https://gitlab.com/multiloop-pku/amflow/-/tree/26005517a288086c4cb4d1b26d829691bc088485),
commit `26005517a288086c4cb4d1b26d829691bc088485`.

- X. Liu and Y.-Q. Ma, [AMFlow: A Mathematica package for Feynman integrals
  computation via auxiliary mass flow](https://arxiv.org/pdf/2201.11669).
- R.-J. Huang, X. Liu and Y.-Q. Ma,
  [AMFlow 2.0](https://arxiv.org/html/2607.08477v1).

Upstream attribution and license terms are preserved in `NOTICE` and
`fixtures/amflow-2.0/LICENSE.md`. Symbolica's separate licensing terms apply to
that dependency.

The complete epsilon-dependent indicial sector is retained during numerical
sampling. For exact sampled IBP systems it is reconstructed uniquely from a
half-integer base and an integer epsilon slope bounded by twice the loop count;
ambiguous samples are rejected. Epsilon values introducing additional indicial
resonances are also rejected by automatic AMF evaluation. Small automatically
generated samples avoid these cases. Standalone systems support integer
resonances and logarithmic solutions.
