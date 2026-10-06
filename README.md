# RustFlow (`symbolica-amflow`)

A native Rust implementation of auxiliary mass flow, using Symbolica for
exact algebra and MPFR arithmetic and RustRed for IBP
reduction. Rust 2024; independent MIT-licensed repository.

RustFlow is being extended with DiffExp-style transport in physical kinematic
variables and a progressively filled native `RustFlowCache`. The Rust crate
retains the name `symbolica-amflow`. See the [active parity plan](docs/parity-plan.md)
and [cache and transport interface](docs/rustflow.md). Full AMFlow/DiffExp parity
remains a development goal. The [Python API and native Higgs+jet notebook](docs/python-notebook-status.md)
generate all sixteen physical boundary configurations from empty numerical
caches and reproduce 4,360 transport coefficients, eight W/Z form factors and
three coherent observables. Published numerical seeds are comparison data only.
The independent reference supports 19 relative comparison digits for the EW
square and 20 for the infinite-top HEFT–EW interference; native uncertainty
estimates and higher-precision checks are recorded separately. General automatic
amplitude evaluation remains open.
The [publication-wheel acceptance](reports/validation/2026-10-06-gg-hg-publication-complete/summary.json)
also passed independent 40-digit regeneration, interruption/resume, binary
restart and nearby-point reuse, with final runtime attestation and process exit zero.

The notebook is `examples/hep/gg_hg.py` in
[symbolica-community PR #19](https://github.com/symbolica-dev/symbolica-community/pull/19),
with staged computation, native diagram displays and reusable binary caches.
See the [checkout and installation instructions](docs/clean-community-build.md)
and [current validation status](docs/python-notebook-status.md) for the exact
published dependency snapshot and remaining generality, packaging and CI work.

The original paper's four two-loop targets pass the required calculation through
epsilon power zero: 20-digit stability under independent sample/precision/order
refinement and agreement with upstream data at its recorded precision. See the
[acceptance report](reports/validation/2026-10-04-paper-two-loop-acceptance.json)
and [implemented coverage and limits](docs/coverage.md).
[Whole-segment Taylor checks](docs/transport-verification.md) now address the
[recorded sparse-polynomial verification failure](reports/validation/2026-10-05-residual-alias-counterexample.json).
The [release gate](reports/validation/2026-10-05-conditioned-transport-release.json) passes 469 tests and one doctest, with 10 long scientific regressions opt-in.
[Endpoint conditioning and bounded precision retries](docs/arithmetic-conditioning.md)
address the separate [large-cancellation counterexample](reports/validation/2026-10-05-coordinate-roundoff-counterexample.json).
Ordinary rational charts also check directed defect enclosures against exact
source coefficients. Optional [regular Möbius charts](docs/local-coordinates.md)
share the same endpoint, path and cache checks. These local checks and independent
refinements are not a general accumulated-error proof.
The [implementation and benchmark comparison](docs/benchmark-comparison.md)
lists completed and pending features with RustFlow and reference timings side by side.
An executable [performance comparison](docs/performance.md) measures the Rust
and original AMFlow 2.0 differential-equation solvers on identical inputs.
The [Bern Wolfram setup](docs/wolfram-setup.md) also enables live checks against
the original Mathematica package with Kira and the university license server.

## Build

Requires Rust 1.96 or newer. The standalone lock pins Symbolica and its companion
crates to community commit
[`c3408e4ba1d3bdd4ea55678fad50e27009be13d4`](https://github.com/symbolica-dev/symbolica/commit/c3408e4ba1d3bdd4ea55678fad50e27009be13d4)
(package version 3.0.1), including the root convergence and exact division fixes.
The shared Python host uses this same numerical and native type graph; see
[embedding and source fingerprints](docs/dependency-embedding.md).

RustRed is a Git dependency, locked to `7c1ed037`; its experimental
reconstruction feature is disabled. The high-level multiloop solver defaults to
exact rational epsilon specialization before IBP elimination;
`sampled_reduction: false` selects fully symbolic epsilon reduction. Upstream
RustRed updates are checked between validation milestones; local work is
preserved. Both reduction adapters now use RustRed’s compiled arity registry
(default 1–16); the factorized path delegates to its native dispatch macro.
`RUSTRED_RUNTIME_ARITIES` configures that registry at build time. Search and
representation limits still apply.
Basis refinement also retains symbolic epsilon, taking precedence over sampled
reduction so that dimension-dependent denominator factors can be identified.

The manifest fetches the public HEPKit/Linnet/Spenso/Idenso owner at
[`a3d1c8a867ac0e89d8f58f9722cf9a7e83338901`](https://github.com/ValentinHirschi/gammaloop/commit/a3d1c8a867ac0e89d8f58f9722cf9a7e83338901),
including its native rendering components. Normal builds need no sibling checkout
or local patch configuration. The community host keeps Vakint at its separate
[`6203c6cbba6ae5e90329ba5081fad55319e678db`](https://github.com/ValentinHirschi/gammaloop/commit/6203c6cbba6ae5e90329ba5081fad55319e678db)
revision and shares the same RustRed and algebra sources. See the
[published community build recipe](docs/clean-community-build.md) and
[embedding requirements](docs/dependency-embedding.md); the community manifest
and lock select its RustFlow runtime revision. The Cargo config example is only
for optional standalone development against a local native owner.

Python features are optional and disabled by default. The bindings are
registered by the community host under `symbolica.community.hep.integration`.
Actual dependency source content is included in persistent-cache identities.
Public source availability and component checks are separate from a completed
clean-machine build and full native notebook acceptance for that graph.
The [native Higgs+jet boundary runner](docs/gg-hg-native-boundaries.md) exercises
the same boundary generation and refinement without Python. Its long acceptance
runs keep completed samples separate from verified boundary values.

```sh
nix develop
# Configure SYMBOLICA_LICENSE in your own shell; never commit a license key.
cargo test --release
cargo run --release --example massless_bubble
cargo run --release --bin rustflow -- graph examples/cli/massless-bubble.json
cargo run --release --bin rustflow -- transport examples/cli/physical-transport.json
cargo run --release --bin rustflow -- transport examples/cli/prescribed-transport.json
cargo run --release --example two_loop_acceptance -- \
  --symmetry-rules --parametric-rules --depth 3 --case-batch 32 \
  --max-targets 262144 --max-exact-frontier 512 \
  --max-backward-frontier 4096 --workers 4
```

The last command runs the required four-target calculation with the validated
reduction settings; the default backend budget is intended for smaller families.
It returns an error
if any stage is unsupported or a limit is exhausted. It does not substitute
stored reference values. The fixtures retain upstream precision metadata.
The [CLI steering format](docs/cli.md) accepts native HEPKit DOT/model files and
supports persistent transport caches across processes.

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
- `RustFlow::evaluate_to(&mut cache, ...)` searches compatible physical
  boundaries, transports epsilon coefficients, and adds verified intermediate
  points and the destination to `RustFlowCache`. Binary snapshots can be loaded
  before the next request. The ordinary interface uses regular straight paths.
  [`with_prescribed_continuation` and `evaluate_prescribed_to`](docs/prescribed-cache.md)
  also support threshold detours with exact local prescriptions, explicit endpoint
  root germs, and caller admission of the global continuation domain. [`RustFlow::new_algebraic`](docs/algebraic-cache.md)
  uses the same bank with explicit root germs and retains their sheets across
  regular rational-complex affine paths and binary restart.
  [Accuracy-aware retries](docs/accuracy-fallback.md) try another admissible source
  when the first cannot meet accuracy, preserving failed-attempt diagnostics.
- `PreparedPhysicalFamily` derives a common physical differential basis from a
  family and appends independently verified AMF/FT seed values to the same bank.
  Original reduction assumptions and physical chart restrictions remain part of
  cache identity and are checked along candidate paths. Retained exact target
  reductions determine required master orders and project cached coefficients
  with propagated uncertainty.
- [`RustFlow::regularize_epsilon`](docs/epsilon-shearing.md) removes epsilon poles
  in rational and registered-square-root systems when one exact diagonal
  rescaling works for every physical partial. Root relations are normalized
  before epsilon expansion, with excluded denominators and sheets retained. Its cache
  adapter preserves basis identity, source orders and errors; the corresponding
  `PreparedPhysicalFamily::with_epsilon_shearing` builder maps automatic AMF/FT
  seeds and target reductions through the same transformation.
- [`algebraic::AlgebraicSystem`](docs/algebraic.md) transports supplied systems
  containing registered square roots on explicit regular contours. Each call
  tracks its own sheets through accepted steps and contour windings. Multivariate
  dlog construction and exact pullback preserve root relations and source domains.
  The [five-point regression](docs/fivepoint-planar.md) checks all 13 one-loop
  masters through epsilon power four against original DiffExp.
- [`CanonicalAlgebraicSystem`](docs/canonical-systems.md) retains canonical dlog
  letters until a path is supplied. `RustFlow::new_canonical` reuses those exact
  data across queries in the same physical boundary cache. The
  [13-master cache validation](docs/full13-canonical-cache.md) checks nearby reuse
  and restart against original DiffExp, preserving supplied-boundary accuracy.
- [`contour::PrescribedContour`](docs/prescribed-contours.md) constructs real-axis
  detours from explicit polynomial prescriptions, retaining exact root identities
  and certified isolation disks from Symbolica.
- `hepkit::GraphIntegral` retains a native HEPKit graph and delegates routing,
  contraction, denominator completion and scalar numerator rewriting to its
  Linnet/Spenso interfaces. `solve_integral_combinations` evaluates the resulting
  exact weighted integral groups before Laurent fitting. Its cut interface
  preserves native `DiagramCut` orientation and denominator slots.
- `FrobeniusBasis::match_constraints` fixes integration constants from partial
  asymptotic coefficients, including resonant powers and logarithms. The full
  [equal-mass banana regression](docs/banana-equal.md) uses an analytic infinity
  boundary and checks every coefficient through epsilon power four. The
  [unequal-mass regression](docs/banana-unequal.md) checks all 75 coefficients
  of 15 masters with independent routes and a live original DiffExp reference.
  The [108-master nonplanar system](docs/diffexp-nonplanar108.md) uses
  `AlgebraicSystem::analytic_origin` for a compatible analytic sector at a singular
  source and validates all 540 coefficients with independent refinement.
- `cuts::CutFamily` preserves distributional and positive-energy cut metadata
  through native IBP reduction. The [phase-space terminals](docs/cuts.md) cover
  unequal-mass two-body and massless N-body final states, using native HEPKit
  normalization, exact routing Jacobians and dimensional gamma factors.

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
accuracy guarantee. Automatic Laurent reconstruction checks coefficient accuracy
using `|new-old| <= 10^-digits max(|new|, |old|)` for resolved nonzero coefficients.
Coefficients smaller than `10^-digits` use an absolute tolerance of
`10^-digits`; relative accuracy is not claimed for numerical zeros. The physical
transport cache uses `error <= 10^-digits max(1, |coefficient|)` and returns component error
estimates; its reported digits use an absolute criterion below unit magnitude.

## References

Target: [AMFlow 2.0](https://gitlab.com/multiloop-pku/amflow/-/tree/26005517a288086c4cb4d1b26d829691bc088485),
commit `26005517a288086c4cb4d1b26d829691bc088485`.

- X. Liu and Y.-Q. Ma, [AMFlow: A Mathematica package for Feynman integrals
  computation via auxiliary mass flow](https://arxiv.org/pdf/2201.11669).
- R.-J. Huang, X. Liu and Y.-Q. Ma,
  [AMFlow 2.0](https://arxiv.org/html/2607.08477v1).
- M. Hidding, [DiffExp](https://arxiv.org/abs/2006.05510), pinned to
  `784c8229bf92369a03f011a48e161522c8c54bbd`.

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
