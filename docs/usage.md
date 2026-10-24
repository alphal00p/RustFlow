# Library usage

The integration measure is `d^D l/(i pi^(D/2))` per loop, without an
`exp(epsilon*EulerGamma)` prefactor. All family coefficients must be exact.

```rust
use symbolica::prelude::*;
use symbolica_amflow::*;

fn tadpole() -> Result<Vec<LaurentExpansion>> {
    let family = IntegralFamily {
        name: "massive_tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator::quadratic(&[1], &[], Atom::num(1), &[])?],
        physical_propagators: 1,
        epsilon: symbol!("eps"),
        dimension: 4,
    };
    solve_integrals(
        &family, &[Integral(vec![1])], &KinematicPoint::default(), 0,
        &FlowOptions::default(), &RustRedBackend::default(), &RunContext::default(),
    )
}
```

`LaurentExpansion::coefficients` is keyed by the epsilon power. Coefficients are
Symbolica `Complex<Float>` values. `verified_digits` records stability under an
independent calculation; `working_bits` records arithmetic precision.
`comparison_errors`, `validation_samples`, and `refinements` expose the evidence
used for that stability estimate. This is a numerical estimate, not an interval
arithmetic proof. The accuracy convention is relative for resolved nonzero coefficients. When
both fits have magnitude below `10^-digits`, only absolute accuracy is claimed.

For repeated evaluations, construct `PreparedFlow::new` and use
`evaluate_samples` with exact rational, nonzero epsilon values and a
`recursive::RecursiveBoundary`. A system prepared using `new_at_epsilon` belongs
to that one epsilon value. `solve_prepared` accepts a custom `BoundaryProvider`.
A `TerminalProvider` can augment recursive terminal integrals through
`RecursiveBoundary::with_terminal`.

For unrelated linear ODEs, construct `DifferentialSystem { variable, matrix }`,
compile at a `Precision` and parameter map, then call `transport` with supplied
`BoundaryData` and a path. `plan_path` constructs pole-avoiding paths with a
specified side of real singularities. `frobenius` constructs singular endpoint
and infinity expansions; `invert_variable` maps infinity to zero.

For several parameter samples or precision refinements, call
`system.prepare_frobenius(&context)` once, then
`prepared.evaluate(precision, &parameters, order, &context)`. The immutable
`frobenius::PreparedFrobenius` retains exact normalization, indicial exponents,
generalized eigenspaces, and leading logarithmic chains. Each evaluation owns
its numerical recurrence workspace, so an `Arc` can share the preparation
between workers. Original denominator restrictions remain visible through
`nonzero_conditions()`: a singular expansion center is allowed, but parameter
samples that annihilate an entire restriction or change the indicial resonance
pattern are rejected. Such nongeneric parameters require a separately prepared,
exactly specialized system.

`PreparedFlow` lazily shares these exact preparations at zero and infinity
across its epsilon samples and refinement attempts. Public changes to its exact
ODE matrix or variable invalidate the cached preparation; supplied declarations
also retain their existing mutation seal. Cancelled or failed preparation does
not install an entry and can be retried with a fresh context. Cancellation is
checked between exact owner operations, matrix rows, eigenspaces, and numerical
recurrence steps; a single Symbolica matrix operation remains indivisible.
Progress callbacks may request cancellation but must not reenter the same
calculation while its exact preparation is in progress.

`ReductionBackend` supports native RustRed and explicit tables. A table must
cover requested derivatives as well as initial targets. Search residuals are
candidates in a spanning basis, not certified masters: the port requires stable
exact differential closure and consistent asymptotic constraints.

Set `cache_directory` to persist prepared systems. Wrap the reducer in
`cache::CachedBackend` to persist intermediate reductions too. Entries include
source digests, exact inputs, backend settings and relevant flow options.
An incompatible or corrupt entry is an error. Cache writes are atomic.

Numerical boundary memoization identifies equivalent denominator permutations,
unused numerator slots, and family display-name changes. Symbolic system caches
retain the caller’s integral ordering.

`MassMode::{Mass, Propagator, Branch, Loop, All, Propagators}` selects the
auxiliary-mass placement. `Mass` groups equal nonzero intrinsic masses; `Branch`
chooses parallel loop-momentum directions; `Loop` uses independent branch spans,
equivalent to completing monomials of the first Symanzik polynomial. The default
`Auto` shifts all physical lines at one loop and a preferred massive line at
higher loops. Explicit modes return an error when no eligible placement exists.

`RustRedBackend::factorized` (enabled by default) selects RustRed’s exact factorized rational-polynomial
field for shared finite-corner replay. It uses the same native sources, modular
search and exact replay, and does not enable experimental reconstruction.
This path emits per-sector reduction progress, completes higher sectors before
searching lower ones, and removes exactly canceled intermediate contributions
before the next search round. For large frontiers, `max_exact_frontier` (512 by
default) instead permits extra dependency searches and postpones coefficient
substitution. Setting it to zero disables intermediate coefficient pruning;
final reductions remain exact, and unsearched nonzero contributions are errors.
Set `checkpoints` on
the backend to resume unfinished reductions after increasing a target budget;
the acceptance examples enable this alongside their ordinary reduction cache.

Use `RunContext::progress` for progress callbacks and its shared
`CancellationToken` for cooperative cancellation. Native IBP calls check
cancellation at their boundaries; the current RustRed adapter cannot interrupt
an in-progress native reduction.

See `docs/coverage.md` for the current completion status and explicit limits.

## Supplied auxiliary-mass equations

`PreparedFlow::from_supplied` accepts an `engine::SuppliedAuxiliarySystem` with
an ordered ordinary-integral basis, its closed differential matrix and target
maps, nonzero conditions, an explicit propagator deformation mask, and source
provenance. The family is specialized once with the given exact kinematic
point; the supplied matrix and maps must already belong to that point. Their
coefficients must be exact rational functions of the declared auxiliary
variable and epsilon. The declared dimension must agree with the evaluation
options.

This interface checks structure and domains, not the mathematical truth of
caller-supplied IBP identities. Unresolved map leaves are rejected. Automatic
recursive boundary generation uses the declared deformation throughout, even
when the evaluation options would otherwise choose another mass placement.
Target maps are applied before extracting the physical endpoint limit. Retained
nonzero conditions restrict the continuation path even if they canceled out of
the differential matrix; a condition annihilated by an exact epsilon sample is
an error.

Prepared supplied declarations are sealed against later mutation, and completed
sample cache keys include their matrix, basis, maps, mask, conditions, and source.
Changed declarations require preparation of a new system. The specialized
`HiggsJetIntegralSystem` exposes its associated mathematical inputs through the
read-only `basis_map()`, `canonical_system()`, and `transport()` accessors so its
boundary identity cannot be detached from those inputs. Its verified boundary
identity also includes the extracted map contents, independently of upstream
file provenance.

## Long acceptance runs

Use `two_loop_acceptance --sample 1/2700` to evaluate all four targets at one
exact nonzero epsilon, with an independent increase in working precision and
series order. This mode prints the change for each target. Omit `--sample` for
the complete Laurent reconstruction, additional sample grid, and upstream
reference comparison; a successful individual sample does not satisfy that
acceptance test.

The recorded four-target sample used the following settings (about six minutes
from an empty cache on the recorded host):

```sh
cargo run --release --example two_loop_acceptance -- \
  --sample 1/2700 --symmetry-rules --parametric-rules \
  --depth 3 --case-batch 32 --max-targets 262144 \
  --max-exact-frontier 512 --max-backward-frontier 4096
```

The default reduction budget is intended for smaller families. The
[sample report](../reports/validation/2026-10-04-paper-four-samples.json)
records the exact executable, input, precision changes, and timings. Keep the
same reduction settings and omit `--sample 1/2700` to run the complete Laurent
acceptance calculation. `--workers N` bounds simultaneous epsilon evaluations.
The [complete acceptance report](../reports/validation/2026-10-04-paper-two-loop-acceptance.json)
records the successful 27- and 31-point fits at 60 and 80 working digits, with
outer orders 80 and 96. All four Laurent expansions pass 20-digit stability and
the upstream recorded-precision comparison. Its elapsed time includes a warm
initial cache and additional exact-preparation workers; it is not a cold timing
for the command above.

The `two_loop_acceptance` example accepts `--cancel-file PATH`. Creating that
file requests cooperative cancellation and saves completed native search work.
Remove the file before resuming the same command. `--native-workers N` controls
concurrent native batches per epsilon sample; budget resources for the product
of this value and `--workers N`.
`--checkpoint-seconds N` sets the periodic save interval (60 by default).
Increasing it can help when individual checkpoints are several gigabytes;
cancellation and search completion still force a save. The `paper_subsector`
diagnostic also accepts `--depth N`, `--max-targets N`, and
`--max-exact-frontier N`.

## External Higgs-plus-jet mathematics

Published Higgs-plus-jet inputs are fetched separately from the library. Before
constructing the published system or form-factor objects in a notebook, use:

```python
from symbolica.community.hep.integration import load_higgs_jet_data
await load_higgs_jet_data(form_factors=True)
```

Native Python and browser/WASM use the same asynchronous API. Downloads use an
immutable Git revision; the library verifies BLAKE3 hashes before accepting any
payload. Validated files are cached by hash. With `form_factors=False`, the large
coefficient document is not downloaded. General reduction and transport do not
need these example inputs.

Rust callers can provide documents with `gg_hg::data::install`, or explicitly set
`RUSTFLOW_HIGGS_JET_DATA_DIR` to a directory containing the matching JSON files.
Mathematical fingerprints retain their existing values; authenticated inputs
protect boundary-cache identities. Tests read repository fixtures; production
builds do not embed them.
