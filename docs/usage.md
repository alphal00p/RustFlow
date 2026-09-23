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
before the next search round. Set `checkpoints` on
the backend to resume unfinished reductions after increasing a target budget;
the acceptance examples enable this alongside their ordinary reduction cache.

Use `RunContext::progress` for progress callbacks and its shared
`CancellationToken` for cooperative cancellation. Native IBP calls check
cancellation at their boundaries; the current RustRed adapter cannot interrupt
an in-progress native reduction.

See `docs/coverage.md` for the current completion status and explicit limits.
