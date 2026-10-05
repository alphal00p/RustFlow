# Regular local coordinates

`FlowOptions::local_coordinate` selects `LocalCoordinate::Identity` (the default)
or `LocalCoordinate::BalancedMobius`. Both use the same Taylor transport,
refinement, physical error propagation and progressive point bank. The opt-in
choice also applies to direct epsilon systems and explicitly registered roots.
The CLI exposes the same selection as `options.local_coordinate`, accepting
`"identity"` or `"balanced_mobius"`.

```rust
use symbolica_amflow::{FlowOptions, LocalCoordinate};
let options = FlowOptions {
    local_coordinate: LocalCoordinate::BalancedMobius,
    ..Default::default()
};
```

Each chart starts at its known physical boundary. On a straight leg from `c` to
`t`, write `x = c + (t-c) a y/(1+b y)`. The nearest real poles choose left/right map endpoints. Nonreal poles
contribute symmetric radial distances about the expansion center, avoiding an
artificial stopping barrier at their real projections. The balanced map takes
these chosen endpoints to `-1` and `1`; a missing side uses the exact one-sided limit. With no poles the ordinary coordinate is retained. Rounded
working coordinates and poles choose the chart, and their binary rationals are
retained explicitly; this does not add precision to the physical inputs.

The exact original connection is substituted into the map and multiplied by its
Jacobian. Original denominator/radicand/domain conditions are substituted
independently before cancellation, and the map denominator is always retained.
The actual transformed complex poles bound a proposed local step. Real
projections alone do not certify convergence or path safety. Exact source
recompilation is deliberate in this first implementation; no speedup is claimed.

Trial endpoints and midpoints remain on the original straight physical leg.
They are mapped back to `y` before evaluating the Taylor polynomial. Derivatives
are divided by `dx/dy` before checking the original physical differential
equation. Taylor tails, endpoint and midpoint defects, branch square/root
residuals, step limits and cancellation retain their existing acceptance roles.
Every trial also requires the shared whole-segment residual estimate built from
the full pulled-back polynomial rows; a sparse zero tail or two vanishing sampled
defects cannot authorize it. Registered-root estimates use the mapped logarithmic
derivative and rational-coefficient owner, while physical branch projection uses
the original radicands. Ordinary rational and epsilon charts form directed
residual enclosures against the exact mapped source. Registered-root assembly
still uses MPFR majorant estimates. Neither establishes a separate global
forward-error proof. Endpoint conditioning and bounded precision retries are
shared with identity charts.
A chart coordinate never selects a new contour or authorizes a singular crossing.

A registered root chart starts with the signed roots of the accepted physical
state. It does not reseed principal roots. Rejected steps do not update branch
state. `TaylorSegment.center` and `.end` remain the traversed physical interval;
its `coordinate` records the polynomial argument. Saved-point evaluation and
independent refinement use that metadata. The point bank still stores physical
coordinates, root germs, route identity and inherited accuracy evidence; its
binary schema does not change because local Taylor segments are not serialized
into the bank.

This implements regular boundary-centered Möbius charts. It does **not** implement
predivision with an expansion center ahead of the known boundary, off-center
fundamental-matrix matching, or singular-center Frobenius charts. Those require
separate matching/conditioning and branch evidence. The selected coordinate is a
solver option, not a new physical basis or branch identity. Both optional
`StepSizeStrategy` policies use the same mapped acceptance predicate and count
every trial against the shared predicate budget. The regular complex-root,
coupled-epsilon, saved-point, cache-restart and refinement regressions cover this
implementation; they do not establish full DiffExp predivision parity or a
scientific benchmark speedup.
