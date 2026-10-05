# Epsilon poles in physical differential matrices

`RustFlow::regularize_epsilon` admits rational physical systems whose epsilon
poles can be removed by a common diagonal rescaling. It reuses the ordinary
physical transport and progressive cache. This extends the direct coefficient
recurrence, which requires matrices regular at epsilon zero; it introduces no
new numerical solver or finite-epsilon approximation.

Write the original master vector as `I = T J`, where
`T = diag(epsilon^w_i)`. For every nonzero entry of every physical partial
matrix, the exact integer weights obey

```text
valuation_epsilon(A_s[i,j]) + w_j - w_i >= 0.
```

The maximum weight is zero. The same transformation is applied to all physical
coordinates before choosing a path; cancellation between different partials on
one path cannot remove an original condition. Symbolica supplies Gaussian
rational polynomial valuations and matrix transformations. The existing
Frobenius integer constraint algorithm is shared, with Linnet graph storage.
The physical derivative of this epsilon-only transformation is zero. Physical
Fuchsian normalization still includes its own nonzero derivative term.

A negative constraint cycle means this diagonal class is insufficient. It does
not prove that a nondiagonal rational change is impossible. Algebraic coefficient
fields, general epsilon regularization and sampled-epsilon fallback remain
separate work. In particular, a supplied equation such as `I' = I/epsilon` has
an essential epsilon singularity; numerical fit stability alone would not
establish a Laurent expansion.

## Boundary orders and errors

The exact mapping is `J_i[n] = I_i[n + w_i]`. Restoring every original component
through order `H` requires a rectangular transformed range ending at
`H - min(w)`. Coefficient errors move with the same entries, without rescaling,
rounding through binary64, or increasing the evidence cap.

The declared original leading bound must hold throughout the intended transport,
not merely describe the first nonzero coefficient at the starting point. For
example, `I_1' = I_2/epsilon` with starting values `(0,1)` develops an epsilon
pole. Ordinary epsilon-independent integral families retain their established
common `-2*loops` bound. Caller-supplied systems need their own valid bound.

Missing positive source orders are errors, never zeros. With weights `[-1,0]`,
restoring through epsilon zero needs `J` through epsilon one. The second original
source component must then also be known through epsilon one. Even when its
finite coefficient is known exactly, the missing next coefficient can change a
finite answer after transport.

## Supplied systems and the growing cache

The returned `EpsilonShearedFlow` retains the original identity and exposes the
transformed connection through `flow()`. Its boundary conversion methods verify
identity, dimensions, source orders and accuracy before producing a new entry:

```rust,ignore
let regular = original_flow.regularize_epsilon(&context)?;
let original_range = EpsilonRange::new(-2, 0)?;
let cache_range = regular.required_range(original_range)?;
let seed = regular.to_sheared_boundary(&original_boundary, cache_range.last)?;
cache.insert(seed)?;
let result = regular.flow().evaluate_to(
    &mut cache, &destination, cache_range, &options, &context, &policy,
)?;
let original_values = regular.to_original_boundary(&result.boundary, original_range)?;
```

The input boundary must already cover the required original source orders.
Exact scaled basis labels distinguish the rescaled values even if another
connection has identical matrices. Original denominator holes, generic epsilon
conditions, normalization, physical prescription and homotopy declarations remain
part of the identity. A prescribed flow still requires the prescribed transport
entrypoint. Original and transformed boundaries can coexist in one bank and
survive its normal binary restart; no cache schema change is needed.

## Reduced integral families

`PreparedPhysicalFamily::with_epsilon_shearing` enables the same transformation
for a reduced family. `basis()` and `target_reductions()` continue to describe
ordinary integrals. `epsilon_shearing()` records their relation to the columns
exposed by `flow()`.

Automatic `seed_cache` evaluates the original integrals with AMF or FT, including
all requested source orders, and only then converts their coefficients and
errors. `required_master_range` and `project_targets` both multiply the retained
target weights by the same exact `T`; epsilon poles in those weights therefore
request extra cached orders before projection. Basis refinement and skipped
initial reduction remain independent preparation choices.
