# RustFlow command-line interface

Build with `cargo build --release --bin rustflow`. The executable accepts one
JSON steering file and writes JSON results to standard output. All paths in the
steering file are relative to that file. Unknown fields and incompatible schema
versions are rejected. Progress, when enabled with `options.progress`, goes to
standard error.

```bash
cargo run --release --bin rustflow -- graph examples/cli/massless-bubble.json
cargo run --release --bin rustflow -- transport examples/cli/physical-transport.json
cargo run --release --bin rustflow -- transport examples/cli/algebraic-transport.json
```

## Native HEPKit graph input

The graph example references a native HEPKit model JSON and an unmodified native
DOT graph. HEPKit reads, validates and routes the graph through Linnet. Idenso and
Spenso contract its numerator; HEPKit completes the denominator basis and rewrites
the scalar numerator. RustFlow evaluates the resulting weighted integrals,
applying all numerator weights before Laurent fitting.

`scalar_products` supplies exact native kinematic assumptions as `{left, right,
value}` entries. Momentum arguments are native labels such as
`gammalooprs::P(1)`. `substitutions` supplies exact scalar parameter values; model
floating-point defaults are never used as integral inputs. Expressions are strings
and can contain exact rational-complex coefficients. Use `"1/3"` rather than a
JSON floating-point number.

Unqualified scalar symbols use the native DOT namespace `feynkit_graph` by
default; `namespace` overrides this. Native model parameters use the explicit
`UFO::` namespace. The default regulator is `eps` in the selected namespace.
`tensor_dimension` names the symbolic tensor dimension (default `D` in the
selected namespace). Use that symbol in native numerator slots, for example
`spenso::mink(D,mu)`. After contraction the CLI substitutes
`D = options.dimension - 2 eps`, with default `options.dimension = 4`.
Surviving tensor slots with another Lorentz dimension are rejected. A trace
explicitly written in four dimensions may already be the scalar `4` after
native parsing and remains that scalar; RustFlow cannot recover a symbolic
dimension from an already contracted expression.

`edge_powers` maps native internal edge IDs to signed powers. It is separate from
DOT graph metadata. `last_epsilon_power` defaults to zero. The result includes
the Laurent coefficients, independent comparison changes, verified digits,
working bits, and sample counts. This evaluates the graph's supplied scalar
numerator; it does not build a full scattering amplitude from Feynman rules.

## Progressive physical transport

The transport example supplies `derivatives`, mapping each kinematic variable to
its exact differential matrix, an ordered `basis`, its `normalization`, and a
`branch_domain` identifying the common physical sheet of its inputs. Matrices
must be regular in epsilon. Without a root registry they have rational dependence
on epsilon and kinematics. Optional registered square roots use the algebraic
native API described below. This command follows regular straight paths and
rejects singularities on the selected path.

`leading_epsilon_power` and `last_epsilon_power` define the common coefficient
range. Coefficient arrays are epsilon-order-major, then basis-component-major.
Each initial `seed` supplies exact coordinates, decimal-string complex values,
absolute error estimates, independently justified `verified_digits`, storage
`working_digits`, and provenance. Working digits alone do not establish accuracy.
The supplied uncertainty is propagated and combined with independent transport
checks.

`cache_directory` contains the binary `RustFlowCache` snapshot. An existing
nonempty bank is loaded before seeds are inserted. Each entry in `destinations`
is evaluated against the progressively enlarged bank. Independently checked
intermediate physical points and the requested destination are retained, and
the snapshot is saved after every successful destination. A failed later request
therefore preserves earlier progress. Repeating an exact compatible destination
uses zero transport steps. On subsequent runs `seeds` may be omitted.

The default search compares arbitrary-precision Euclidean distances, resolving
rational distance ties exactly. Optional exact positive `distance_scales` change
the relative importance of coordinates. The branch-domain assertion must be
valid for the supplied seeds and requested straight paths; more specialized
admissibility and cost policies are available through the native library API.
Cache format, source and dependency fingerprints prevent incompatible reuse.

Both commands accept `options` for requested digits, guard digits, series order,
worker count, step/refinement limits, AMF/FT recursion, `+i0`/`-i0`, basis
refinement, initial-reduction skipping and symbolic caching. The default values
match `FlowOptions`. Persistent numerical snapshots and symbolic reduction caches
serve different purposes and use separate paths.

The Python API will live in HEPKit; this executable uses only the native Rust
libraries.

Transport requests may include `nonzero_conditions`, a list of exact expressions
carried from IBP reduction or a kinematic chart. These restrictions participate
in cache identity and are checked at stored points and along candidate paths,
including restrictions that cancel out of the final differential equations.
Epsilon remains formal: `epsilon*s` excludes `s=0`, while a pure epsilon factor
does not exclude a physical point. Use the request's actual epsilon symbol.

## Registered square roots

An optional `roots` map registers exact radicands by symbol, for example
`"roots": {"r": "s"}` for `r²=s`. Derivative matrices use the declared symbol,
such as `"derivatives": {"s": [["eps/r"]]}`; do not put an unregistered `sqrt`
function into a matrix. Each algebraic seed supplies an explicit `root_germ`,
and each destination uses the following structured form:

```json
{
  "coordinates": {"s": "4"},
  "root_germ": {"r": "principal"}
}
```

The allowed sheet names are `principal` and `opposite`. Every registered root
must occur exactly once, using the same namespace rules as the matrices. The
principal root of a negative real radicand is positive imaginary; the opposite
sheet reverses the sign. These are local discrete sheet choices, not approximate
root magnitudes or new accuracy evidence. Missing, foreign or multiply named
root symbols are errors. A rational request without roots cannot supply a germ.
Existing plain coordinate-map destinations remain unchanged for rational input.

The example transports `Y'=epsilon*Y/sqrt(s)` with `Y(1)=1` to `s=4`, then `9`,
then `9` again. It checks progressive boundary reuse, and the repeated point is
an exact cache hit. Remove `seeds` to restart from its binary bank. The first
three coefficients at `s=9` are `[1,4,8]`. To request the opposite sheet, supply
an independently justified seed on that sheet; a principal-sheet cache entry
cannot answer an opposite-sheet request, even at identical coordinates.
Algebraic result entries include `root_germ` and `starting_root_germ`; rational
output keeps its existing format.

This steers the existing `RustFlow<AlgebraicKinematicSystem>` and `RustFlowCache`.
The CLI adds no numerical solver or second bank. Its current cache scope is
exact real rational coordinates and straight paths whose radicands remain real
and nonzero. Root-sum denominators use the native registered-root normalization
and conservative formal norm domain. The global `options.prescription` does not
choose germs implicitly. `branch_domain` must still describe the intended
integral/logarithmic branch; local root signs alone do not establish it.
Threshold-crossing contours and complex physical paths require the explicit
native transport API. See [the algebraic cache documentation](algebraic-cache.md)
for uncertainty estimates, synchronous native-call limits and schema-4 binary
compatibility. The steering JSON schema remains version 1.
