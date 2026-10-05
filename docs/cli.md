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
must be regular in epsilon, or admit the optional exact diagonal epsilon
shearing described below. Without a root registry they have rational dependence
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
The CLI adds no numerical solver or second bank. Its regular cache scope is
exact rational-complex coordinates and affine paths with certified nonzero,
finite radicands and compatible endpoint germs. Root-sum denominators use the native registered-root normalization
and conservative formal norm domain. The global `options.prescription` does not
choose germs implicitly. `branch_domain` must still describe the intended
integral/logarithmic branch; local root signs alone do not establish it.
Threshold crossings require the explicit prescribed continuation declaration
below (or the equivalent native transport API); they are never inferred from a
germ or a branch-domain string. See [the algebraic cache documentation](algebraic-cache.md)
for uncertainty estimates, synchronous native-call limits and schema-5 binary
compatibility. The steering JSON schema remains version 1.

## Canonical logarithmic input

`canonical` is an alternative to `derivatives`. Exactly one must be present:

```json
"canonical": {
  "variables": ["s"],
  "letters": ["s+r+1"],
  "matrices": [[["1"]]]
},
"roots": {"r": "s"}
```

This represents `dY = epsilon sum_a M_a dlog(L_a) Y`. Variables and letters
are ordered; each exact constant matrix has the same square size as the basis.
The CLI constructs `CanonicalAlgebraicSystem` and passes it to the existing
native `RustFlow` cache orchestrator. Letters remain separate until the selected
physical path is substituted. Original letter/radicand holes and additional
`nonzero_conditions` still constrain the admissible source and path.

Run `rustflow transport examples/cli/canonical-transport.json` for progressive
transport of `((s+sqrt(s)+1)/3)^epsilon` from `s=1` through `4` and `9`, followed
by an exact repeat. Omitting `seeds` restarts its binary bank. Root registries
and germs follow the same rules as dense algebraic input. A canonical system
without roots uses plain coordinate-map destinations and supplies no germ.

Canonical output adds `representation: "canonical"` and `identity`, the native
content key of that exact canonical connection and its basis/domain metadata.
Canonical and dense identities remain distinct even for equivalent equations;
entries are never silently transferred between them. Existing dense rational
and algebraic input/output is unchanged. The steering schema remains 1 and the
same schema-5 binary bank can retain both representations.

## Accuracy-aware source retry

`options.max_boundary_attempts` controls how many compatible cached sources may
be tried when propagated uncertainty or numerical refinement fails the requested
accuracy. Its default is 8 and it must be positive. Sources retain the existing
cost-policy order. Only accuracy failures trigger another source; cancellation,
invalid inputs, unsupported domains and other numerical errors propagate.

A successful result that needed more than one source includes `boundary_attempts`.
Each entry records the exact restart coordinates, decimal-string MPFR cost,
source verified digits and `accepted` or `accuracy_rejected` outcome. Rejected
outcomes include their accuracy reason, and registered roots include their germ.
Ordinary single-source results, including exact hits, keep their existing JSON
fields. Failed attempts do not alter the bank; completed earlier destinations
remain saved. See [Accuracy-aware cached boundary selection](accuracy-fallback.md)
for the native diagnostics and failure semantics.

## Prescribed threshold continuation

`continuation` opts the existing transport command into the native prescribed
contour planner. It is available with rational derivatives, dense registered-root
systems and canonical logarithmic input. Without this field, the existing strict
straight-path behavior and output shape are unchanged.

```json
"continuation": {
  "kind": "prescribed_affine",
  "domain": "lower physical s half-plane; no additional winding",
  "prescriptions": [{"polynomial": "s", "side": "-i0"}],
  "unprescribed_side": "-i0",
  "homotopy_admission": "all_planner_routes_in_declared_domain"
}
```

Every field in this declaration is required. The prescription sides are exactly
`+i0` or `-i0`. Polynomials use the request namespace and physical variables only;
root symbols, epsilon, unknown variables and nonpolynomial expressions are
rejected. The native planner determines each local bypass from the restricted
polynomial and checks the exact source domains along the resulting complex
route. No second solver or cache is introduced.

`homotopy_admission` is an explicit caller assertion that **every candidate-to-
destination route generated by this policy belongs to the declared integral
homotopy domain**. It does not certify that assertion. Local root sheets do not
determine logarithmic/integral monodromy, and the CLI cannot infer global sheet
membership from endpoint coordinates. If a problem needs a more selective
route predicate, use the native `HomotopyAdmission` callback. Arbitrary loops and
additional windings are outside this affine planner interface.

Run `rustflow transport examples/cli/prescribed-transport.json`. The example
continues `Y'=epsilon*Y/sqrt(s)`, `Y(1)=1`, below `s=0` to `s=-1`, then `-4`, then
repeats `-4`. Seeds and destinations still require explicit registered-root
germs. The first three coefficients at `s=-1` are `[1, -2-2i, 4i]`; the endpoint
root is opposite to the principal square root. The next destination reuses a
verified physical cached point and the repeat is an exact hit. Removing `seeds`
restarts from the same binary bank.

Routed output additionally reports the native `identity` and canonicalized
`continuation` descriptor. The exact ordered polynomials, prescription sides,
default side and domain are part of that typed identity; changing the route
policy does not silently reuse an ordinary or differently prescribed boundary.
Canonical input also retains its existing representation field. The steering
schema remains 1 and binary banks use schema 5. Only independently checked
physical endpoints are retained; contour detours are not silently treated as
physical source points. A rejected germ or route cannot mutate that destination
in the bank, while previously completed destinations remain saved.

## Diagonal epsilon shearing

For a dense `derivatives` input, including registered square roots, `"epsilon_shearing": true` asks the existing
native adapter to find integer weights with
`I_i = epsilon^weight_i J_i` such that every physical derivative matrix for J is
regular in epsilon. Weights have maximum zero. This is one restricted gauge
class: a negative constraint cycle is a typed failure, not a claim that every
nondiagonal transformation is impossible. Explicit shearing is rejected for
canonical inputs, which are already epsilon regular. Omitting the flag keeps existing behavior.

The top-level `basis`, `leading_epsilon_power`, `last_epsilon_power`, seed values
and output coefficients continue to refer to the **original I basis**. The leading
power must be a justified global pole bound along the declared domain, not merely
the first nonzero coefficient at the source. Epsilon-singular derivatives can
generate a pole in a component that vanishes at the source.

A seed may specify `last_epsilon_power` separately; it defaults to the requested
output last order. Its arrays cover every original power from the top-level
leading bound through this explicit seed last order. To recover original outputs
through H, the transformed bank needs coefficients through
`H - min(weight)`. Supply the corresponding extra original source orders. Missing
high coefficients are rejected, including when the bank already contains a
compatible endpoint; they are never guessed to be zero. Extra supplied rectangular
source orders are retained for later compatible queries. The optional seed range
also works without shearing.

Run `rustflow transport examples/cli/epsilon-shearing-transport.json`. This example
has `I1'=I2/epsilon`, `I2'=epsilon*I1` and exact initial values
`I(0)=[0,1+epsilon]`. It declares the global bound -1 and requests original output
through zero, while supplying source coefficients through +1. The cache basis is
`J=[epsilon*I1,I2]`, with weights `[-1,0]`. The returned pole and finite coefficients
of I1 both equal sinh(s), and I2's finite coefficient equals cosh(s).

The binary bank stores the explicitly rescaled identity and ordinary native
transport results; outputs are restored to the original requested range before
a destination is saved. JSON `epsilon_shearing` metadata records both identities,
weights, original output range and requested cache range. The top-level `identity`
field names the original basis of the returned coefficients; the nested
`epsilon_shearing.cached_identity` names the binary bank basis. Original
denominator restrictions, normalization, branch domain, prescribed-contour
metadata and source uncertainty remain attached.
Original-basis entries cannot become rescaled entries by relabeling; conversion
checks identity, range and evidence. Source accuracy caps are unchanged.

A restart uses the same flag, original equations/basis/domain and requested range;
`seeds` can be omitted when the rescaled bank already covers it. The optional
`continuation` declaration applies unchanged. No new contour or numerical solver
is used by this CLI option.


Registered-square-root `derivatives` can also set `epsilon_shearing: true`.
Source and destination `root_germ` declarations remain mandatory. Source
coefficients and returned values still describe the original basis; the nested
`epsilon_shearing.cached_identity` names the transformed bank, while the top-level
`identity` names returned original values. Root sheets are not redefined by the
epsilon gauge. Raw denominator/norm restrictions and prescribed continuation
metadata remain enforced on transport and exact cache hits. `canonical` input
continues to reject explicit shearing rather than silently ignoring the option.

## Taylor step proposals

`options.step_size_strategy` accepts `"halving"` (the default) or `"bracketed"`.
Both use the same Taylor-tail, midpoint/endpoint defect, domain and branch checks.
Bracketed transport can test up to two larger steps after a first successful
halving, reusing that Taylor chart. It is opt-in while broader performance
coverage is gathered; its measured improvement on one five-point leg is not a
general speed guarantee.

`options.max_steps` limits all acceptance-predicate evaluations, including
successful candidates later superseded by a larger successful step. Transport
output includes `predicate_evaluations`, `steps` (committed segments),
`rejected_steps`, and `superseded_successes`. An exact cache hit has zero transport
counts. Cache identity and accuracy requirements are unchanged; either strategy
can use a mathematically compatible, sufficiently accurate cached boundary.
