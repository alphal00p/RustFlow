# Native HEPKit integration

RustFlow accepts HEPKit's model JSON and native compact or serialized DOT directly through `FeynmanDiagram::from_dot`. HEPKit/Linnet own parsing, edge identities, graph validation, routing, denominator construction and partial fractions. Spenso/Idenso own scalar numerator contraction. This crate translates the resulting complete denominator basis into `IntegralFamily` and exact `LinearCombination` values. The application owns one loaded UFO model at a time and shares its `Arc<Model>` with all graph inputs; RustFlow does not maintain a competing model registry. Python bindings are optional modules in this crate, registered by the shared community extension under `symbolica.community.hep.integration`; the library and CLI remain usable without Python. See [dependency-embedding.md](dependency-embedding.md) for the current shared dependency graph and source fingerprint configuration.

Current builds select the public [native owner revision
`a3d1c8a867ac0e89d8f58f9722cf9a7e83338901`](https://github.com/ValentinHirschi/gammaloop/commit/a3d1c8a867ac0e89d8f58f9722cf9a7e83338901),
based on upstream HEPKit `6c707c6b77a437256eb1180da13d4d327b371d13`. The standalone
and community manifests contain the full shared-owner patch tables. The
[published-input build recipe](clean-community-build.md) uses those pins
directly, with no private checkout or patch application. Symbolica, Numerica
and Graphica resolve from official community `c3408e4ba1d3bdd4ea55678fad50e27009be13d4`;
RustRed resolves from official main `7c1ed03722b8c05daf60c89ba4ecc79457ed2ada`.
The community host preserves its separate Vakint implementation at
`6203c6cbba6ae5e90329ba5081fad55319e678db`.
The native owner changes are proposed upstream in
[GammaLoop PR #125](https://github.com/alphal00p/gammaloop/pull/125).

The native owner supplies:

- Literal momentum-square substitutions, so scalar and dimension names ending
  in `_` remain data during graph-family construction.
- Closed single-vector angular moments through even rank 32. General,
  external-basis and block projectors retain rank 20; mixed rank 22 and
  single-vector rank 34 return typed errors.
- Concurrent model-parameter registration that preserves existing printers,
  invalid-name errors, and independent LaTeX/Typst labels.
- Native directed strongly connected components for differential-system blocks.
- Shared model and diagram borrows for the Python adapter, retaining upstream
  APIs and complete-diagram selection checks.
- Exact transitive analytic parameter expansion, with native dependency ordering,
  literal parameter names and typed cycle errors. Value-only model defaults
  remain symbolic.

FeynKit graph/model/kinematics/tensor are version 0.1.0, Linnet 0.17.0,
Spenso 0.6.0 and Idenso 0.3.0. Ordinary Rust builds keep Python features disabled.
The owner revision preserves the newer native external-wavefunction and
rendering APIs. RustFlow's build fingerprint includes the resolved source
identities, consumed source and manifest contents, and embedded model data.
Cache identity follows the implementation that was built, rather than only
package version numbers.

The six `scripts/patches/hepkit-*.patch` files preserve the earlier owner setup
and its isolated regression history. Early runs used HEPKit `9d086ce` with local
fixes; still earlier native validation used `8f834d9c62ae06fb327e4ef0b14abffda755b610`
and Symbolica `75f8350094b90254ee71dc2a391fde0d14b0204a`. Those recorded checks
retain their original source provenance. Current installation uses the public
owner pin above; historical numerical evidence does not by itself certify a
new dependency graph's full notebook acceptance.

## Rust interface and CLI conventions

```rust,ignore
let model = std::sync::Arc::new(feynkit_model::Model::from_json(&model_json)?);
let dimension = symbolica::parse!("feynkit_graph::D");
let kinematics = feynkit_kinematics::Kinematics::in_dimension(&dimension)?
    .with_mass_squared(&feynkit_graph::symbols::external_momentum().call(1),
                       symbolica::Atom::num(-1))?;
let graph = symbolica_amflow::hepkit::GraphIntegral::from_dot(
    model, &dot, &kinematics,
)?.with_powers(&edge_powers)?;
let groups = graph.integral_groups(&exact_point, epsilon, 4, 10_000, &context)?;
let result = symbolica_amflow::solve_integral_combinations(
    &groups, &Default::default(), 0, &options, &backend, &context,
)?;
```

`with_powers` accepts a map from native `EdgeId` to signed `i16` powers. Powers are integral metadata because the native DOT schema does not define propagator powers. Graph numerator, vertex and edge numerators, projector and overall factor are included exactly once; graph automorphism factors are not divided out implicitly.

Native tensor dimensions must be symbolic, for example `feynkit_graph::D`. Use that same symbol in the original DOT tensor slots (`spenso::mink(D,mu)`). After contraction the adapter replaces it by `D0 - 2*epsilon`. Surviving Lorentz slots with another dimension are rejected rather than implicitly changed. Native parsing can contract a metric trace immediately: an explicitly four-dimensional trace becomes the scalar `4` and remains that scalar. Its original tensor dimension cannot be recovered afterwards.

Native `Kinematics` provides `with_mass_squared` and `with_scalar_product` setters for unindexed momentum names or labeled calls. Canonical momentum heads are `gammalooprs::P`, `gammalooprs::K` and `gammalooprs::Q`; the massless bubble fixture uses `gammalooprs::P(1)`. These setters do not solve constraints on sums of momentum labels. Loop scalar products must remain unconstrained integration variables.

Unqualified DOT expression symbols belong to `feynkit_graph`; model formulas belong to `UFO`. Fully qualified names remove ambiguity. Exact scalar substitutions use symbol keys. Model `f64` default values are never used as exact inputs; required parameter values must be provided explicitly. The adapter specializes the supplied point before partial fractions, so pass an empty `KinematicPoint` when solving its returned groups. This preserves simultaneous substitution semantics even for maps such as `a -> b, b -> 3`.

HEPKit and RustRed order scalar products differently when there is more than one loop. The adapter matches exact Atom identities: loop-loop products first, followed by loop-external products in RustRed order. Native complex affine algebra uses a rational coefficient field; RustFlow reuses its formal imaginary-unit bridge for exact complex kinematics, then restores Gaussian-rational coefficients before evaluation. No arbitrary-precision value is routed through `f64`.

## Tensor projection and supported scope

`TensorProjector` retains its Gram-matrix interface for boundary recursion, but now labels those vectors in Spenso notation and calls native `TensorReducer`. Only the translation and a cache of native projection results remain here. HEPKit owns contraction orbits, pairing counts and coefficient inversion. A mixed rank-20 case with multiplicities `[2,18]` reduces to two invariants; the former RustFlow implementation stopped at mixed rank 14.

The adapter supports scalar numerators and quadratic model denominators with exactly zero widths. It rejects free tensor indices without a scalar projector, nonstandard model denominator formulas, inexact substitutions and loop kinematic constraints. Cut metadata is preserved by the native graph. Ordinary `integral_groups` rejects it; `cut_integral_group` converts one explicitly selected native `DiagramCut` with oriented momenta and unchanged denominator slots. Dependent cut denominators are rejected until their decomposition can preserve measure metadata. The supported two-body phase-space terminal and native cut IBPs are described in [cuts.md](cuts.md); a cut is never silently treated as an uncut propagator. Native compact `is_cut` attributes identify matched dangling legs, rather than Boolean labels on arbitrary internal edges.

The fixtures in `fixtures/hepkit` are small model and graph inputs, not copied computational implementations. FeynKit and Spenso declare `MIT OR Apache-2.0`; Linnet declares MIT. Some companion manifests, including Idenso, do not declare a separate SPDX license; the repository's top-level notice states that GammaLoop has no usage restrictions. Symbolica retains its independent license. RustFlow's MIT license does not relicense these dependencies.

## Native graph algorithms

Gaussian Wick connectivity now uses Linnet's existing `UnionFind`; only the physical endpoint/Gram interpretation remains in RustFlow. Differential-system blocks are a native `HedgeGraph` whose directed edges are nonzero matrix dependencies. The pinned public Linnet owner supplies `strongly_connected_components` and its subgraph variant; `scripts/patches/hepkit-linnet-strongly-connected.patch` records its earlier implementation history. It traverses native half-edge incidence directly, includes isolated nodes, respects underlying/superficial direction and returns sink components first with stable node-order ties. RustFlow no longer carries its dense transitive-closure graph algorithm. The standalone native-API probe covers all 512 directed three-node graphs, reversed/undirected/split/dangling edges, a 20,000-node cycle and both supported node-storage implementations; integration tests retain the solver's block order.

The adapter assumes one active UFO model, as HEPKit does. Native symbol registration is process-wide; the concurrency fix permits simultaneous imports of that model and preserves previously registered printers. It does not claim isolation between incompatible UFO models sharing parameter names. Momentum bases continue to come from native `LoopMomentumBasis`; no second graph or routing representation is maintained here.

## Automatic physical boundary seeds

`PreparedPhysicalFamily` closes all declared physical partial derivatives on one ordered master basis and can fill `RustFlowCache` with automatic AMF or FT master values. Its identity retains canonical original family coefficients and symbol names, so anonymous reducer parameter names cannot merge different normalizations. Every remaining physical parameter must be declared or specialized first. Seed accuracy comes from independent epsilon grids and refinement, and subsequent transport carries that uncertainty forward.

Reduction nonzero conditions and original matrix denominator/branch restrictions travel with cache identity and persistence. A pure epsilon guard is valid generically; `epsilon*s` excludes `s=0`, while the reduction guard `s-epsilon` does not by itself exclude that point. A matrix denominator `1/(s-epsilon)` separately requires a regular Laurent chart and therefore excludes `s=0`. Exact polynomial gcd and native root analysis filter unsafe straight paths before comparing cached starting points, including poles hidden by stationary path coordinates. Rational-complex guarded paths are supported; conditions requiring unsupported algebraic coefficient fields return a typed error. Physical boundary cache codec version 2 rejects snapshots lacking these restrictions.

## Python interface and current coverage

The optional `python` feature supplies bindings for the community host's existing native Symbolica extension. `python_stubgen` adds the inventory used by that host to generate type stubs. This crate does not build a second Python extension, and these bindings are excluded from browser builds. All new Python classes and exceptions are registered under `symbolica.community.hep.integration`; existing HEPKit and Hyperbolica APIs retain their namespaces. Build/source-identity requirements and the coordinated dependency graph are documented in [dependency-embedding.md](dependency-embedding.md).

```python
from symbolica.community.hep.integration import (
    IntegralEvaluator, EvaluationOptions, BoundaryCache, KinematicTransport,
)

evaluator = IntegralEvaluator(
    options=EvaluationOptions(digits=20), reduction_batch_size=32,
)
```

The evaluator accepts native HEPKit `IntegralFamily`, `FeynmanDiagram` and `Kinematics` objects; `HiggsJetAmplitude` reuses the supplied native `Model` for generation, diagram display and scalar assembly. Inputs and outputs use the shared Symbolica `Expression`, `Float` and `ComplexFloat` classes. Exact kinematic substitutions remain simultaneous, and arbitrary-precision values are never converted through binary64.

`IntegralEvaluator` exposes automatic evaluation, shared linear projections, prescribed nonzero rational epsilon samples and physical-family preparation. `PreparedIntegralFamily` exposes automatic boundary generation, physical transport, basis definitions, reduction conditions and target projection. `ReductionTables` are scoped to native families. `LaurentExpansion.fit` performs a supplied-sample fit, whose evidence is explicitly unverified until independent checks are provided by a higher-level evaluation. Exact and completed-sample caches are configured separately through `EvaluationOptions`; verified physical boundary values belong to `BoundaryCache`.

`BoundaryCache.load`, `save`, `entries` and `extend` preserve native binary entries, identity and provenance. Passing the same cache through successive transports retains intermediate points for subsequent source selection. `extend` merges seed data into an existing growing cache without discarding those points; self-merges are idempotent and source/destination locks are never held together. Transport results expose the selected source, attempted sources, timing, propagated uncertainty and exact cache-hit diagnostics.

The standalone Python `DifferentialSystem` currently exposes rational regular-point continuation with supplied `BoundaryData`, block decomposition and rational basis changes. Its fixed-precision transport reports diagnostics and retains input accuracy/provenance, but does not certify a global output error bound. Native Frobenius initialization, infinity expansions and dimensional endpoint selection are not yet separate Python methods. Generic `KinematicTransport` and prepared-family transport currently expose admitted straight paths with native singularity/branch checks; custom physical route objects are not wrapped.

`evaluate_diagram` covers the ordinary scalar graph workflow. A graph carrying cut metadata raises `UnsupportedInputError`; the Rust `cut_integral_group` workflow has not yet been exposed to Python. This preserves the distinction between a cut measure and an ordinary propagator.

Computation releases the GIL. `ComputationControl.poll()` drains native progress events and `cancel()` requests cooperative cancellation. The Python reducer defaults to batches of 32 targets, configurable with `reduction_batch_size`; cancellation is checked between reducer batches and native algorithm stages. A running HEPKit/Idenso tensor operation finishes before acknowledging cancellation. No Python callback runs while a cache lock is held.

The Higgs-jet scalar kernel delegates exact internal-parameter expansion to native `Model::expand_parameters`, including the owner's dependency ordering and typed cycle detection. This API is included in the public owner pin; [hepkit-exact-parameters.patch](../scripts/patches/hepkit-exact-parameters.patch) preserves the earlier isolated-owner change. Model numeric defaults and reference integral boundaries are not inputs to this kernel. Form-factor kinematics and W/Z masses must match the explicit model parameters; the notebook constructs both from one set of exact physical inputs.

## Remaining numerical specialization

Exact matrix determinant, multiplication and inversion already use Symbolica's matrix facilities. At official Symbolica/Numerica revision `c3408e4ba1d3bdd4ea55678fad50e27009be13d4`, native `Float::{add,sub,mul,div}_round` accepts an explicit precision and rounding direction. `numeric.rs` uses those primitives with `RoundingDirection::Nearest` for the real components of complex arithmetic, retaining the order and working precision of each intermediate operation. Output accuracy continues to come from independent checks, not from the working precision.

Floating matching and reusable Frobenius block solves still require numerical policies absent from the native matrix APIs. Numerica's `lib/numerica/src/tensors/matrix.rs::partial_row_reduce`/`solve` selects the first nonzero pivot, whereas these solves use magnitude pivoting. `FloatField` delegates to ordinary `Float` arithmetic, whose `fixed_precision()` remains false; it does not select the explicit rounding primitives above. The sparse native row reducer exposes LU/RREF factors, but no magnitude-pivot policy or reusable numerical right-hand-side solver. Replacing these routines therefore needs an owner-level field/pivot policy and reusable-factor solve API, followed by ill-conditioned matching and resonance validation. A generic field adapter alone would not preserve their behavior.
