# Native HEPKit integration

RustFlow accepts HEPKit's model JSON and native compact or serialized DOT directly through `FeynmanDiagram::from_dot`. HEPKit/Linnet own parsing, edge identities, graph validation, routing, denominator construction and partial fractions. Spenso/Idenso own scalar numerator contraction. This crate translates the resulting complete denominator basis into `IntegralFamily` and exact `LinearCombination` values. The application owns one loaded UFO model at a time and shares its `Arc<Model>` with all graph inputs; RustFlow does not maintain a competing model registry. Python bindings remain in HEPKit; the RustFlow library and CLI are native Rust.

The sibling dependency checkout is `/common/dev/hepkit`, based on upstream [commit 8f834d9c62ae06fb327e4ef0b14abffda755b610](https://github.com/alphal00p/gammaloop/commit/8f834d9c62ae06fb327e4ef0b14abffda755b610). It uses the same Symbolica/Numerica/Graphica revision as RustFlow, `75f8350094b90254ee71dc2a391fde0d14b0204a`. FeynKit graph/model/kinematics/tensor are version 0.1.0, Linnet 0.17.0, Spenso 0.6.0 and Idenso 0.3.0. Native dependencies have Python features disabled.

Three local dependency fixes are maintained as reproducible patches under `scripts/patches`:

- `hepkit-literal-substitution.patch`, local commit `3f392e22bfd13e4647cf106440dfd22df5607a20`, makes the momentum-square substitution literal. Without this change, a valid scalar or dimension name ending in `_` is interpreted as a wildcard and graph family construction can panic. The isolated same-Symbolica probe reproduced the panic and then passed after the fix.
- `hepkit-isotropic-rank32.patch`, local commit `3311f81d36ff39d35de22bc8c409bdf462bff8c4`, permits the existing closed native single-vector angular moment through even rank 32. General, external-basis and block projectors retain rank 20. This preserves RustFlow's previously supported single-vector boundary numerators while moving projection to HEPKit. The independent checks include the exact rank-32 angular moment and typed rejection of mixed rank 22 and single-vector rank 34.
- `hepkit-parameter-registration.patch`, local commit `cb1098330687c3def94c8394c5380e1044ec83c5`, reuses a parameter symbol when a concurrent model import creates it between lookup and registration. This preserves existing user printers and invalid-name errors. An independent 16-thread stress probe reproduced six failures in 1,024 registrations before the fix and zero afterwards; the parallel native-DOT integration tests also pass.

For a fresh checkout, clone the upstream repository into the sibling `hepkit` directory, check out the base commit, and apply these three fixes and the Linnet extension below with `git apply` from that directory. RustFlow's build fingerprint includes each consumed path crate's source and manifest, the HEPKit workspace manifest, and the embedded model data. Cache identity therefore follows the actual implementation, including local fixes, rather than only path-package version numbers. No RustRed modification is required.

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

Gaussian Wick connectivity now uses Linnet's existing `UnionFind`; only the physical endpoint/Gram interpretation remains in RustFlow. Differential-system blocks are a native `HedgeGraph` whose directed edges are nonzero matrix dependencies. A small owner-level Linnet extension supplies `strongly_connected_components` and its subgraph variant; the reproducible patch is `scripts/patches/hepkit-linnet-strongly-connected.patch` (local HEPKit commit `a33254faf0ac960c0589891ffdd2d8f38d334b05`). It traverses native half-edge incidence directly, includes isolated nodes, respects underlying/superficial direction and returns sink components first with stable node-order ties. RustFlow no longer carries its dense transitive-closure graph algorithm. The standalone native-API probe covers all 512 directed three-node graphs, reversed/undirected/split/dangling edges, a 20,000-node cycle and both supported node-storage implementations; integration tests retain the solver's block order.

The adapter assumes one active UFO model, as HEPKit does. Native symbol registration is process-wide; the concurrency fix permits simultaneous imports of that model and preserves previously registered printers. It does not claim isolation between incompatible UFO models sharing parameter names. Momentum bases continue to come from native `LoopMomentumBasis`; no second graph or routing representation is maintained here.

## Automatic physical boundary seeds

`PreparedPhysicalFamily` closes all declared physical partial derivatives on one ordered master basis and can fill `RustFlowCache` with automatic AMF or FT master values. Its identity retains canonical original family coefficients and symbol names, so anonymous reducer parameter names cannot merge different normalizations. Every remaining physical parameter must be declared or specialized first. Seed accuracy comes from independent epsilon grids and refinement, and subsequent transport carries that uncertainty forward.

Reduction nonzero conditions and original matrix denominator/branch restrictions travel with cache identity and persistence. A pure epsilon guard is valid generically; `epsilon*s` excludes `s=0`, while the reduction guard `s-epsilon` does not by itself exclude that point. A matrix denominator `1/(s-epsilon)` separately requires a regular Laurent chart and therefore excludes `s=0`. Exact polynomial gcd and native root analysis filter unsafe straight paths before comparing cached starting points, including poles hidden by stationary path coordinates. Rational-complex guarded paths are supported; conditions requiring unsupported algebraic coefficient fields return a typed error. Physical boundary cache codec version 2 rejects snapshots lacking these restrictions.

## Remaining numerical specialization

Exact matrix determinant, multiplication and inversion already use Symbolica's matrix facilities. Floating matching and reusable Frobenius block solves still use `numeric.rs` because the current native APIs have different numerical contracts. At Symbolica revision `75f8350`, `numerica/src/tensors/matrix.rs::partial_row_reduce`/`solve` select the first nonzero pivot, whereas these numerical solves require magnitude pivoting. `FloatField` delegates to ordinary `Float` arithmetic, whose `fixed_precision()` is false; ODE recurrences instead round each operation to a chosen working precision. The sparse native row reducer owns exact LU/RREF infrastructure, but does not expose the needed magnitude-pivoted reusable numerical solve. Replacing these routines therefore needs an owner-level field/pivot policy and reusable-factor API, followed by ill-conditioned matching and resonance validation. A generic field adapter alone would not preserve their behavior.
