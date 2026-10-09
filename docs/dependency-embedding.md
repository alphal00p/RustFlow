# Embedding the Python bindings

The optional `python` feature exports native PyO3 classes for registration by a
shared host extension. `python_stubgen` additionally enables stub metadata. This
crate does not build a second Python extension. PyO3 0.28 and the Symbolica
`python_export` API are shared with the community host.

For a fresh consumer checkout, follow the [published-input setup
recipe](clean-community-build.md). The checked-in manifests and locks fetch the
native owners directly from public Git pins. They require no local patch
application or sibling checkout. The community manifest and lock select its
RustFlow runtime revision.

The standalone root and community host each own their Cargo patch tables;
Cargo ignores a dependency's patch table. Both source HEPKit, Linnet, Spenso,
Idenso and the lean Typst components from the canonical
`https://github.com/alphal00p/gammaloop` repository, branch `feynkit`. The current
standalone lock records `635d3a1feb44067583a668c7d18f67405fe144e1`; the former
fork-based package overrides have been removed. Historical validation reports
retain their original owner revisions.

The current lock selects Symbolica, Numerica and Graphica 3.0.1 from official
community revision `ed2374f1d880d52c3a7ca48cd7c22f4baad5c020`, and RustRed revision
`78969aab524b7d6a2eec36f59b01e9af1e04cc04`. This descendant of the former
`acc92b0dad27b11fd194a4c284765fb6a93cbc94` lock adds guarded distribution-source
reduction. The update leaves the Symbolica, Numerica, Graphica and HEPKit source
identities unchanged; RustRed's standalone vendored patch table is not imported
into RustFlow. The manifest pins this revision because the guarded interface
is on RustRed's `fermi` branch; unpinned `main` is not an equivalent embedding
source. Embedding workspaces that also depend on RustRed must select this same
Git source identity. RustRed's core Cargo package is
`rustred`; its experimental `reconstruction` feature remains disabled. Native
and Pyodide consumers select one shared arithmetic backend. `python` enables
the native automatic evaluator; `python_wasm` enables its serial portable
counterpart. The core `automatic` feature is independent of `native` and can
be paired with `wasm`. See [WASM support](wasm.md).

The host additionally pins Vakint separately at
`6203c6cbba6ae5e90329ba5081fad55319e678db`; it remains native-only. Run the
root-scaling regression when changing the numerical dependency graph:

```sh
cargo test --locked --test symbolica_root_scaling
```

## Shared upstream Symbolica revision

Physical nonplanar boundary preparation exposed an exact-zero-policy bug in
Symbolica's checked division and inversion. The [minimal reproducer and
validation](symbolica-exact-division-mre.md) explain the failure. The correction
is now in official community commit `c3408e4ba1d3bdd4ea55678fad50e27009be13d4`.
The current shared revision retains that fix and the earlier polynomial
root-convergence correction, and adds the generic C++ complex-constant export
fix. The owning manifests keep all three algebra crates on the same source:

```toml
[patch.crates-io]
symbolica = { git = "https://github.com/symbolica-dev/symbolica", branch = "community" }
numerica = { git = "https://github.com/symbolica-dev/symbolica", branch = "community" }
graphica = { git = "https://github.com/symbolica-dev/symbolica", branch = "community" }
```

Each checked-in lock selects its recorded immutable commit above. Once dependencies are
fetched, `cargo metadata --locked --offline --format-version 1` checks the graph
without rebuilding. Exactly one instance of each of these three packages must
resolve from that upstream revision. The two-line local patch and isolated
`1fbdb0a` checkout are retained only as diagnostic history; neither is an active
dependency override. Normal builds use the public manifest patch tables and
checked-in locks. Remove stale development overrides from a publication build;
otherwise Cargo may select a different source than the declared public pins.

A dependency/source change creates a new numerical cache identity. Earlier
snapshots retain their original provenance and are not silently reused with
the updated graph.

Large exact coefficients exposed a separate problem in Symbolica's convenience
coefficient collector. The [standalone reproducer](../tools/mre/symbolica-coefficient-list/README.md)
shows a nonzero coefficient of order `10^309` being discarded by its statistical
zero test. RustFlow therefore selects Symbolica's exact expression field and
native polynomial grouping explicitly. This uses the same official dependency;
no local Symbolica patch is needed for the workaround.

Before this conversion, the library bounds total polynomial degree at 100,000
using Symbolica's expression traversal. This also prevents the native converter's
unsigned-to-signed exponent cast from silently wrapping very large powers.
Function arguments and opaque power coefficients retain their literal meaning.
Native conversion can expand polynomial coefficients in unrequested variables;
the degree bound is not a general bound on the number of multivariate terms.

## Public native owners

The required HEPKit extensions have been consolidated into the canonical
`feynkit` branch. It provides the graph, parameter and tensor operations
documented in [hepkit-integration.md](hepkit-integration.md). RustFlow reuses
those implementations. Historical fork pins and upstream review links in older
validation reports describe those tested artifacts, not the active source.

Native borrowing uses `PyIntegralFamily::as_family`,
`PyKinematics::as_kinematics`, `PyModel::as_model`, and the checked
`PyFeynmanDiagram::as_shared_diagram`. The latter returns a borrowed
`Arc<FeynmanDiagram>` so the evaluator can clone the shared owner before
releasing the GIL. It preserves complete-selection checks and leaves upstream
`as_diagram() -> PyResult<&FeynmanDiagram>` unchanged. No native model, diagram
or kinematics is serialized and reparsed at this boundary.

`Model::expand_parameters` performs exact transitive substitution using
Linnet's dependency ordering and rejects cyclic definitions. External and
value-only parameters remain symbolic. Both sides of each substitution are
literal Symbolica patterns, including names ending in underscores. Numeric
model defaults are not inputs to exact loop-integral evaluation.

The host's [Vakint revision `6203c6cbba6ae5e90329ba5081fad55319e678db`](https://github.com/ValentinHirschi/gammaloop/commit/6203c6cbba6ae5e90329ba5081fad55319e678db)
is based on `8d6c8f7b14f2438126819328e20a10e6caf9e0b7` and aligns its RustRed
references with the host. Its FeynKit dependencies resolve to the shared owner
through the host's patch table. Replacing Vakint with the crate from the HEPKit
owner revision would discard that distinct implementation.

For a custom embedding workspace, carry this crate's crates.io patch table
into the owning workspace and use the same canonical FeynKit branch for every
companion crate. Verify the resolved graph and lock contain one source and
revision per owner. The community host's package checks validate native and
browser dependency ownership, including disabled RustRed reconstruction.

`.cargo/config.example.toml` is optional and intended only for standalone
RustFlow development against a local native-owner checkout. It contains no
Vakint override and is not part of the public build recipe. Do not copy it over
the community host's source-fingerprint configuration. Local source overrides
require a deliberate lock update and a new cache identity; keep machine-local
paths out of published manifests and locks.

The six historical `scripts/patches/hepkit-*.patch` files retain the earlier
`9d086ce`/`fc9ee6aa5` reconstruction history. Current public builds fetch the
ported owner revision directly. Component validation and complete native
notebook acceptance remain separately recorded results.

## Source-sensitive cache identity

A Cargo build script knows its own manifest but cannot reliably discover the
workspace of the executable embedding it. Therefore the **host** must set
`RUSTFLOW_WORKSPACE_MANIFEST` to its owning manifest, for example in the host's
`.cargo/config.toml`:

```toml
[env]
RUSTFLOW_WORKSPACE_MANIFEST = { value = "Cargo.toml", relative = true }
```

Cargo resolves this path relative to the directory containing `.cargo`. For
standalone builds the current crate manifest is used automatically. Hosts must
activate the binding feature in their dependency declaration. If the outer Cargo
invocation changes host features, forward those flags explicitly; for example,
a native stub build using `--no-default-features --features python_stubgen` sets
`RUSTFLOW_WORKSPACE_FEATURES=python_stubgen` and
`RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES=1`. The build verifies the resolved
RustFlow feature set against its actual Cargo feature environment and rejects a
mismatch. This prevents a stub build from fingerprinting the default host graph. The build script
runs offline, locked Cargo metadata against this manifest and rejects a graph
that does not contain the actual current crate. Dependencies must already be
fetched, as with any offline build. Source overrides must be declared in Cargo
configuration or the manifest, not exclusively in transient `--config` command
arguments, so the nested metadata invocation sees the same graph.

The fingerprint incorporates the owning workspace manifest and lock, Cargo
configuration (including `$CARGO_HOME/config.toml`), resolved package identities and features, and actual Git/path dependency
contents, including workspace manifests and model data. Registry package
checksums come from the owning lock. The existing `RUSTRED_RUNTIME_ARITIES` build setting is also recorded without
changing RustRed's solver registry. Build directories and Git administration
are excluded. Directory change tracking catches added source files as well as
changes to existing files. Package contents are hashed with relative labels;
embedded `fixtures/gg-hg` data is included in the owned code identity, and
source edits invalidate caches even when the package version and Git HEAD are
unchanged. No fallback to guessed sibling checkouts exists.

The build also records `dependency-metadata.json` in its Cargo `OUT_DIR` for
inspection. A dependency update changes the cache identity; existing snapshots
remain subject to the normal incompatibility checks rather than being silently
reinterpreted under a different algebra kernel.
