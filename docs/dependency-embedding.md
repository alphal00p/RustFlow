# Embedding the Python bindings

The optional `python` feature exports native PyO3 classes for registration by a
shared host extension. `python_stubgen` additionally enables stub metadata. This
crate does not build a second Python extension. PyO3 0.28 and the Symbolica
`python_export` API are shared with the community host.

The public dependency declarations use the host's source identities: released
Symbolica 3.0.1, GammaLoop's `feynkit` branch, and RustRed's `main` branch with its
experimental reconstruction feature disabled. The workspace root owns Cargo's
patch table; dependency patch tables are ignored by Cargo. The standalone lock
pins Symbolica's community revision `942bd2c0cd2ef16414d69c176fc9eff7b21c0ab2`
and RustRed `b3cecd6ae9b7683ae639204e43fa05685c6a5802`. Run the root-scaling
regression when changing the numerical dependency graph:

```sh
cargo test --locked --test symbolica_root_scaling
```

## Native owner patches

The native owner fixes documented in [hepkit-integration.md](hepkit-integration.md)
are still required, including Linnet's strongly connected components. The Python
integration uses the HEPKit source at `9d086cec7971005ec7244b43e8fcea2a403c18ef`
with those fixes and native borrowing accessors. These changes belong to HEPKit;
this crate does not duplicate graph or tensor implementations. Until they are
upstream, configure Cargo to use the patched owner checkout. Copy
`.cargo/config.example.toml` to `.cargo/config.toml` and replace the checkout
paths. The template lists all shared HEPKit packages, avoiding duplicate native
Rust types when a host also uses Spynso3, Vakint, or the RustRed bridge.

The borrowing accessors are supplied by
`scripts/patches/hepkit-python-borrow-access.patch`. They expose references to the
existing model, diagram and kinematics; no Python objects are serialized or
reparsed when passed to an evaluator. The isolated validated owner checkout is
commit `0bf1cd991` on `codex/python-borrow-access`, based on `9d086ce`. The original
HEPKit checkout remains unchanged.

Local source overrides require updating the lock once with `cargo metadata` or
`cargo check`; subsequent checks can use `--locked`. Do not commit machine-local
configuration or claim an unpatched upstream checkout supplies these owner
extensions. No special sibling directory layout is required.

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
checksums come from the owning lock. Build directories and Git administration
are excluded. Directory change tracking catches added source files as well as
changes to existing files. Package contents are hashed with relative labels;
embedded `fixtures/gg-hg` data is included in the owned code identity, and
source edits invalidate caches even when the package version and Git HEAD are
unchanged. No fallback to guessed sibling checkouts exists.

The build also records `dependency-metadata.json` in its Cargo `OUT_DIR` for
inspection. A dependency update changes the cache identity; existing snapshots
remain subject to the normal incompatibility checks rather than being silently
reinterpreted under a different algebra kernel.
