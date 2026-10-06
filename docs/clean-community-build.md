# Building the community extension from published inputs

Use a community checkout whose `Cargo.toml` and `Cargo.lock` contain the public
integration dependencies below. Those files select the RustFlow revision;
keep them together when rebuilding or sharing an environment. Activate the
community repository's native Rust/Python build environment and run these
commands from its root. Python 3.11 or newer is suitable for the setup and
notebook tools.

The API and notebook are on community `main` after
[squash-merging PR #19](https://github.com/symbolica-dev/symbolica-community/pull/19).
To obtain the validated source snapshot in a new directory:

```sh
git clone --branch main \
  https://github.com/symbolica-dev/symbolica-community.git symbolica-community-integration
cd symbolica-community-integration
git checkout --detach 66a67bb0bc7b01ae53757c274aba92fc2ef234f8
```

This snapshot aligns the host with RustFlow `92cfc9d`, official Symbolica
community `6defcca9` and native HEPKit owner `b96600b`, including regenerated
integration stubs. The earlier wheel and controller completed full cold,
independent-refinement, restart and nearby-point acceptance on `8758b93d`;
that evidence retains its original dependencies and numerical cache identity.
See [current acceptance status](python-notebook-status.md) for exact provenance
and remaining upstream/CI checks. Merging the source does not publish a package.

Normal builds fetch the required native owners directly from public Git pins.
They require no sibling checkout, private owner commit, patch application, or
machine-local Cargo source configuration.

| Component | Published source selected by the manifest and lock |
| --- | --- |
| RustFlow | The community manifest's `symbolica-amflow` Git revision and matching lock entry |
| HEPKit, Linnet, Spenso, Idenso and native rendering | [ValentinHirschi/gammaloop `b96600b0085d9ddfa9e6acbc11fa72ec6163253c`](https://github.com/ValentinHirschi/gammaloop/commit/b96600b0085d9ddfa9e6acbc11fa72ec6163253c) |
| Vakint | [Separate revision `6203c6cbba6ae5e90329ba5081fad55319e678db`](https://github.com/ValentinHirschi/gammaloop/commit/6203c6cbba6ae5e90329ba5081fad55319e678db), based on upstream `8d6c8f7b14f2438126819328e20a10e6caf9e0b7` |
| Symbolica, Numerica and Graphica | [Official community revision `6defcca968ca8411977fb1f641a9dee49ee7b7a7`](https://github.com/symbolica-dev/symbolica/commit/6defcca968ca8411977fb1f641a9dee49ee7b7a7) |
| RustRed | [Official main revision `7c1ed03722b8c05daf60c89ba4ecc79457ed2ada`](https://github.com/alphal00p/rustred/commit/7c1ed03722b8c05daf60c89ba4ecc79457ed2ada) |

The Cargo package for RustRed's core is `rustred`. Its experimental
`reconstruction` feature stays disabled. The community bridge retains
`campaign-api` and uses the same RustRed source as Vakint and RustFlow.

The native owner revision is based on upstream HEPKit `6c707c6b77a437256eb1180da13d4d327b371d13`.
The manifests contain the complete 22-package patch table from the original
GammaLoop source to that public fork revision. The community host also declares
FeynKit's Python bindings and Spynso3 directly at that revision. Its lean Typst
library and SVG crates use the same pin. Vakint retains its separate revision;
its dependency alignment does not replace its implementation with the HEPKit
checkout's Vakint crate.

## Native installation

Provide `SYMBOLICA_LICENSE` through your environment when required. Keep license
material outside the checkout. With `maturin` and `pytest` available in the
activated Python environment:

```sh
cargo fetch --locked
python .github/scripts/check_integration_package.py
RUSTFLOW_WORKSPACE_FEATURES=pyo3/extension-module RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES=0 \
  maturin develop --release --locked --extras notebook-display
python -m pytest -q tests/test_loop_integration.py tests/test_loop_integration_reductions.py
```

The package check inspects the active native, stub-generation and browser
dependency graphs. It checks shared native owners and browser exclusions; it
does not execute a numerical acceptance calculation. Preserve the host's tracked
`.cargo/config.toml`, which declares `RUSTFLOW_WORKSPACE_MANIFEST` for the build
fingerprint. A published build should resolve RustFlow from its declared Git
revision, without a development path override in the checkout or `CARGO_HOME`.
The feature forwarded above matches the extension-module feature enabled by
Maturin in `pyproject.toml`, so the cache fingerprint includes the build's
actual Python linkage features.

The tracked Python stubs can be regenerated from the same dependency graph:

```sh
RUSTFLOW_WORKSPACE_FEATURES=python_stubgen RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES=1 \
  cargo run --release --locked --no-default-features --features python_stubgen \
  --bin stub_gen -- --hepkit-only
```

These host features matter to the source fingerprint. See
[dependency embedding](dependency-embedding.md) for custom hosts and optional
standalone development overrides. Ordinary library builds keep Python disabled.

## Notebook and acceptance

The notebook is `examples/hep/gg_hg.py` in the community repository. Install
Marimo in the same Python environment and launch it from that checkout:

```sh
python -m pip install 'marimo>=0.24.2,<0.25'
marimo edit examples/hep/gg_hg.py
```

The `notebook-display` extra above provides the native display widget. Routine
controller and notebook tests are separate from the long calculation:

```sh
python -m pytest -q tests/test_hep_gg_hg.py tests/test_hep_gg_hg_anchors.py tests/test_hep_notebooks.py -k gg_hg
python examples/hep/gg_hg_anchor_acceptance.py --directory /path/to/new/anchor-bank \
  --workers 4
python examples/hep/gg_hg_acceptance.py --directory /path/to/new/acceptance-bank \
  --anchor-report /path/to/new/anchor-bank/anchor-acceptance.json \
  --workers 4 --boundary-workers 1 --seed-digits 30 --interrupt-after-samples 1
```

Use a new output directory for each cold run. The anchor runner independently
generates the planar and nonplanar Euclidean boundaries before opening either
comparison fixture. Its report checks all 545 coefficients against their
recorded reference uncertainties. The full runner verifies that the supplied
anchor report matches the installed extension and steering sources before
starting the physical calculation. Keep that environment frozen between runs.

The long acceptance checks native
boundary generation, all 4,360 reference transport coefficients, form factors,
coherent observables, independent refinement, and cache restart. The optional
interruption occurs during the later forced refinement stage, preserving a
complete cold-stage timing. Worker settings are resource choices; they do not
change the requested accuracy.

Freeze the installed extension and steering sources during a run. The acceptance
report records their hashes and the native source-sensitive identity. A new
dependency graph produces a new cache identity, so earlier caches and numerical
reports keep their original provenance. The commands above describe the required
validation; availability of public dependencies and successful component tests
alone do not establish a completed clean-machine build or full native notebook
acceptance for this graph.

## Historical patch reconstruction

The six `scripts/patches/hepkit-*.patch` files record the earlier owner setup.
On 2026-10-05 they were applied from RustFlow `6e86afd` to a disposable archive
of HEPKit `9d086ce`, reproducing twelve files from the then-validated local owner
`fc9ee6aa5` byte for byte. That historical source audit is recorded in
`target/published-hepkit-patch-audit.json`. Current builds fetch the public owner
revision above, which carries the ported fixes and preserves the newer HEPKit
APIs. The historical patch sequence is not an installation step.
