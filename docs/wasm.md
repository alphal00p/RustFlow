# Automatic evaluation in WebAssembly

The automatic integration code can use RustRed's portable backend. Reduction,
differential-equation closure, recursive boundaries, finite-epsilon sampling and
Laurent reconstruction use the same Rust implementation as native builds.
Symbolica supplies Malachite integers, Astro arbitrary-precision floats and Gamma
evaluation. No new numerical backend or Gamma approximation is implemented here.
HEPKit, Linnet, Spenso and Idenso continue to own graph and tensor operations.

| Cargo selection | Capability |
| --- | --- |
| Default features | Native GMP/MPFR arithmetic, RustRed and automatic evaluation |
| `--no-default-features --features wasm` | Portable supplied-boundary transport |
| `--no-default-features --features wasm,automatic` | Portable transport and automatic evaluation with RustRed |
| `--no-default-features --features python_wasm` | The portable automatic API registered in the community host's shared Symbolica extension |
| `--no-default-features --features python` | The native automatic Python API |

`automatic` selects reduction and boundary algorithms independently of the numeric
backend. `native` and `wasm` remain mutually exclusive. The CLI requires
`automatic`; enabling that feature does not turn a CLI executable into a browser
application. The host's Pyodide/Emscripten pipeline builds and packages the shared
extension. `python_wasm` does not build a second standalone Python module.

The intended browser runtime is **Pyodide, targeting
`wasm32-unknown-emscripten`**. Selecting `wasm` on Linux is useful for testing the
portable arithmetic implementation but is not evidence of browser execution.
This support does not claim a working `wasm32-unknown-unknown` or WASI host.
Source fingerprints must describe the embedding host's actual feature selection,
using `RUSTFLOW_WORKSPACE_MANIFEST`, `RUSTFLOW_WORKSPACE_FEATURES` and
`RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES=1` when applicable.

## Scheduling and cancellation

Browser evaluation requires `FlowOptions.workers = 1` and
`RustRedBackend.native_workers = 1`. Requesting additional workers returns a typed
unsupported-input error. Work executes inline; there is no automatic browser
thread pool. Releasing Python's GIL does not yield to JavaScript's event loop.

Rust progress callbacks run synchronously through `RunContext`. A callback can
cancel its context's token, which is checked at cooperative cancellation points.
Python's `ComputationControl` can cancel before execution and retains queued
progress for `poll()`. A UI cannot deliver a new cancellation event while its
single event-loop thread is blocked inside a synchronous computation. Hosts that
need responsive buttons must schedule work in their own supported worker or
divide it into calls and process events between them; do not introduce a Python
`ThreadPoolExecutor` in the Pyodide notebook.

## Caches and practical scope

Symbolic caches, reduction checkpoints, completed epsilon samples and
`BoundaryCache.save()` use Pyodide's virtual filesystem. Atomic cache writes and
compatibility checks still apply there. Files survive subsequent calls in the
same runtime, but do not automatically survive a page reload. Export the cache
files or arrange host-backed persistent storage, then restore those files before
`BoundaryCache.load()` in a new runtime.

Cache identities remain sensitive to their algorithm, source dependencies and
numeric backend. Do not assume a native binary cache can be opened in the browser.
For cross-backend boundaries, use the checked exact-ratio scientific bundles
described in [portable-higgs-boundaries.md](portable-higgs-boundaries.md).

Automatic boundary generation can take substantially more time and memory than
transport from supplied boundaries. Availability of the automatic API is not a
claim that a multiloop boundary calculation is suitable for a live, single-core
notebook. The gg→hg demonstration should continue to load its documented supplied
boundaries and compute physical transport interactively. Search, recursion and
precision limits remain explicit, with typed failures for incomplete work.

The regression suites `gamma`, `portable_backend`, `portable_automatic`,
`end_to_end` and `native_physical` cover owner-supplied Gamma evaluation, worker
validation, automatic tadpoles/bubbles, physical prescription and the HEPKit graph
to AMF seed to growing transport-cache path. Record actual target and runtime
results separately; these test names alone do not establish that a particular
browser wheel passed them.
