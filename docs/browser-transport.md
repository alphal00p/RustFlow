# Supplied-boundary transport in browser builds

The optional browser build uses Symbolica's official Malachite integer and Astro floating-point backends and the existing HEPKit, Linnet, Spenso and Idenso browser features. It retains the same canonical transport, singularity and uncertainty checks, supplied-boundary validation, intermediate-point cache, Higgs-plus-jet form-factor projection and native-object amplitude assembly algorithms.

| Features | Capabilities |
| --- | --- |
| Default | Full native library, GMP/MPFR and automatic RustRed reduction and boundary generation |
| `python` | Existing complete native Python API, including when default features are disabled |
| `python_stubgen` | Existing native API and generated stubs |
| `native`, without defaults | GMP/MPFR transport and supplied boundaries; no automatic reduction |
| `wasm`, without defaults | Portable transport and supplied boundaries; no automatic reduction |
| `wasm,automatic`, without defaults | Portable transport, RustRed reduction and automatic boundaries |
| `python_wasm`, without defaults | Portable automatic Python API and Higgs-plus-jet classes for the shared community extension |

`python_api` is the backend-neutral binding implementation used by the two Python entry features. Native and browser backends cannot be selected together. `automatic` selects RustRed independently of the numeric backend. The CLI requires `automatic`. See [wasm.md](wasm.md) for automatic browser evaluation, serial execution and persistence constraints.

`symbolica.community.hep.integration.automatic_boundary_generation_available` and the property of the same name on `HiggsJetIntegralSystem` report the selected capability. Builds without `automatic` do not export `IntegralEvaluator`, `PreparedIntegralFamily`, `ReductionTables`, or `HiggsJetIntegralSystem.generate_boundary`. Integral-family validation still belongs to RustRed: asking for it without that feature returns a typed unsupported-input error. Gamma evaluation uses Symbolica's selected numeric backend; it is not needed by the supplied canonical boundary, transport, form-factor or amplitude path.

Import precomputed values through `HiggsJetIntegralSystem.kinematic_transport().add_boundary()` after checking `mathematical_fingerprint`, as described in [portable-higgs-boundaries.md](portable-higgs-boundaries.md). Exact ratios and their recorded precision preserve supplied values across numerical backends. Runtime-specific binary caches retain their original compatibility checks; a native binary cache is not a portable input file.

The browser feature enforces one numerical worker. HEPKit's existing browser graph generation is serial. Python hosts must schedule the stages without creating a `ThreadPoolExecutor`; releasing the GIL does not create browser threads or yield control to its event loop. Progress and cancellation remain cooperative, and a host must arrange event-loop access between synchronous calls. Binary files in Pyodide's virtual filesystem support restart within that filesystem; durable browser storage is a separate host concern.

For an embedding community build, forward the actual host selection to source fingerprinting: `RUSTFLOW_WORKSPACE_FEATURES=wasm,pyo3/extension-module` and `RUSTFLOW_WORKSPACE_NO_DEFAULT_FEATURES=1`. The host uses its existing Pyodide/Emscripten build pipeline and registers these bindings in its single Symbolica extension. This feature does not create a second numerical runtime or a standalone Python wheel.

Cache schema 8 encodes every numerical payload through Symbolica's existing binary `Encode`/`Decode`, including checkpoint coordinates, boundary values, error matrices, and singular-endpoint values and matching data. This avoids a demonstrated Astro serde roundtrip error without modifying Symbolica; the [minimal owner reproducer](../repros/symbolica-astro-float-serde/README.md) is included. The versioned envelope, source fingerprints, payload checksum, exact-consumption checks, and atomic replacement remain in force. Earlier cache schemas are rejected; supplied exact-ratio bundles provide the explicit portable import path.

The [feature milestone validation](../reports/validation/2026-10-06-browser-transport/report.json) passed 189 regular release tests with Malachite/Astro on Linux, including the complete planar 48-master transport, coherent amplitude assembly through HEPKit, rooted Python continuation and singular endpoints. A separate supplied-data gate checks all 4,360 complex starting coefficients and all 4,360 errors at their recorded precision through import, sixteen exact source hits and binary reload; it also checks metadata/provenance and rejects the earlier native binary bank. Native validation passed 261 regular tests plus the scientific bundle gate, and strict Clippy passed for both feature selections. These are host runtime checks of the portable backend. Actual Emscripten/Pyodide execution and full sixteen-configuration browser timing remain separate host acceptance gates.
