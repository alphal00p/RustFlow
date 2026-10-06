# Explicit physical propagator slots in native interfaces

This milestone exposes the generic partial-cut implementation from `fa26f780712c154c1fb4fb273e0bda84a954e842` through the existing `EvaluationOptions` and CLI settings. `MassMode::from_selection` owns the shared name/list decoding. Native family and cut owners retain the physical slot, ISP, cut-line and supported-measure checks. No numerical algorithm, persistent cache schema, graph parser, model registry or dependency revision was added.

`deformed_propagator_slots=[2]` denotes zero-based physical slot 2. The native connected-cut fixture has propagator edges `[EdgeId(1), EdgeId(2), EdgeId(3), EdgeId(4)]`, so it selects native edge 3. Cut conversion keeps these positions; ordinary partial-fraction families have their own denominator inventories. The documentation explains that distinction. Explicit selections can be roundtripped with `mass_mode="explicit"`; conflicting named modes are rejected, while duplicates retain the native mask owner's existing treatment.

The final release gate passed 33 tests with `python_stubgen` enabled:

- 29 across `massless_phase_space`, `mixed_cut_flow`, `native_cut_interfaces`, `partial_cut_boundaries`, and `python_cut_interfaces`.
- Three exact ordinary placement tests in `cache`, including option roundtripping and ISP rejection.
- One ordinary Euclidean bubble comparison across alternative placements in `end_to_end`.

The Python and native DOT CLI tests compare partial and all-uncut evaluations for a connected two-loop cut diagram. They retain exact native expressions and arbitrary-precision outputs, check conflicts and cut/ISP admission, and exercise generated stubs, progress and GIL-free cancellation. Explicit selection is also checked before terminal dispatch and when evaluating an already prepared terminal, including a zero projection. Samples continue to report working precision only; Laurent reconstruction retains independent-fit evidence. The existing three-loop nonfactorized test independently checks 20 relative digits and higher-precision stability; its original upstream comparison remains in the preceding milestone's archive.

Commands:

```text
cargo test --locked --release --features python_stubgen --test python_cut_interfaces --test native_cut_interfaces --test partial_cut_boundaries --test mixed_cut_flow --test massless_phase_space -- --nocapture --test-threads=1
cargo test --locked --release --features python_stubgen --test cache mass_ -- --nocapture --test-threads=1
cargo test --locked --release --features python_stubgen --test end_to_end alternative_mass_placements -- --nocapture --test-threads=1
cargo clippy --locked --all-targets --features python_stubgen -- -D warnings
cargo fmt --all -- --check
git diff --check
```

`validation.json` records the final command results, compiled fingerprints and frozen artifact hashes. `dependency-ownership.json` confirms single ownership for Symbolica, Numerica, RustRed, HEPKit, Linnet, Spenso and Idenso. `compiled-dependencies.json.gz` retains the resolved feature graph. Source hashes were frozen before the final build and checked after validation.

The first strict Clippy pass found an unused connected-diagram helper in another test that imports the simple graph fixture. The connected fixture was moved to its own support module, imported only by the two interface suites. Those 12 tests were rerun after this test-only change; the other 21 release tests had already passed against the unchanged final library. The initial Clippy log is preserved separately.

The initial feature run passed 20 tests, but its build script had captured the source fingerprint before the final terminal evaluation admission guard was added during dependency compilation. It is retained separately as `initial-feature-tests.log.gz`; the final run rebuilt the library with the completed frozen source. No initial-run binary is the final validation artifact.

The earlier measured `fa26f78` binaries remain frozen under `/common/dev/amflow/target/partial-cut-release-fa26f78/`. The tested interface library, CLI and Python/native interface executables are separately preserved under `/common/dev/amflow/target/explicit-cut-slots-release-a6e7f926-final/`. Their hashes are included in the report. Timing lines in the regression logs are diagnostics, not a new performance-parity claim.

Automatic scope remains the complete positive-energy final-state cut classes documented in `docs/cuts.md`, with uniform +i0 virtual directions and certified real-domain poles. Dependent soft-pole partial fractions, incomplete final states, mixed-sign pinches and general massive many-body terminals remain explicitly unsupported. Supplied cut tables require a backend that owns `reduce_cut`; ordinary reduction tables do not authorize this measure.
