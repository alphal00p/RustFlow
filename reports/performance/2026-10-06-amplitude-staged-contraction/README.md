# Staged native amplitude contraction

The first-use gg→hg amplitude constructor takes **8.53–8.60 seconds instead of
33.00–33.08 seconds** on one CPU: about **74% less time, or 3.85–3.87× faster**.
All three final symbolic kernels are identical. Native observable values,
absolute errors, arithmetic-change bounds and binary precisions are also
identical to the frozen baseline.

The change runs Idenso's existing algebra simplifier before arithmetic expansion,
then retains expansion and the final native simplification. Both passes use the
same `Dots`/HEP settings. The first pass reduces connected tensor work while
preserving factorization; expansion then exposes deferred closed scopes for the
final pass. An explicit `contraction_complete()` check prevents an unresolved
closed tensor expression from reaching scalar extraction. No owner implementation,
coefficient algorithm, model, scientific input or precision setting changes.

| Native constructor measurement | Original (s) | Staged (s) | Reduction |
|---|---:|---:|---:|
| Initial frozen comparison | 33.0786 | 8.5978 | 74.01% |
| Fresh-process confirmation pair | 33.0010 | 8.5275 | 74.16% |

Each process is pinned to CPU 43 on a shared AMD EPYC 9754 host. Diagram
generation remains configured for one thread. Other compilation continued on
the host. Model loading and subsequent numerical evaluation are excluded from
constructor timing. Numerical evaluation takes 0.854 and 0.845 seconds in the
two staged runs. These measurements establish a constructor improvement; they
are not end-to-end notebook, physical-transport, cold-AMF or browser/WASM timings.
Two samples do not provide a statistical confidence interval.

## Scientific checks

The frozen original and staged kernels use the same effective UFO model and
reuse its native owner in every generated diagram. All three canonical kernel
strings have the same combined BLAKE3 hash,
`2a909a82184597995afe8d5774d4ca045e0949032b9c7e213575f354d4489410`.
The harness also requires exactly three reference expressions and verifies each
exact expanded difference is zero. Native coupling-order reconstruction and
form-factor-degree checks remain active and exact; statistical zero testing is
not substituted for them.

Numerical evaluation uses the unchanged W/Z form factors and admitted errors
from the earlier full sixteen-configuration native transport run. Their exact
dyadic values and precisions are restored without rounding. The model,
kinematics, masses and couplings match the native amplitude regression test;
requested accuracy remains 20 digits with 40 guard digits.

Both staged runs reproduce every original observable value, error and
arithmetic-change bound exactly, including precision. Exact rational comparison
also checks differences against the combined admitted errors and the requested
20-digit relative tolerance. The EW square, interference and HEFT square retain
**35, 36 and 47 verified relative digits**. The separate archived reference's
19-digit EW evidence limit remains unchanged.

Validation passed:

- The native amplitude integration test, including the archived reference,
  native model-owner identity and kernel reuse at two changed kinematic points.
- Exact three-kernel and observable/error comparisons in both measured runs.
- Strict release Clippy for all targets, with native features and again with
  `python_stubgen`; formatting checks.
- Independent read-only review of the production change by two agents.

## Rejected probes and diagnosis

Removing expansion entirely and moving expansion after a single contraction
were rejected. They reached scalar preparation quickly but stalled in native
polynomial conversion; their processes were stopped after 321 and 66 seconds,
respectively. Neither reached the kernel or observable equality gate.

The factorized first pass explicitly reports `Deferred`, with
`contraction_complete() == false`, despite `is_scalar() == true`. The delayed
single-pass intermediate retained closed tensor scopes and grew to 154.6 MB,
versus the original 5.6 MB. A second native reduction reports `Complete` and
produces the original intermediate byte for byte. Its SHA256 is
`56b032e0f554de42bf359e9e71d281cd1d86fd51b19353d1744cbc5c80cf9428`.
The failed probes' patches, timings, profiles and explicit rejection records
are retained in the evidence archive; they support no scientific equivalence
claim.

## Provenance and reproduction

The production change is based on commit `bf6b865`. Its source digest is
`b6e7b8db752f1f9d008094c1a80ca7233809de2be77cda4fafa24abee0cd3812`;
the original digest is
`a8e17584df583d823059db180f72a08fd18082d057578b16a3b3e4348bad76a2`.
Both use dependency digest
`9ee6eced2bc39b05bb333b1b021ab556d2a1bb790a255205d3f141fc64cc96ff`,
Symbolica/Numerica/Graphica `6defcca`, HEPKit `b96600b` and RustRed `7c1ed037`,
with the unchanged compiled arity registry 1–16. The frozen staged executable
SHA256 is `b753ab5cd3dd13b19fd0329f0b5c20471ec3096982ddaa1cc8c3e7cf9e9d1e32`.

`report.json` contains exact comparisons, both timing pairs, source/input/binary
hashes and validation results. `evidence.tar.gz` contains the frozen baseline
kernel strings and form-factor input document, harnesses, exact-rational checker,
diagnostic source, build/test/profile logs, patches and a SHA256 manifest.
Executables and raw perf samples remain outside the archive; their hashes or
textual profiles are included. Archive SHA256:
`ca8eeff1ef2e9fbb7a555db352dc34189d074ef5ab0d9f28b18cc4616bf43f39`.

Extract the archive in an isolated checkout containing this change. Copy
`candidate/amplitude_staged_check.rs` into `examples/` and use the same native
development environment and a valid Symbolica license:

```bash
cargo build --release --locked --example amplitude_staged_check
taskset -c 43 target/release/examples/amplitude_staged_check staged.json \
  baseline/baseline-kernel-expressions.json baseline/full16-high.json
python3 candidate/compare_observables.py \
  baseline/full16-high.json staged.json comparison.json
cargo test --release --locked --test gg_hg_amplitude
```

The input document contains results of the already verified physical transports
solely to compare amplitude evaluation. The production optimization does not
install a precomputed amplitude or bypass physical transport.
