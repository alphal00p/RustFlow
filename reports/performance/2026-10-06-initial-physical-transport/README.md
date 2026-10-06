# All sixteen initial physical transports and the native amplitude

Keeping the requested accuracy at **20 digits**, using 40 initial guard digits
and starting series order 64 reduces measured native physical transport time
from **270.213 to 163.129 seconds (39.63%, 1.656×)** across all sixteen gg→hg
configurations on one CPU. The comparison profile uses guard 60/order 96.
Both retain independent precision/order refinement, residual and conditioning
checks, supplied-error propagation and normal source-capped cache admission.

The same frozen executable evaluates both profiles. Its production source is
commit `9815c98` on base `699691c`, with source digest
`a8e17584df583d823059db180f72a08fd18082d057578b16a3b3e4348bad76a2`.
That commit reuses the final successful algebraic numerical object during
uncertainty admission; its separate matched comparison found no material
performance gain. The improvement measured here comes from the smaller starting
profile, not from that object-lifetime change.

| Stage | Guard 60/order 96 (s) | Guard 40/order 64 (s) |
|---|---:|---:|
| All 16 checked physical transports | 270.2131 | 163.1292 |
| W/Z form-factor contraction | 1.6879 | 1.6670 |
| Three coherent observable evaluations | 0.8733 | 0.8547 |
| Initial systems and supplied seeds | 0.6943 | 0.6840 |
| One form-factor projector construction | 1.3117 | 1.2903 |
| One native amplitude kernel construction | 33.6351 | 33.5720 |

Setup is measured separately. First-use native amplitude construction remains
substantial; it is not included in the transport figure. These are one run of
each profile on CPU 43 of a shared AMD EPYC 9754 host, not repeated statistical
trials. Other compilation and workloads continued. There is no cold AMF,
upstream implementation, browser/WASM or end-to-end notebook timing claim.

## Scientific gate

All sixteen exact native starting boundaries are preloaded into one shared,
growing cache. Two immutable canonical systems, one native projector and one
HEPKit amplitude kernel are reused. The original native seed values and errors
remain exact 415-bit dyadics with 40-digit supplied caps; their original identity
and provenance are retained as origin metadata. The importer authenticates the
fixed export hash and checks each current native configuration's dimension,
epsilon range, exact start and root germ before supplied-boundary admission.
It does not reuse a foreign binary cache or precomputed destination values.

The run covers W/Z, planar/nonplanar and all four crossings. In both profiles:

- All sixteen final boundaries retain 37 verified digits and a 40-digit input
  cap. Each solve has one accepted step and no rejections. Two verified
  boundaries are inserted per solve, growing the cache from 16 to 48 entries.
- The smaller profile finishes at 282 bits/order 96; the original finishes at
  349 bits/order 128. Requested accuracy remains 20, and the adaptive solver's
  refinement and acceptance requirements are unchanged.
- Exact rational comparisons check all **4,360 final complex master
  coefficients**, every progressive cache boundary, all eight form factors and
  all three final observables. Differences lie within the sum of the independently
  admitted error bounds. Both profiles' admitted errors and differences meet
  the requested tolerance; form factors and observables use relative tolerance.
  No binary64 rounding is used in these inequalities.
- Starting coordinates, branch germs, cache coordinates, source caps and
  provenance agree. Progressive checkpoints retain their own 37-digit cap.
- The EW square, interference and HEFT square report **35, 36 and 47 verified
  relative digits**, respectively, in both profiles. The harness fails unless
  every observable admits at least 20 digits.

Form-factor evaluation stays at 100 working decimal digits. Amplitude evaluation
stays at requested 20/guard 40 and uses the exact model, masses, kinematics and
couplings from the existing native amplitude test. Archived reference values
are used only after evaluation for comparison. Their separate EW input evidence
supports 19 relative digits; this run does not upgrade that reference claim.

## Evidence and reproduction

`report.json` contains the full comparison, per-case timings, scientific
diagnostics, source/dependency/input hashes and frozen executable hash.
`evidence.tar.gz` contains both complete exact output records, the standalone
harness, exact-rational comparison script, build and run logs, and a manifest.
The input is the unchanged seed export archived with the preceding
[two-case milestone](../2026-10-06-physical-compiled-reuse/README.md), SHA256
`984b64a1617eb111543fc3b97e1cd36c84c610cfd8a02b11e8dfa861809acee4`.

Extract the evidence, place `full_initial_transport_profile.rs` in an isolated
checkout's `examples` directory, and decompress the seed export to
`boundaries.json`. Use the same native development environment and valid
Symbolica license for both invocations:

```bash
cargo build --release --locked --example full_initial_transport_profile
INITIAL_TRANSPORT_GUARD=60 INITIAL_TRANSPORT_ORDER=96 taskset -c 43 \
  target/release/examples/full_initial_transport_profile boundaries.json high.json
INITIAL_TRANSPORT_GUARD=40 INITIAL_TRANSPORT_ORDER=64 taskset -c 43 \
  target/release/examples/full_initial_transport_profile boundaries.json low.json
python3 compare_full_profiles.py high.json low.json comparison.json
```

This report changes no library or notebook defaults. The measured physical
transport profile must be kept distinct from native automatic boundary-generation
options, which are outside this comparison.
