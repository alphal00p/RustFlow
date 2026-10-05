# Pinned DiffExp benchmark inputs

The [ancillary provenance report](../reports/diffexp/ancillary-provenance.json)
records archive and member SHA256 hashes, matrix dimensions, reference precision,
and comparisons between paper versions. The inputs belong to DiffExp commit
`784c8229bf92369a03f011a48e161522c8c54bbd`'s example notebooks. The [full 13-master planar one-loop benchmark](fivepoint-planar.md) has
completed live original and native PH1→PH6 transport through ε⁴, including
precision/order refinement. The remaining larger five-point benchmarks are
prepared but have not completed numerical acceptance.

Reproduce the downloads and verify the extracted scientific files with:

```sh
python scripts/fetch_diffexp_ancillaries.py
```

The default output is `target/diffexp-oracle/ancillary`; `--output PATH` selects
another directory. The script verifies existing files, refuses changed bytes,
and leaves the original checkouts untouched. It executes no Mathematica code.

| Example | Input version | System | Reference points |
| --- | --- | --- | --- |
| Massless nonplanar five-point | [1812.11160v2](https://arxiv.org/abs/1812.11160v2) | 108 masters, 31 letters | X0 and X1, epsilon orders 0–4 |
| Planar one-mass five-point | [2005.04195v2](https://arxiv.org/abs/2005.04195v2) | 13, 75, 74, and 86 masters; 58 letters | Seven points per family, epsilon orders 0–4 |

The nonplanar notebook explicitly selects v2, even though a later paper version
exists. Its supplied boundary normalization includes
`epsilon^4 exp(2 EulerGamma epsilon)`. The stored X0 numbers have approximately
65–67 precision digits and X1 approximately 50–58; this does not independently
certify all those digits.

The planar notebook prints a v1 URL but downloads an unversioned URL. The
unversioned source archive currently matches v2 byte for byte. Both versions are
preserved for comparison: their alphabet and differential-equation matrices are
identical, while the reference files and `determinants.m` differ. Matching
reference rows agree within their recorded accuracy, around 132 absolute digits.
The fixture exports retain each point's supplied error estimate and distinguish
finite-accuracy numerical zeros from exact values.

**Known upstream data defect:** v2's `mzz` point 7 has only 13 entries in its
epsilon-zero row, although this family has 74 masters. Its other four rows have
74 entries; v1 has 74 in every row. The exported metadata marks this v2 point
unusable as a complete boundary. It is not repaired or replaced silently.
The notebook's PH1→PH6 and PH1→PH2 benchmarks use valid points and are unaffected.

The target directory also contains sparse planar coefficient tensors, the
alphabet with explicit algebraic-root definitions, and nonplanar matrix/reference
JSON exports. They preserve Wolfram precision metadata and require explicit
branch handling before native physical transport. Only scientific data was loaded
for these exports; ancillary helper implementations were not evaluated or copied
into Rust.

The [equal-mass banana](banana-equal.md) now has a complete four-master
infinity-boundary acceptance test, and the [unequal-mass banana](banana-unequal.md)
validates all 15 masters through ε⁴ on two native paths against a complete
original-DiffExp endpoint oracle. Both use analytic gamma-function boundaries;
no numerical reference seeds the integration. The unequal original full route
was separately censored by its memory cap, while its alternate route completed.
The separate [planar one-loop five-point result](fivepoint-planar.md) validates
all 65 coefficients using its supplied ancillary boundary. The 75-, 74-, and
86-master planar families and the 108-master nonplanar system remain pending.

## Complete 74-master MZZ family

The full MZZ system from the pinned 2005.04195v2 ancillary data now passes supplied-boundary PH1→PH6 transport through epsilon 4. All 370 coefficients agree at 20 absolute digits between independent 60-digit/order-56 and 70-digit/order-64 native runs. The maximum refinement difference is 1.9624e−48; the higher run differs from the ancillary endpoint by 2.4385e−58 and the live original endpoint by 3.0959e−26. These are observed residuals, not claims of 58 guaranteed digits.

The exact canonical letters and constants are pulled back before materializing physical derivative matrices. Original source-domain conditions and all three root sheets are retained. Both native profiles start independently from the supplied boundary. The fixture records finite-accuracy zeros, component Accuracy/Precision and inherited source-error metadata. This low-level benchmark does not certify propagation of an arbitrary input-error box, generate its own boundary, or evaluate a scattering amplitude.

Original transport took 1147.09 s (working precision 150, initial order 40, goal 15, 137 segments, CPU 28). Native transport took 3524.86 s at 60/56 (1191 accepted, 1533 rejected) and 4541.57 s at 70/64 (1214 accepted, 1591 rejected), on CPU 29. These single-run timings perform different work and establish no performance ratio. The full raw outputs, declarations, caps, precision metadata and hashes are in [the validation report](../reports/diffexp/fivepoint-mzz-validation.json).

The offline regression runs both native profiles with `cargo test --release --test fivepoint_mzz -- --ignored --nocapture` and needs no Mathematica. Its source passes strict Clippy compilation; production-library numerical execution remains pending; the completed evidence uses immutable library SHA `5be94f7cf6d87188b8fc1e64831126a34849e6b14ccfc499bb936981cf369754`.

The shared canonical helper also passes a [current-library full13 check](../reports/validation/2026-10-05-fivepoint-canonical-helper.json): all 65 coefficients agree with both references at 20 absolute digits, and both endpoint root signs match. This checks its aggregate waypoint budget and root-state handoff without claiming a new full74 run.
