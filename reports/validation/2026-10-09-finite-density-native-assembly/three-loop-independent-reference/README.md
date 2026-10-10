# Independent three-loop chain reference

`reference.json` contains independently derived raw Euclidean references for
both targets in `examples/finite_density/massless_three_loop_chain.json`.
It was saved before any successful native full-amplitude prediction for this
case. No supplied oracle answers or native AMF predictions are read.

The scalar and raised original mixed medium numerator are evaluated at
D=31/3, D=13/2 and D=15/4, plus Laurent orders [-3,-2,-1,0], at 50 and 80
working digits. The complete sector list is vacuum, cuts [0], [3], and [0,3].
The first three are zero under the independently justified dimensional
prescription. The double cut includes the original mass derivative and its
separate nonzero moving upper surface. The exact Gamma/Beta derivation is in
`docs/finite-density-three-loop-reference.md`; it is validation-only.

`validation.json` records 14 precision comparisons, with largest relative
change 4.52e-54. Direct Gauss-Legendre integration at D=31/3 uses orders
24,40,64,96 on the original radial/angular mass-jet kernels, retaining five
separate scalar/raised components. The largest final relative discrepancy
is 8.14e-15. These are observed checks, not rigorous interval enclosures.
The standalone run passed in 0.132 seconds with 15,428 KiB child peak RSS;
compilation, Nix startup and pre-launch executable hashing are excluded.

`generate.rs` imports only native arithmetic and elementary Gamma evaluation.
The numerical production evaluator must use native weighted AMF for every
occupied sector in this gate. The external reference generator is not a
subloop formula inside that evaluator. Source, input and executable hashes
are retained in `source-sha256.json` and `resources.provenance.json`.

The saved-output comparator is `tools/finite_density/compare_three_loop_reference.py`.
It requires four fixed-dimension or five Laurent profiles, all requested
cuts/targets/orders, exact input and normalization, 1e-12 relative agreement,
and 1e-25 absolute agreement for expected zeros. It also checks successful
native closure metadata and source/compile/executable identities. Current
Laurent harness records omit direct occupied-construction reports, so the
comparator requires a successful matching fixed-dimension report with the
same executable, source snapshot, native closure and source options. That
provenance limit is stated in the comparison report. Native comparisons
remain pending until complete predictions have been saved.

`exact-identities-validation.json` additionally records five exact Symbolica checks of the original mass jet, separate bulk and upper-surface Beta reductions, their complete rational ratio, and its epsilon expansion. This validation-only standalone executable imports symbolic arithmetic, without the finite-density evaluator or AMF. The saved numerical references were unchanged.
