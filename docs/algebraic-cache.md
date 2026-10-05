# Registered roots in the physical boundary cache

`RustFlow::new_algebraic` constructs `RustFlow<AlgebraicKinematicSystem>`. It uses
the same mutable `RustFlowCache` as the rational engine; existing rational
`RustFlow::new` calls retain their API. The algebraic `evaluate_to` additionally
requires a `RootGerm` for the destination:

```rust,ignore
let germ = RootGerm {
    sheets: BTreeMap::from([(root_symbol, RootSheet::Principal)]),
};
let source = CachedPoint::Exact(source_coordinates).with_root_germ(germ.clone())?;
// Insert a CachedBoundary carrying independently established coefficients and
// BoundaryAccuracy, together with flow.identity().clone() and this source point.
let result = flow.evaluate_to(
    &mut cache, &destination_coordinates, &germ, range,
    &options, &context, &admissible_cost_policy,
)?;
```

Each germ records one discrete sign relative to the principal square root at
that exact point. Root magnitudes are recomputed from the exact radicands at each
working precision. A germ supplies no new numerical accuracy. The point wrapper
preserves whether the coordinates originated as exact input, a derived exact
path image, or rounded numerical input. Root registries, the ordered basis,
normalization, prescription, domain and original nonzero conditions belong to
cache identity. Entries on opposite root sheets may coexist at identical
coordinates; exact-coordinate hits and candidate selection require the requested
germ.

This entrypoint accepts exact rational-complex physical coordinates and regular
affine paths. Every registered radicand must remain finite and nonzero along the
whole segment. For a radicand `R=N/D`, native Symbolica Gaussian polynomial
arithmetic and real/imaginary polynomial gcds identify its real-parameter zeros
and poles. Native real-root intervals and refinement check the complete parameter
interval. Unresolved certificates are rejected conservatively.

A germ can change along a regular segment even though the continued root is
analytic: the principal square-root convention jumps across the negative real
axis. RustFlow certifies these intersections using `N*conj(D)`. Odd crossings of
the negative axis change the endpoint germ; tangencies do not. One-sided endpoint
signs handle a segment starting or ending on the cut, where the principal root
is the positive imaginary root. Exact radicand substitution precedes numerical
seeding and germ classification, preserving the side under cancellation. This
check runs before source ranking, so a nearby source on an incompatible sheet
does not prevent reuse of a farther compatible source.

Native isolation/refinement calls are synchronous; cancellation is checked
around them and during bounded certificate refinement, but cannot interrupt an
individual native exact-algebra call. Root, matrix and original reduction
conditions are checked before candidates are ranked. Generic leading epsilon
coefficients are retained before point substitution, so a special point cannot
silently change the assumed Laurent structure, including on an exact cache hit.

Accepted intermediate endpoints retain their independently checked local root
germ and the exact coordinate image of the rounded path parameter. The caller's
existing `TransportCost` admissibility
contract still controls the integral's physical branch and logarithmic monodromy;
root signs alone do not determine those. A true radicand zero or pole on the
segment is rejected. General singular endpoints and sheet-specific removable
norm poles remain outside this entrypoint. Use `evaluate_prescribed_to` with a typed physical prescription identity and a
mandatory homotopy-admission callback for supported threshold detours; see
[Prescribed physical routes](prescribed-cache.md). Low-level explicit algebraic
contour transport remains available for other paths.

Rational and algebraic flows share candidate selection, independent
precision/order refinement, checkpoint comparison, inherited-error propagation,
stronger-evidence checks and transactional insertion. Algebraic amplification
uses the existing disk-subdivision/Gronwall estimator with rational coefficient
bounds and `|sqrt(R)| = sqrt(|R|)`, in the same weighted infinity norm as rational
transport. These are finite-precision uncertainty estimates, not directed
interval-arithmetic error certificates. The source's verified accuracy and
per-coefficient errors limit every new entry; repeating a low-accuracy source at
higher working precision does not promote its evidence.

Binary cache schema 5 stores the exact root registry, discrete germ and
explicit dense/canonical connection representation. Loading snapshots from
schemas 1 through 4 is an explicit incompatibility error. The existing source
and dependency fingerprints and payload integrity digest remain enforced;
older snapshots are never silently reinterpreted as algebraic data. Saving is
atomic and loading reconstructs the existing insertion index.

Focused regressions exercise progressive analytic transport, binary restart,
nearest compatible-source selection, opposite-sheet exact hits, inherited
accuracy refusal, original and generic-epsilon domain guards, imaginary roots,
complex mass transport, principal-cut crossings and tangencies, winding along
successive affine segments, exact endpoint cancellation, and the analytically
known weighted matrix norm. These branch tests do not authorize arbitrary
logarithmic monodromy: their source policies select the intended solution branch.

Canonical logarithmic forms use the same cache while keeping ordered letters
and constant matrices separate until a path is selected. See
[Canonical forms in the physical point cache](canonical-cache.md).
