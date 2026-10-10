# Native moving-shell one-loop experiment

This is an exploratory exact-source experiment, not a production deformation or a multiloop acceptance result. The original scalar coordinates are `g=q²` and `u=q0`. Set `A=m²+eta>0`, `B=mu²−m²>0`, and use the measure `C_n(g−A) u^(−p) H_h(B−u²+g) H_l(u)`, with polynomial energy insertions `p<=0`. The upper spatial ball stays fixed when eta changes. Positive-A real support makes all positive lower-contact labels vanish; that statement is bound into the source context, separately from the generated polynomial identities.

Four polynomial vector fields generate the source corpus through the existing `WeightedMeasure` API: spatial Euler, temporal translation, full Euler, and temporal Euler `u*U`. Multiplication sources are included. Discovery and all reductions use native RustRed, including native ordering and exact replay. There is no evaluated period, hand-written reduction, chosen terminal list, or alternate elimination backend.

The automatic frontier searches both its labels and their auxiliary derivatives. After each discovery it reduces the old frontier before forming new derivative obligations. Only after the frontier stabilizes does it declare those labels as explicit terminals, then replay the original targets and every basis derivative.

The resulting two labels are `I=C1 H0` and `S=C1 H1`, with exact system

```text
d_eta I = (D−2)/(2A) I − (B/A) S
d_eta S = −S/[2(A+B)]
```

The original `(g,u)` scalar run uses 14 generated sources and closes in three rounds with seven rules. The broader fixed-original-numerator run starts from the bulk, u² bulk, and u² raised-cut targets, closes in three rounds with 15 rules, and still discovers precisely the same two labels. It checks native rational coefficients exactly for `N=g+u²`:

```text
N*C1 = A(D+1)/D I + 2B(A+B)/D S
N*(C2 H0 − C1 H1) = (D+1)/2 I − 2(A+B) S
```

The second identity is the physical mass derivative, `partial_A−partial_B`, with the original polynomial held fixed. Auxiliary differentiation is only `partial_A`. The upper-contact term is retained. Native conditions on A, A+B, and any target-reconstruction denominators are retained; see the exact program and condition audit.

The independent spatial-coordinate control uses `(g,s=|r|²)`, seven generated sources, and closes to the same matrix in four rounds with nine rules. That representation is restricted to scalar/even-energy polynomials; it is not a generic completion for odd energy insertions or products of energies from different loops.

Two earlier bounded diagnostics remain preserved: derivative-only search did not close after eight rounds, and adding frontier searches without retiring its reduced labels still did not close after eight rounds. Their exit status zero means the diagnostic completed and wrote an explicit null `closed` result, not that closure succeeded. The succeeding change combines temporal-Euler source presentation with retirement of reduced frontier labels; these two changes have not been separately benchmarked.

Files with `resources` in their names record only the prebuilt process, excluding compilation and Nix startup. `source-binding.json` identifies the canonical dependency libraries and source snapshot. The `full-targets` and final `roundtrip` directories preserve native programs, exact coefficient checks, and final replay evidence. Separate parent-directory radial integration and shared RustFlow transport checks are independent numerical evidence; this directory itself makes only exact-source and closure claims.

For an occupied one-loop tadpole the existing fixed-shell deformation has no uncut denominator and is eta independent. Therefore this experiment establishes a valid nontrivial moving-shell flow at one loop; it does not demonstrate that moving shells improve multiloop reduction. No production admission, physical-mass convention, or generic numerator policy was changed.
