# Correlated Lorentz sources for finite-density reduction

Two exact source constructions improve selected small reductions, but neither
has closed the three-loop family in its bounded native pilot. No production
source, admission, ordering or numerical acceptance changed.

The [correlated boost](derivation.md) uses the same infinitesimal Lorentz
transformation on every integrated momentum. It annihilates all Gram products,
so physical scalar propagators and cut shells are invariant. Its divergence is
`(D−1) E_a`; occupation derivatives retain the upper and lower surface terms.
Scalar-factor conversion can lower cut powers and those terms remain present.
Independent exact checks cover 56 Cartesian vector fields, 340 Gram actions,
160 energy actions, 90 bulk moments and 720 upper-surface moments.

The [raw polynomial controls](native-source-probe/README.md) add four generated
rows to the original 52-row corpus. Three of the 19 reductions change: a raised
target shortens from 13 to 4 terms and an odd surface term from 2 to 1, while
another surface result grows from 1 to 8. The total changes from 192 to 189.
Prepending or appending the rows gives identical small-control results.
The [full pilot](native-active-pilot/README.md) reaches frontiers
44, 152, 467 and 1009 before its 300-second timeout during the following round.
It does not reach a final original-target and basis-derivative closure audit.

The [narrow normalized identity](normalized-polynomial-c1-proof.md) follows from
a separate tangent-shell derivation with an origin cutoff. It applies to exactly
one simple massless cut, zero lower occupation index, and zero energy-completion
indices; virtual factors and numerators may be Lorentz scalars in the admitted
family. It gives `(D−2) I_0 − mu I_1 = 0` and the corresponding positive-upper-
index recurrence. No reciprocal-energy integral is admitted. In particular,
this is not replay of an old source at a forbidden reciprocal-energy seed.
Thirty exact radial checks pass, and fourteen explicit counterexamples confirm
the need to exclude raised cuts, energy numerators and positive cut mass.

The [narrow-source controls](polynomial-normalized-source-probe/README.md) retain
all original physical targets. One scalar derivative shortens from 20 terms to
5, and the other 18 probes are unchanged, for a total of 177 terms. Both source
orders and both upper-index presentations agree; all native programs round-trip
and replay. Source-only tests exercise the surface recurrence and reject labels
outside its guards. The [full narrow-source pilot](polynomial-normalized-active-pilot/README.md)
still reaches the frontier limit, with 43, 190, 610 and 1025 unresolved labels.
Its final original output rows contain 10 and 1413 terms, versus 19 and 1506 in
the original-source baseline. Runtime was 228.999 seconds and peak child RSS
51,580 KiB; shared-host timings are not a controlled performance comparison.

All pilots use fresh native rules, the unchanged production integral order,
original output sums, explicit queried-point history and the same bounded
closure scheduler. Smaller intermediate expressions do not establish a closed
differential system, a boundary or an amplitude. Three-loop numerical validation
remains pending and continues to gate four-loop numerical scale-up.

The stored source, build, proof and resource records retain failed attempts and
the exact conditional domains. Verbose completed-round payloads are losslessly
gzip-archived with original/compressed hashes and byte-restoration checks.
`manifest.json` binds this complete frozen report. The subsequent combination
with polynomial temporal sources is a separate report.
