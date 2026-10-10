# Explicit measure-zero domains during native source discovery

Native guarded discovery can now remove an instantiated source column when its
whole integer-domain image is contained in a caller-supplied measure-zero box.
Previously these boxes skipped zero search targets and pruned zero terms during
concrete rule application, but their columns could still enter discovery pivots
and exact row reduction. That difference could hide a direct zero identity or
produce a less useful recurrence. This change does not infer new zero sectors.

The weighted measure owner remains responsible for each supplied zero theorem.
Examples include the existing certified free virtual polynomial directions and
lower energy contacts with their explicit support/origin prescriptions. Ordinary
vacuum scalelessness is not applied to occupied integrals. An occupation index
zero means a theta factor; it is never discarded merely because its index is
zero.

The projection uses the complete recorded discovery box and exact source shifts.
A partial intersection with a zero box is insufficient. Every surviving
nonzero instantiated term must first have a valid occupation-index image; a
zero projection cannot hide an invalid image. The original source guard and all
replayed nonzero conditions remain part of the proof. There is no assumption
that zero boxes are monotone or closed under arbitrary translations.

Exact rule replay invokes the same projection on the same original rows and
discovery domain. When a pivot is recentered, the final rule domain retains its
intersection with the translated pullback of the discovery domain. Thus a zero
identity proved on one index face cannot silently extend to another face. A
one-original-row direct-zero search can now succeed after uniformly certified
zero columns are removed; it still performs no row elimination of its own.

This is a native guarded proof change with a deliberate persistence boundary:
`rustred.guarded-source-program.v2` replaces version 1. Version-1 guarded proof
programs are explicitly rejected with a rediscovery diagnostic. Their saved
RHS columns and extra pivot conditions may differ from the new quotient proof,
so the loader does not silently normalize or migrate them. Ordinary native
program persistence is unchanged. RustFlow's guarded source identities also
bind the native source digest, so its old guarded caches are invalidated.

Historical validation archives remain intact. A diagnostic may explicitly
extract their original source equations, domains, conditions and zero evidence,
verify the archived hashes, and run fresh discovery using the current native
library. Such a run must record the new library hash and must not import the old
rules as version-2 proofs.

The focused regressions cover pointwise zero discovery, nonzero-condition
retention, persistence roundtrip and context mismatch, partial-box rejection,
occupation bulk preservation, invalid occupation images, recentering guards,
and coexistence with required-cut zeros. Test and numerical outcomes are
recorded separately in validation reports; this implementation does not by
itself establish closure or numerical acceptance for a new graph.
