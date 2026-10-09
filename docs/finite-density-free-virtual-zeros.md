# Measure-aware free-virtual-loop zero domains

The source option `WeightedSourceOptions::free_virtual_zero_sectors` enables
these domains and defaults to false. Its factory requires the matching sealed
massless endpoint/origin evidence; setting the option without that proof is an
error. Existing source behavior is preserved when the option is false. Closure
and runtime comparisons remain separate from this exact zero proof.

## Physical zero and its limits

Use the occupied family's compact-first coordinates and only the virtual columns
of each uncut physical momentum row R_e. Required shells and every upper/lower
occupation depend solely on the compact coordinates. Under the sealed massless
permit, completions have nonpositive indices, so they and all nonpositive
physical factors are finite polynomials in the virtual momenta.

If the positive ordinary factors span a proper virtual row space W, choose a
nonzero exact rational v in its nullspace. An invertible rational change of the
virtual loop basis makes one virtual vector t run along v. Every positive
ordinary denominator is independent of t. The complete remaining t dependence
is polynomial, even when it contains eta, compact vectors, or other virtual
vectors as polynomial coefficients. Every tensor moment of the unrestricted
integral over t is dimensionally scaleless and vanishes. This proof acts only
on unrestricted virtual directions; a compact loop cannot be used in the rank
matrix or declared scaleless.

The statement holds at fixed positive eta before the physical endpoint. For any
finite cut/occupation order, use the same high-D origin prescription to pair the
remaining compact distributions and continue the zero meromorphically. Negative
powers do not introduce virtual poles, and no evaluated virtual period is used.
The proof is compatible with the existing active-rank label classification.
It is distinct from lower-contact zero jets and from an ordinary RustRed sector
census applied to a weighted measure.

## Finite exact boxes covering all rank-deficient positive supports

The ground set is the list of uncut physical virtual rows R_e, of full total
rank h. For every independent (h-1)-row subset, close its span under all ground
rows. Deduplicate the resulting rank-(h-1) flats F. Equivalently, retain a
exact rational null vector v_F normalized to first nonzero entry one and define F={e:R_e v_F=0}.

For each flat attach the index box:

- every uncut physical slot outside F has n_e<=0;
- physical slots inside F retain the original admitted bounds;
- every required cut has n_cut>=1 (the existing required-cut zero covers n<=0);
- occupation indices retain their nonnegative admitted bounds;
- every completion retains n<=0;
- every native storage-tail coordinate remains fixed to zero.

Each box is sufficient because every allowed positive row lies in F. The union
is also complete for rank-deficient supports: any independent set of rank<h
extends, within the full-rank ground matroid, to rank h-1, and hence lies in one
of these flats. Zero virtual rows lie in every flat and never force a physical
slot nonpositive. For h=1 the unique flat is rank zero; for h=0 there are no
free-virtual-loop zero domains.

Domain construction must intersect the already admitted box, never construct a
looser replacement. Invalid completion, occupation, cut, or storage coordinates
must fail admission before these physical zeros can apply. Source identities
and persisted replay must bind the actual full routed row list, canonical flat
closures, normalized rational null witnesses, deformation, family, and theorem version.
No new symbolic nonzero condition arises from these exact assigned rational
rows; any pre-existing source/target conditions are retained normally.

## Bounded reusable algorithm

The sealed helper
`free_virtual_zero_supports(family,shifted,budget)` returns
private-constructed certificates with allowed-positive slots, forced-nonpositive
slots and a normalized rational virtual null vector. It validates family/deformation
binding before exposing any support certificate.

Compute binomial(number_of_uncut_rows,h-1) with checked arithmetic before
enumeration; stop at an explicit budget. For the independent-block variant,
the h coordinate hyperplanes give exactly h flats directly, avoiding subset
enumeration. Canonicalize closures by physical-slot order. Verify every retained
null vector against all ground rows and verify the flat rank h-1.

The source factory can translate these certificates into the boxes above using
its existing admitted bounds. The same explicit budget bounds candidate subsets and retained domains.
Exhaustion returns a resource error and supplies no zero domains; the
implementation does not claim partial coverage. The option and every canonical
flat/null witness enter the persisted source identity. The boundary identity
also retains the selected source options through its family-bound proof.

## Independent E7 enumeration check

The independent prototype reads only the input graph definition,
constructs the compact-first rational routing, enumerates primitive integer null-vector
flats, and independently exhausts every positive active support. The Python witness
uses primitive integers; the Rust helper normalizes the same rational direction
to first nonzero coefficient one. These normalizations define identical flats. It checks that
box membership is equivalent to rank<h for all 160 supports across the three
occupied cut sets. Its exact enumeration gives:

- cut [0]: 6 flats among 64 supports;
- cut [4]: 6 flats among 64 supports;
- cuts [0,4]: 2 flats among 32 supports.

For the double cut, the two flats are exactly the two cases where one entire
virtual quadratic block has nonpositive indices. The pure transfer slot stays
unrestricted in both. These are source-support algebra checks only, not a
closure or numerical-flow result.

The Rust certificate unit test independently enumerates all active supports and
checks equivalence with the flat boxes. The guarded-source regression compares
policy off/on, persisted replay identities, restored full-rank nonzeros, invalid
inverse completions, negative occupation indices, storage-tail escape and budget
failure. Their pass status is recorded only after the coordinated native build.
