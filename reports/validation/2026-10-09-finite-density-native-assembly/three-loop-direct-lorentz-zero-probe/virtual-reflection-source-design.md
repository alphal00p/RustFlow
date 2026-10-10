# Bounded virtual-reflection source proposal (not implemented)

This is a source-extension proposal, not an evaluated-period shortcut. All
relations must enter the native guarded corpus and pass native discovery and
original-source replay. Existing source/contour admission remains authoritative.

## Exact change of variables

Choose an unoccupied loop `k_v`, and a rational linear combination `b` of the
other loops. Set `k_v -> b-k_v`; leave every other loop fixed. Its full routing
matrix squares to the identity and has determinant `-1`. In occupied-first
coordinates no required shell or upper/lower occupation factor depends on
`k_v`, so each distribution factor is preserved individually. Polynomial
medium insertions transform with the energy coordinate, not as scalar labels.

Use the HEPKit/native exact scalar-coordinate map to transform every full
deformed factor, including physical masses and eta. Verify the factor identities
exactly. The real virtual change of variables preserves the full Feynman
integration domain at finite eta; its use still requires the caller's admitted
common distribution/contour prescription. It does not independently admit an
unsafe occupied family or assign a value to a singular product.

Generic finite candidates are `b=0`, the reflection center `b=-2r/a` for a
physical momentum `a*k_v+r`, and pair centers
`b=-r_i/a_i-r_j/a_j`. Enumerating candidates never proves a symmetry: every
full factor and distribution must pass the exact map check. Coefficient,
candidate, polynomial-degree and expanded-term budgets must be checked before
allocation, and insufficient budgets produce explicit unsupported/limit output.

## Native index representation constrains the first implementation

Current guarded source terms encode `a_i + constant`, not an arbitrary index
permutation. It is incorrect to represent a denominator exchange as fixed shifts
while both exchanged powers remain independent symbolic indices.

A small valid initial option therefore retains arbitrary powers only for
pointwise-invariant denominator factors. Every non-invariant factor is either
fixed to a nonpositive integer (and expanded polynomially), or is assigned fixed
positive powers in a bounded source slice whose transformed denominators occur
in the same full factor list. Every non-invariant numerator-completion index is
fixed to its bounded polynomial degree; invariant indices, cut indices, and
occupation indices remain symbolic on their existing admitted domains.

For a fixed slice, form the exact polynomial numerator image using the existing
inverse-propagator basis. Convert every image term to a fixed native shift and
emit `I - transformed(I) = 0`. The guard must fix all coordinates whose concrete
values were used in expansion/permutation. Require completion powers <=0, tails
exactly zero, original C/H role domains, and all raw coefficient conditions.
The emitted row and its exact map/guard are bound into the measure identity.

This includes odd moments under `k_v -> -k_v` without parity assumptions on
unbounded symbolic numerator indices: e.g. a fixed degree-one numerator gives
the native source `2 I=0`, valid for every invariant propagator/cut/H power.
It also includes an external-minus-loop reflection on a finite denominator-power
slice; it does not claim a generic arbitrary-power denominator permutation.

## Corpus lifecycle

Generate a complete bounded collection before constructing the immutable native
context. Adding sources later changes persistence/replay identity. An adaptive
collection would need an explicit native source-extension/replay API; the
existing exact-corpus `union_replayed` cannot silently accept it. The simplest
controlled experiment is an upfront finite degree/power budget, false by default,
with source ordering and all proof metadata in the cache identity.

## Required tests and honest performance experiment

1. Exact Cartesian/symbolic substitution reproduces each emitted polynomial row;
   shell and H factors are unchanged, determinant absolute value is one.
2. Odd virtual numerator with one invariant denominator produces a native replayed
   zero; even numerator does not; an occupied-loop reflection is rejected.
3. Deformed/mass-mismatched factor exchanges reject; unequal symbolic swapped
   powers reject unless fixed by the guard; inverse completion powers reject.
4. All changed indices fixed by guards, all padding zero, raw conditions retained,
   explicit degree/term/source limits, and cache mismatch when maps change.
5. Bounded native off/on runs use identical targets/closure budgets and report
   actual final derivative closure, historical gaps and application failures.
   An earlier direct hit alone is not evidence of global improvement.

## Existing-source diagnostic comes first

The three-loop odd labels already have a direct Lorentz source: in occupied-first
coordinates `lorentz/2/0` is `P1 . d/dP2`. At
`[1,1,2,1,0,0,...]`, all terms except
`-2 I[1,2,2,1,0,-1,...]` should vanish. Thus a new reflection source is not required
to make this particular relation true. First specialize the actual serialized
row, then compare native discovery with the unchanged corpus, a stable row-order
change, and a diagnostic corpus containing only that existing row. Every result
must pass existing native replay; no bespoke elimination or analytic zero is
installed. This distinguishes missing identity from first-rule/seed selection.
