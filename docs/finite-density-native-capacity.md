# Guarded native storage capacity

Status: production storage adapter implemented and tested, including native
closure, integrated boundaries and endpoints at padded capacities. The implementation does not broaden contour, endpoint,
or distributional admission. The complete massive sunset has already passed
its fixed-dimension and five-profile Laurent reference comparisons. Mandatory
four-loop evaluations remain separate, incomplete gates.
The native dependency inspected is RustRed commit
`78969aab524b7d6a2eec36f59b01e9af1e04cc04`.

## Representation and current bottleneck

The physical weighted arity is the number of completed scalar/medium slots
plus two occupation slots per cut. For the supplied four-loop inputs the
completed basis has 14 slots, so one through four cuts require arities
16, 18, 20, 22. This count is a representation size, not a graph identity or
an admission criterion. The public assembly preserves exact capacities 7 and 9
and otherwise chooses the smallest available capacity from 12, 16, 20, 24, 32.
A physical arity above 32 requires a separately compiled const-generic Rust
capacity. Compact terminal and strict empty-support paths run before this
storage dispatch. See [assembly.rs](../src/finite_density/assembly.rs).

[WeightedMeasure](../src/finite_density/measure.rs) stores only the actual
factors in a vector. `from_physical` accepts these factors and their real roles;
the existing array-taking `new` remains an exact-arity compatibility wrapper.
Numerator-basis construction and factor/occupation differentiation iterate over
this physical vector. Native index and guard arrays alone have capacity N.
[geometry.rs](../src/finite_density/geometry.rs) passes its original factor
vector directly; no constant, zero, cut, occupation or ordinary factor is
invented for storage. Source identifiers and expressions at exact p=N are
unchanged.

## A sound embedding using the current public guarded API

For physical arity p and a compiled capacity N>=p, define
`embed(a)=(a_0,...,a_(p-1),0,...,0)`. First derive the p-coordinate sources
from the unmodified physical measure. Only then lower them to native arrays:

1. Append zero to every source shift and target; append `Ordinary` roles.
   Do not append a propagator, occupation, cut, multiplication identity, or IBP.
2. Append the exact bounds `[0,0]` to every source guard, discovery domain,
   admitted domain, measure-zero domain and subsequent refinement domain.
   Padding is not a vanishing integral: a nonzero tail is an invalid label.
3. Supply fresh unused dummy index symbols solely to satisfy the native
   coefficient-map arity. Reject dependence on these symbols in coefficients
   and all retained nonzero conditions, before and after reduction.
4. Reject any nonzero tail in terminals, application inputs, generated RHS
   terms, unresolved residuals, reconstructed target weights and exported
   basis labels. Strip a verified zero tail before physical boundary,
   normalization, differential-system and public result interfaces.
5. Bind p, N, the exact embedding and its version to the frontend measure/cache
   identity. Keep physical identities and zero certificates separately named;
   do not describe a dummy `Ordinary` role as a physical denominator.

This uses symbolic zero offsets in the original native source templates,
with their dummy indices constrained to zero by guards. It does **not** pass
numeric source coordinates to the existing guarded constructor, which rejects
them. Native `domain_case` specializes fixed bounds to numeric zero for search;
source image domains and exact replay enforce each source's fixed-zero guard.

The mathematical invariant is straightforward. Every admitted source
specialization has a zero dummy tail, and every source offset has zero dummy
tail, so every admitted row lies in the embedded physical integral space.
Exact linear elimination and source replay preserve that space. All native
ordering priorities are unchanged by adding identical trailing ordinary
zeros: denominator counts, sector comparisons, absolute/numerator degrees,
physical coordinate ties and occupation ties receive no new contribution.
Thus projecting a certified rule gives the same physical identity and strict
descent. This establishes soundness of the embedding, not equal heuristic
search output or equal performance at a given budget.

Relevant pinned native sources:

- [guarded/model.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/guarded/model.rs#L54):
  `GuardedSourceSystem::new` calls `SourceSystem::new`, with no prepared fixed
  coordinates or physical-arity parameter. `valid_indices` checks role bounds,
  so the frontend must additionally reject nonzero dummy tails.
- [guarded/search.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/guarded/search.rs#L28):
  source-domain image checks precede row instantiation; `domain_case` near
  line 346 turns a `[0,0]` bound into a numeric fixed case.
- [guarded/replay.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/guarded/replay.rs#L71):
  replay requires the order's arity to equal storage N, checks the fixed case,
  original-source guards, multiply-back, and role-aware descent.
- [index.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/index.rs#L417):
  the current order compares all N coordinates. The zero-tail argument above
  applies to this implementation; arbitrary coordinate permutations must
  retain the physical order and leave dummy coordinates trailing.
- [guarded/lifecycle.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/guarded/lifecycle.rs#L91):
  terminal admission and application enforce occupation/cut validity, but do
  not know the frontend's physical arity. Checking tails only on successful
  rule RHS terms is insufficient.

## Cost and the native capacity alternative

The guard-only embedding has a real search cost. A fixed numeric zero in a
coordinate case can still be seeded to a negative index within its nonpositive
sector. Such a dummy seed is rejected by every fixed-zero source guard, but is
still enumerated. [seed.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/seed.rs#L42)
and `SectorSolver::seed_frozen` in
[search.rs](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/search.rs#L227)
show that only prepared padding/removed coordinates are excluded from this
enumeration. Dummy axes can therefore increase work, change discovery
chronology, and reach native compact-power limits at large search depth.
Native powers remain in [-64,63]; padding does not widen that limit.

RustRed already has the lower-level pieces for efficient capacity support:
[SourceSystem::new_with_fixed and from_family_with_capacity](https://github.com/alphal00p/rustred/blob/78969aab524b7d6a2eec36f59b01e9af1e04cc04/crates/rustred-core/src/solver/source.rs#L45)
validate numeric fixed source coordinates, retain `active_arity`, and reject
coefficient dependence on them. Its private `resize` and public
`physical_coefficient` preserve/strip the exact coefficient map. A native
guarded capacity constructor could reuse those pieces and freeze dummy seeds.
It would also need:

- Physical-arity-aware order/replay and guarded persistence, replacing the
  current `physical_arity()==N` assumptions and binding active arity/fixed
  coordinates in saved context identity.
- Zero-tail validation at native domain, terminal and application boundaries,
  with no ordinary scaleless or sector-zero assumptions introduced.
- A capacity contract independent of the ordinary application registry.
  `fits_storage` currently admits only exact arity unless the optional
  `capacity-dispatch` feature selects the matching configured capacity;
  this RustFlow dependency does not enable that feature. Merely setting an
  order's physical arity is therefore not a general arbitrary-capacity API.

That owner-level API extension is cleaner and avoids wasted dummy seeds, but
it requires a reviewed dependency change and updated native transport tests.
The guarded-only embedding can be an interim sound adapter with its resource
cost stated explicitly.

## RustFlow implementation and validation

`WeightedMeasure<N>` keeps p physical factors and emits native index/guard
arrays of capacity N. Every emitted source has zero tail shifts and fixed-zero
tail guards. Source vector fields cannot depend on unused dummy index symbols.
Normal-source shifting and domain pullbacks act only on actual physical slots.

[preparation.rs](../src/finite_density/preparation.rs) pads already-converted
physical target labels and fixes every admitted and support-zero tail domain.
[guarded.rs](../src/finite_density/guarded.rs) validates the embedding at source,
program discovery/decoding/application and terminal boundaries, and records
`zero-tail-storage-v1:physical=p:capacity=N` only when p<N. Exact p=N measure
identity is preserved. [reduction.rs](../src/finite_density/reduction.rs) requires
the same true arity in its deformation, retains zero tails during discovery
and refinement, and verifies then trims every exported physical master, target
weight and candidate reduction. Physical region, boundary, normalization and
transport owners receive unpadded labels.

The adapter retains the seed-enumeration cost described above. Compiled bucket
bounds and bounded search budgets are backend resource limits, not loop-order
or graph-name assumptions. Native exact source replay and closure are still
required before numerical evaluation; the wider storage dispatch does not
promise closure for every represented family. Massless flowing graphs remain
unadmitted even if their index arrays fit a bucket.

The fresh finite-density library suite passed 36/36, including equivalence of
all five physical source presentations at p=4 versus capacity 6 and invalid
physical-factor construction. The guarded adapter suite passed 3/3, covering
invalid tails/domains/coefficients, exact replay and decoding, and identical
physical closure exports. The physical occupied-flow test passed p=7 versus
capacity 12 and p=9 versus capacity 12: the closed bases have 6 and 11 members,
respectively, and all four target values agree at the test's 15-digit criterion.
The saved physical bases, connection matrices, target weights and retained
nonzero conditions are also exactly identical; see
[`native-capacity-physical-comparison.json`](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-physical-comparison.json).
That prebuilt test process took 245.7088 seconds and peak RSS 55332 KiB; native
storage padding can change search cost even when results agree. Another 33
boundary, terminal, input, normalization, measure and independent-chemical-
potential regression checks passed.

Evidence is retained in the native-assembly report directory as
[`native-capacity-unit-resources.json`](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-unit-resources.json),
[`native-capacity-adapter-resources.json`](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-adapter-resources.json),
[`native-capacity-physical-resources.json`](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-physical-resources.json),
and [`native-capacity-fast-gates.json`](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-fast-gates.json),
tied to `native-capacity-source-hashes.json`. These gates compare storage
representations of the admitted massive problem; they contain no numerical
four-loop prediction or supplied-oracle comparison.

The later final-source owner audit also exercised genuine four-loop sources
for the E7 definition at physical arities 16 and 18, using capacities 16 and 20.
Its depth-3 single-cut pilot derived 84 source rows and 1650 native rules,
with 5 actual rule applications, but reached a frontier of 147 above its budget
128. Its double-cut pilot derived 328 source rows and 874 native rules, with
one application, and remained unresolved after three rounds at frontier 88.
Closure discovery took 21.611 and 108.605 seconds respectively. These are
bounded algebraic pilots, with retained native proof gaps and zero contour or
numerical admission. See
[`native-capacity-final-owner-gates.json`](../reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-final-owner-gates.json),
which preserves the lower-depth attempts and final-source snapshot as well.

The original design-stage probe below is historical evidence, with its source
snapshot preserved separately from the current implementation.

A probe in [finite_density_capacity.rs](../tests/finite_density_capacity.rs) derives an
unpadded occupation multiplication identity from `WeightedMeasure<2>`, then
embeds only its native labels/guards into capacity 4. It compares target
reduction, exact native replay and decoded-program application, tests a bulk
point where the surface identity is invalid, rejects escaped tails through a
frontend wrapper, and rejects a different embedding identity on cache load.
The fresh native test passed 1/1: direct execution took 0.06316 seconds with
peak child RSS 12388 KiB, excluding compilation. The exact command, binary
hash and source snapshot are recorded in
[`frontier-sector-capacity-gates.json`](../reports/validation/2026-10-09-finite-density-native-assembly/frontier-sector-capacity-gates.json).
This verifies the two-to-four coordinate algebraic probe only; no production
capacity dispatch or numerical amplitude was tested.
The current tests above extend this historical probe through the production
physical owner, guarded closure and actual occupied numerical evaluation.
