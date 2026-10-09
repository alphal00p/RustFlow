# Guarded native capacity: source-inspected design

Status: design only. No production padding, new runtime arity, or four-loop
numerical admission is implemented by this note. The massive sunset's requested
native closure now passes, but its full endpoint numerical comparison and the
mandatory four-loop evaluations remain incomplete.
The native dependency inspected is RustRed commit
`78969aab524b7d6a2eec36f59b01e9af1e04cc04`.

## Representation and current bottleneck

The physical weighted arity is the number of completed scalar/medium slots
plus two occupation slots per cut. For the supplied four-loop inputs the
completed basis has 14 slots, so one through four cuts require arities
16, 18, 20, 22. This count is a representation size, not a graph identity or
an admission criterion. The public assembly currently dispatches occupied
flows at exactly 7 and 9; its compact terminal and strict empty-support paths
run before that dispatch. See [assembly.rs](../src/finite_density/assembly.rs).

[WeightedMeasure](../src/finite_density/measure.rs) owns actual factors and
generates their distributional identities. It is currently const generic:
its constructor builds the numerator basis from its `[Atom; N]` factors and
source generation loops over all N factors. Adding constant or zero factors
here is not storage padding: it changes the physical basis and source corpus.
[geometry.rs](../src/finite_density/geometry.rs) also converts its factor
vector to an exact-size array. These owners must keep the original factors.

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

## RustFlow implementation choices and required probe

The smallest short-term dispatch change is to compile exact N=16,18,20,22
beside 7/9. It needs no padding and no new mathematical source rules, but adds
large reduction/flow monomorphizations and leaves runtime support tied to a
finite exact-size list. It does not admit any currently excluded contour.

For a capacity backend, separate dynamic physical source generation from
const-generic reduction storage. A vector-backed physical owner can preserve
the current array-taking public wrappers, emit physical shift/guard vectors,
then lower once to the selected backend capacity. A transitional alternative
is exact physical source generation at several small const arities followed
by one larger guarded reduction capacity; it reduces solver monomorphizations
but still needs physical source dispatch.

Specific frontend integration points are
[preparation.rs](../src/finite_density/preparation.rs) for source/target lowering,
[guarded.rs](../src/finite_density/guarded.rs) for checked embedding and cache
identity, and [reduction.rs](../src/finite_density/reduction.rs) for admitted
domains and final projection. Its `discovery_domains` broadens ordinary
nonpositive indices, so the existing intersection with the deformation's
admitted domain must retain `[0,0]` for all dummy axes. Its current
`native_integral` exports every N index and must instead project a checked
tail before `OccupiedFlowBoundary` receives the physical master labels.

A draft external probe at `/tmp/finite_density_capacity_probe.rs` derives an
unpadded occupation multiplication identity from `WeightedMeasure<2>`, then
embeds only its native labels/guards into capacity 4. It compares target
reduction, exact native replay and decoded-program application, tests a bulk
point where the surface identity is invalid, rejects escaped tails through a
frontend wrapper, and rejects a different embedding identity on cache load.
The draft has not yet been compiled or run; it is not acceptance evidence.
Follow-up tests must also compare a physical shell/occupation IBP corpus,
guard refinements, nonzero conditions, closure/exported basis and ordinary
nonpadded behavior before enabling any capacity dispatch.
