# Independent audit: verified guarded-program reuse

Status: read-only design review of the current production owners, 2026-10-10. No prototype or performance result is claimed. No production files were edited.

## Verdict and invariant

No soundness blocker was found in the proposed narrow design. A live `GuardedProgram` can serve as evidence that every contained rule has already passed exact original-source replay, provided its fields become private to `guarded::lifecycle`, no mutable getter or unchecked constructor is added, and all arbitrary-rule/deserialization entry points still replay. The invariant is conditional algebraic consequence under the supplied measure/source/zero-domain assumptions. It does not certify the physical source identities, complete domain coverage, a master basis, closure, or numerical admission.

Construction establishes the invariant once. Verified terminal replacement changes no rules or source context. Verified union concatenates two already-valid rule sequences under an exactly equal context, so each retained rule keeps the same proof meaning. These operations must still check the structural compatibility conditions below. No stored Boolean or digest should be accepted as an alternative to replay at the persistence boundary.

## Inspected owners and paths

* `vendor/rustred/crates/rustred-core/src/solver/guarded/lifecycle.rs`: `GuardedProgram::new` currently checks common coordinate permutation, disallows order programs, calls `sources.replay_rule` on every rule, and checks all replacement terminals. `with_terminals_replayed` and `union_replayed` ultimately call `new`. Those replayed APIs and their corruption regressions should stay unchanged.
* `guarded/persistence.rs`: `decode_generated` binds the full expected source context, decodes candidate/source/domain/condition/order records, and finishes with `Self::new`. It must keep that final replay. `encode_native` needs only immutable source/rule/terminal getters after fields become private. Its format is not itself an authority boundary.
* `guarded/model.rs`: discovery outputs owned `GuardedRule` values in `GuardedSolution`; these remain arbitrary inputs until the program constructor replays them. `GuardedSourceSystem` has immutable public getters and a consuming `with_zero_domains`; it has no interior-mutability owner. Its internal field visibility is wider than lifecycle, but code outside lifecycle cannot remove or mutably borrow the program's Arc once that field becomes private.
* `solver/source.rs`: the guarded source constructor uses `SourceSystem::new`, with full native arity, no fixed coordinates and no inherited ordinary-family conditions. `same_source_context` compares every currently variable guarded context component: measure, roles, index map, coefficient variable map, zero domains, source IDs/domains/conditions and ordered polynomial rows. If guarded construction later admits prepared ordinary `SourceSystem` variants, the equality contract must additionally bind their active arity, fixed-coordinate pattern and inherited conditions; the present factory does not expose that variation.
* No `Clone`, `Default`, `Deserialize`, mutable program accessor, or other production `GuardedProgram` struct-construction path was found. The current `pub(super)` fields permit sibling-module access; lifecycle-private fields deliberately remove that opening. Lifecycle child tests still have privileged mutation access for testing the replayed APIs.
* RustFlow's `src/finite_density/guarded.rs` wrapper additionally binds physical arity and dummy symbols, validates terminal labels, and checks all output labels/conditions for escaped storage tails. Native verified reuse cannot replace those separate wrapper invariants.

## Required behavior of the new methods

1. Consume valid programs; never accept a raw rules vector as a verified input.
2. Compare source contexts by Arc identity or the complete exact equality check. A measure name, hash, rule count or common variable names alone is insufficient. In the distinct-Arc case, equality must include original source ordering, because proof source ordinals are positional.
3. Use checked addition for union counts and apply the cap before extending vectors. Duplicate rules count. Do not deduplicate, reorder, minimize conditions or silently substitute a different program.
4. Preserve `self` then `fallback` precedence. Retain each rule's exact domain, discovery domain, source trace, candidate, ordering and nonzero conditions. Application failures and numerical guard checks remain unchanged.
5. Recheck cross-program common coordinate order, including the empty/nonempty cases. Two individually valid programs can have incompatible coordinate permutations. Sector-specific fields need not be identical if the existing constructor permits them; use the existing common-coordinate contract.
6. Replace both former terminal sets with exactly the supplied current set. Check invalid occupation indices, missing required cuts and declared zero domains exactly as `new` does. Removing a terminal can expose an unresolved label and does not prove extra coverage.
7. Keep `new`, decode and explicitly replayed methods as real replay paths. Do not silently redirect an API named `*_replayed` to the trusted path, because its current privileged-corruption regressions deliberately exercise fresh validation.
8. Keep the proof object immutable for its lifetime. An `&Arc<SourceSystem>` getter is safe here: cloning the Arc does not grant unique mutable access while the program retains its own reference. No `Arc<Mutex<_>>`, mutable coefficient access, `into_parts` followed by unchecked reconstruction, or public forged-program deserializer should be introduced.

## Persistence and cache consequences

No proof-schema change is mathematically required if the encoded fields, replay semantics and source context are unchanged: a verified union should encode identically to the replayed union of the same ordered rules and terminals, and decode must replay both identically. Keep schema-v1 rejection and all existing v2 corruption tests. Native source/build fingerprints must change normally for the implementation revision; do not reuse an old executable/library hash for the prototype. A persistent cache still crosses the replay boundary even when its bytes were originally emitted from a verified live object.

The optimization only removes repeated in-memory proof reconstruction. It does not optimize fresh discovery/sealing, cold decoding, first program construction, exact context comparison for distinct allocations, storage validation, or rule application. Measure those costs separately before claiming a benefit.

## Focused test requirements

* Exact equivalence of verified versus replayed union and terminal replacement: rule order, duplicate count, current terminals, reduced coefficients, unresolved statuses and all nonzero conditions. Include complementary domains, a vanished index condition, a symbolic parameter condition, and a removed provisional terminal.
* Both Arc-pointer equality and separately built equal contexts succeed. Reuse the existing mismatch matrix: measure, roles, index positions, coefficient variables/order, source ID/domain/condition/row/row order and zero domains.
* Incompatible valid coordinate permutations fail; empty+nonempty works in either direction; empty+empty respects terminal validation.
* Duplicate-inclusive exact cap, zero-cap cases, and count-overflow helper behavior where practical. No rule is dropped to fit a cap.
* Every terminal rejection class: occupation below zero, required-cut power at or below zero, and a declared zero domain. Duplicate valid terminals may be set-deduplicated as in `new`.
* Encode/decode a verified union, replay on load, and compare the result. Preserve corrupt RHS, source guard, discovery/usable domain, source ordinal/condition, legacy-schema and foreign-context rejection tests.
* External compile-fail or privacy coverage for mutable access/struct construction, where the existing test framework supports it. Do not test that a deliberately corrupted lifecycle-private object is rejected by `*_verified`: the premise of that API is precisely that such an object cannot be constructed outside its trusted module.
* A test-only replay counter can show that new/decode/replayed APIs do replay, while verified operations do not. This is stronger evidence for the intended mechanism than a noisy timing assertion.

## Remaining limitation

This is approval of the design invariant, not of an unseen implementation. A source diff audit should check every new `Self { ... }`, helper visibility, persistence getter conversion and caller transition. Public claims should describe reuse of already replayed immutable native proofs, not omitted validation of external proof data.
