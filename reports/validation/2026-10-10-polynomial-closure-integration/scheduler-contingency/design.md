# Contingency: bounded requested-ray and exact-point discovery

Status: design only. No production or native change, compiled experiment, or recommendation to enable this schedule. First run the frozen production PolynomialClosure policy. This contingency is justified only if the source-equivalent production attempt fails to transfer the already proved singleton closure.

## Evidence and the smallest proposed interface

The unchanged-native pilot in `../../2026-10-10-native-source-portfolio-experiment/active-resume/pilot.rs` closed the cut-[0] three-loop family with seven spanning labels after 16 rounds. Its initial requests were the three original labels and twelve raw derivative labels. Each uncovered requested label received a one-domain symbolic-ray search followed, when necessary, by a one-domain exact-point search. The saved baseline progress records 3,153 exact-point fallbacks, including all fifteen first-round requests. These counts establish a scheduling difference, not that every fallback was necessary.

Production `reduction/active.rs` starts with roots only. Its residual schedule fully reduces requests and searches descendant leaves. Its exact-point search currently handles `ConditionVanished`; ordinary `NoApplicableRule` leaves receive rays. Existing options therefore do not reproduce the successful schedule.

Add one optional schedule value to `WeightedClosureOptions`, preserving the existing branch when absent:

```rust
pub struct RequestedRayPointOptions {
    pub max_domains_per_ray: usize,
    pub max_domains_per_point: usize,
}
// Proposed field; None by default.
pub requested_ray_point: Option<RequestedRayPointOptions>,
```

This option selects a cohesive schedule, including initial raw derivatives, original-request applicability probes, and exact fallback. It does not add another reduction engine, source presentation, integral order, source theorem, or basis selector. Reuse the current `GuardedContext::discover`, native one-step `apply`, `union_replayed`, `with_terminals_replayed`, actual-sum reduction and final audit.

For the first implementation, require active-target closure, frontier searches, positive reuse/ray/point caps, requested-index rays, and positive total discovery budget. Reject simultaneous residual allocation (`max_domains_per_residual > 0`), interval refinement (`guard_refinement.max_passes > 0`) or priority hints instead of silently ignoring those controls. This gives one independently reviewable scheduling branch. Existing validation must distinguish the new branch from its current requirement that active closure implies residual discovery. No existing default or caller changes automatically.

## One bounded round

1. Seed round zero with nonzero original target labels and every raw derivative label returned by the existing `FixedShellDeformation`. Validate admission and storage tails before any search. Apply the existing cumulative historical-request cap to this enlarged set.
2. Rebind retained rules to an empty stopping set. Optionally run the existing bounded direct-zero prepass, keeping its successful proofs first, retained rules next. Preserve its completed-point memoization; this differs from the pilot's redundant zero searches but changes no proof or ordinary rule order. The successful pilot added no direct-zero rules.
3. Freeze that round's retained program for scheduling. Visit the actual requested set in deterministic lexicographic order. Do not replace requested labels by the descendants of full reductions, and do not let earlier fresh programs suppress later requests within this batch.
4. Probe each requested label with native one-step application. `Applied` or `Zero` means no new search is needed; record conditions from the application and validate output labels. A `Terminal` with this empty stopping set is an invariant failure. Only `NoApplicableRule` and `ConditionVanished` permit discovery. Propagate invalid occupation, unsupported powers, coefficient poles, non-descent, work limits and native errors explicitly. An applicability probe is not a completed reduction or a basis certificate.
5. For an eligible label, form the existing anchored orthant: zero coordinates fixed zero, positive coordinates bounded below by the requested value, negative coordinates bounded above by that value. Intersect with the admitted domain, including polynomial-completion bounds and fixed-zero storage tails. Run native discovery with the configured ray cap, depth and seed, and replay its program against the unchanged source context.
6. Probe the ray program at the original label. If it returns `NoApplicableRule` or `ConditionVanished`, run a second native discovery on the exact singleton box with the point cap. Preserve both ray and point programs and every returned native gap. A still-unresolved exact point remains explicit. Other failures stop this branch; they are never turned into fallback requests or cancelled.
7. Combine fresh programs in request order, ray before point, using a balanced tree of existing native replayed unions. Then union them after the frozen retained program. Zero-first/retained-first precedence matches the pilot. Count duplicate rules toward the existing rule cap. Do not export a provisional terminal set during discovery.
8. Rebuild reachability from the actual original weighted sums using the existing active closure traversal. Every newly exposed prospective basis label and each raw derivative label must have been explicitly considered in a completed scheduling transaction before it is differentiated. Record labels deferred by allocation exhaustion separately; do not add them to completed history. Completion of a failed bounded search permits later provisional exploration, but never certifies coverage.
9. Closure still requires rebinding one final native program to the proposed spanning set and reducing every original target sum and every complete basis derivative sum to that set without any residual or failure. Preserve the existing exact coefficient coalescing within each sum, admission checks, condition collection and historical-retirement reporting. Never promote an unresolved derivative to a master.

For `ConditionVanished`, the failed condition is a recorded applicability failure, not an asserted nonzero fact. The old rule remains in the retained program; a later valid rule may cover the exact label through existing native selection. Reapply after union and retain the failure if no valid replacement exists. Preserve every actual nonzero condition returned by successful native applications and all original target/source conditions. A later cancellation does not erase an acquired condition.

## Allocation, cancellation and determinism

`discovery.max_domains` is the total conservative allocation cap for this round, including rays and fallbacks; `max_rounds` bounds the full preparation. Existing application, historical-request, frontier, rule-count and direct-zero-attempt limits remain independent. Report allocated caps, not an invented measured native visit count.

Before starting a new point transaction, reserve `ray_cap + point_cap` with checked arithmetic. Charge the ray cap when called; release the point reservation if the ray applies. If the reservation does not fit, leave that label explicitly deferred and do not mark it attempted. A deterministic prefix with explicit remaining obligations is preferable to running a ray and silently losing its required point fallback. The pilot profile (history cap4096, ray1, point1, round-domain cap8192) can fund every round's worst-case transactions without reducing its historical bounds. Smaller profiles may return a bounded incomplete outcome; allocation gaps are scheduling records, not fabricated native proof gaps.

Check cancellation before every applicability probe, native discovery call, replayed union and active-sum reduction, and before checkpoint publication. A cancellation during a native call is observed after that existing bounded call returns; do not claim mid-call interruption unless the native API supplies it. Never publish an interrupted round as complete. Balanced unions avoid quadratic repeated replay; they preserve exact left-to-right rule precedence. All calls use the actual native source system and its existing common order.

The new helper should return the replayed discovery result plus the subset of completed requested-point transactions and explicit deferred transactions. Active closure must update history with the completed subset, not unconditionally extend it with every requested label. A label successfully applicable in the frozen retained program counts as explicitly considered without consuming domain budget.

## Source identity, checkpoints and conditions

The physical source epoch must stay unchanged: this schedule adds no source identity or origin theorem. Native programs remain bound to the exact source corpus, measure/support/deformation identity, roles, index/coefficient variables and ordering, zero domains, physical arity and capacity. Use the existing compiled source digests and native replay checks; do not pretend a schedule version is a new physical proof.

Add a separate `requested-ray-point-v1` schedule tag to checkpoint metadata and record both caps, total budget, depth/seed, application limits, rule limit, raw-derivative initialization, empty stopping set, rule precedence, options and source bindings. Each transaction records label, initial status, ray and point allocations, statuses, returned gap IDs, successful conditions, and completion/defer reason. Keep historical native gaps and final actual-sum conditions even when a finite connection closes.

Publish program bytes and round metadata atomically with hashes or a final commit marker binding both; metadata also records the exact completed-request history, next-needed set, frontier and cumulative counters. Existing root checkpoints are evidence rather than a complete resume API. Do not add resume implicitly. A future reader must verify the same source and schedule bindings, decode/replay natively, and recheck saved actual target/derivative rows before restoring state.

## Meaningful required regressions before use

- A small native source system whose bounded ray misses a concrete request but the exact box produces a replayed recurrence. Assert both native gaps survive, the fallback covers the original label, and no hand-written rule or chosen basis was installed. Derive/minimize such a corpus from the saved pilot instead of mocking discovery.
- The existing coupled-guard fixture at `[2,2]`: the retained rule's index condition vanishes, exact discovery produces an independently replayed replacement, and the physical parameter guard remains. With no replacement or inadequate allocation, retain the exact failure and report incomplete.
- A physical compact-ball source family with two weighted outputs and native auxiliary differentiation. Verify final target reconstruction and derivative closure after encode/decode, including actual-sum cancellations. Record that initial raw derivatives were considered; do not merely assert option values.
- A budget boundary with two eligible requests: enough allocation for the first ray/point pair but not the next. The second stays deferred and cannot enter completed history or a closed basis. Separate tests cover no remaining allocation, checked counter overflow and cancellation between transactions without a completed checkpoint.
- Failure preservation: `CoefficientPole`, invalid occupation, non-descent/unsupported-power where constructible, and reduction work-limit failures must not trigger successful fallback or cancel across terms. Padded-tail and positive polynomial-completion admission must fail before native discovery.
- Precedence and context binding: an existing replayed rule remains selected ahead of a fresh alternative; current terminals are empty; source/arity mismatch and rule-budget overflow fail. All successful guards survive union and serialization. A changed schedule/configuration cannot be resumed from a mismatched round marker.
- Default-off regression plus the existing native guarded suite and ordinary/cut numerical gates. Finally, run the actual public source-equivalent three-loop singleton with the saved pilot bounds, then the other occupied sectors and full amplitude. A singleton closure alone is not complete three-loop numerical validation.

No implementation is included: production transfer has not yet been tested after the current fixture corrections, and widening the scheduler surface before that result would be premature.
