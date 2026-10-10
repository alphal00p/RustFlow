# Optional terminal-promotion memo: isolated integration draft

This patch is an unapplied, default-disabled draft. It reuses complete native unit reductions within one immutable program, then carries those results through a narrowly certified terminal addition. It does not replace descendant tails, discover identities, choose masters, specialize epsilon, or alter any reduction ordering.

The separate `2026-10-10-terminal-promotion-cache-experiment` passed eight semantic prototype tests. This revised encoded-storage integration has only been formatted and checked with `git apply --check`; it has not been compiled or executed. Its ten native tests and three RustFlow tests are proposed gates, not reported passes. The first formatting attempt failed on deliberately omitted unchanged child modules; the corrected invocation formats every copied source with `skip_children=true`.

## Exact reuse contract

`GuardedProgram` owns the memo privately. Entries are keyed by the unit input and both native reduction limits; the owner binds the exact source corpus, ordered rules, coefficient/index maps, common order and stopping set. Only completed reductions with no unresolved reason except `NoApplicableRule` can enter the memo. Native errors and all other failures remain explicit and are never cached. A hit restores the exact result, ordered conditions and logical application count. Entry decoding errors propagate rather than triggering silent recomputation.

Ordinary construction, decoding, verified/replayed union, verified/replayed terminal replacement, and resetting the memo discard every entry. Changing reduction limits selects a different key. There is no mutable insertion API or cross-round cache transfer.

The only exception is `with_promoted_terminals_replayed`: old terminals must remain, and each added label must previously apply to an empty, conditionless `NoApplicableRule` result. It then calls the existing full source replay. Only after replay succeeds are the corresponding cached residual records reclassified as terminals. Coefficients and ordered conditions remain untouched. Work-limited results are ineligible: a new terminal can bypass an application limit that blocked the old nonterminal call, so that case cannot be transformed.

This is complete-call reuse. In general, substituting already reduced descendant tails is unsound for the native applicability contract: live pending coefficients can cancel a guarded descendant before that descendant is applied. The draft does not perform that substitution.

## Retained-data limits

The defaults, used only when explicitly enabled, are 64 entries, 4,096 output records, 32,768 condition records, 262,144 polynomial terms and 64 MiB owned payload bytes. All limits are aggregate across the cache. A result too large for the quotas is returned normally and is not retained; it does not become a reduction failure. Invalid metadata allocations fail when enabling the memo. Encoder/decoder or native arithmetic errors propagate.

The precise byte charge is:

```
size_of::<Option<Mutex<UnitMemo<N>>>>()
+ slot_count * size_of::<Option<Entry<N>>>()
+ sum(entry.records.len() * size_of::<Record<N>>()
    + entry.conditions.len() * size_of::<CoefficientId>()
    + entry.state.len() + entry.atoms.len())
```

Every retained variable-length allocation is a boxed slice. The fixed entry layout includes the input label, reduction limits, application count and payload metadata. Only the native generated coefficient encoding is retained; no decoded coefficient object or decoded table remains in an entry. Polynomial-term counts conservatively charge every result numerator/denominator and condition, including repeated occurrences. Integer-bit and record-count checks precede encoding; exact encoded-byte admission follows encoding.

The charge excludes allocator bookkeeping, the immutable program, returned reductions and transient encoding/decoding/CAS scratch. It is an enforceable retained payload metric, not a total process-memory or arithmetic-work bound. Codec limits are intentionally unrestricted only for internally generated, privately owned records; they do not authorize accepting untrusted external encodings.

## RustFlow integration

`WeightedClosureOptions::unit_reduction_memo` defaults to `None` and requires active-target closure when set. Each discovery round begins with a fresh memo. Existing cancellation checks before unit calls remain. Incomplete rounds clear the memo before retaining their program. The final terminal addition still performs full independent source replay, and both existing loops remain: every original weighted target and every final basis derivative is reconstructed. Successful results clear the cache before the program is returned or persisted. Cache contents are not checkpointed.

The original `native_rule_applications` count remains the logical count. `memoized_rule_applications` records the portion supplied by exact hits; their difference is the performed native reduction-application count. These counters do not count rule replay, encoding, decoding or the one-step terminal-admission checks. Hit/miss counters apply only while a memo is enabled. Peak quota counters describe retained data, not scratch.

The runtime diagnostic harness accepts `RUSTFLOW_WEIGHTED_UNIT_MEMO=true`, with optional `RUSTFLOW_WEIGHTED_UNIT_MEMO_ENTRIES`, `_TERMS`, `_CONDITIONS`, `_POLYNOMIAL_TERMS` and `_BYTES` caps. Its full and occupied closure reports expose the counters and configured caps. Existing environment/default behavior leaves the option disabled. No public CLI/Python option is added by this draft.

## Required gates before promotion

- Compile native library tests and RustFlow all-target/Python tests; run the existing full guarded and finite-density gates plus all thirteen new tests.
- Native tests compare ordered terms, residuals, guards and logical counts, including 90 work/pending/capacity combinations, canceled unused guards, retained used guards, work-limit exclusion, exact cap decline, invalidation, source mismatch, corrupt owned payload errors, and mandatory replay rejection of a corrupt proof.
- RustFlow tests compare the compact fixture's complete basis, matrix, weighted target rows, conditions, persisted program bytes and logical counts with/without memo; exercise disabled-active rejection and cancellation before final replay.
- Run matched lower-loop numerical controls and, only after those gates, the actual expensive target/closure with the option explicitly enabled. No speed or real-target memory claim has been established by this draft.

`base-source-map.json` binds the copied originals; `source-map.json` binds the draft and patch. No production files or shared build artifacts were changed to create this report.
