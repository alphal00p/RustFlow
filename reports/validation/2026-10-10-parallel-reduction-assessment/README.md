# Bounded native unit reduction batches

This report develops an isolated native batching API for independent complete unit reductions. Production source and the shared target directory are unchanged. No timing or closure improvement is claimed.

- `send-sync-build.json` records the successful compile-time Send/Sync gate for actual immutable native and RustFlow owners.
- `j3-completed-checkpoint-observations.json` binds the completed J3 discovery-stage observations. Those one-step application statuses are not a final reduced basis, and checkpoint silence is not a measured attribution of reduction cost.
- `design.md` states ordered cache semantics, explicit non-fail-fast outcomes, actual versus logical work, resource limits, and the narrow proposed frontend callsite.
- `independent-batch-review.json` is a static review of the original draft; it records no execution. Its source hashes identify the reviewed revision.
- `implementation-draft/revision-01/` preserves the first compile input and patch. It adds nine meaningful batch regressions, including both review suggestions.
- `build-isolated.py` compiles the complete actual native crate in `/tmp` against six hash-verified copied dependencies. It never attaches private new fields to an old native type and never writes Cargo's shared target.

The build/test records are the execution authority. A separate frontend draft is owned by the boundary agent; it is not applied or validated by this native report. New batched calls are opt-in, and existing serial APIs are unchanged. Final source replay, original weighted target reconstruction and all final basis derivatives remain mandatory at the RustFlow boundary.

Final native validation: `tests-build-02.json` records the complete isolated crate build (281.24 s), and `guarded-tests-02.json` records **111 passed, 1 pre-existing ignored, 0 failed** in the full guarded filter (0.13 s harness). This includes ten new batch regressions; the earlier nine-test run is historical, not additional unique coverage. The final five-file patch is `implementation-draft/native-batch.patch`, SHA256 `0f19e0018e8a8a39644839f0885daf2723935e4fad096af5c1f755221d57f01d`, and passes `git apply --check`. It has not been applied to production.

The controlled tests establish exact native outputs, ordered conditions, failures, logical application counts, cache-entry/LRU equality, per-item usage snapshots, and enforced caps for their fixed process state. They do **not** establish universal identical cache bytes, cache hits or usage peaks against frontend executions with different process-global Symbolica registrations. `frontend-read-only-review.json` records the separate source review and this qualification.

The same frozen native revision also built successfully as an isolated rlib in 183.22 s (`library-build-02.json`). The user requested a safe pause as it completed. `safe-pause-state.json` binds the library, patch and completed gate; no agent-owned process remains active and frontend compilation has not started.
