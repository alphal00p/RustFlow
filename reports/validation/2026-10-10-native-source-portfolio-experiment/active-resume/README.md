# Bounded timeout continuation

This standalone wrapper resumes only a completed exit-124 timeout from the same native engine and complete source corpus. It preserves the original targets, search settings, source guards, zero domains, measure identity, global 16-round limit, 4,096-request history, 1,024-label frontier cap and 65,536-rule cap. The continuation receives an additional bounded 600 seconds; prior elapsed time, including unfinished work, remains in the cumulative resource report.

The last entry in `progress.json` is the commit marker. Its matching round JSON and native binary must exist and agree. Native `decode_generated` replays the proof against the exactly reconstructed full source system. The wrapper restores the frontier, attempted exact points, next requested labels, encountered conditions and cumulative history. It rechecks saved target maps and processed derivative maps before starting the next round. It imports no rules from the original corpus file; only this engine's newly discovered parent checkpoint rules are restored.

Incomplete trailing round files are ignored and recorded. A corrupted or partially written progress marker is rejected; there is no speculative recovery. A closed result, terminal-bearing program, frontier/request/round exhaustion, changed configuration or different native engine is rejected. Original reports remain unchanged. A copy of the previous commit marker and complete checkpoint pair makes a timeout before the first continuation round resumable without fabricating progress.

The baseline and proposal use one source with cfg-gated portfolio statistics. The existing native discovery, reduction, exact derivative sums and final one-program closure audit are unchanged. `checkpoint-controls.json` records 25 metadata checks; these are separate from the native replay performed when a continuation starts. An initial compile failure due to an unspecified `Vec<serde_json::Value>` type is preserved in `baseline-build-attempt-1.*` and `failed-inference-source.rs`.

Run only after the parent has completed and continuation is authorized:

```text
python active-resume/run.py baseline active-pilot/polynomial62-baseline polynomial62-baseline
python active-resume/run.py proposal active-pilot/polynomial62-proposal polynomial62-proposal
```

Paths in these examples are relative to the source-portfolio report directory. This is a validation-only closure diagnostic. It does not compute a physical amplitude or establish three-loop or four-loop acceptance.
