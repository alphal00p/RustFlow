# Native source-portfolio experiment and closed polynomial-source controls

**Both the unchanged native reducer and the isolated portfolio prototype close the same occupied singleton channel of the three-loop massless chain on a seven-element spanning basis.** Their final exact reductions of both original weighted outputs and all seven basis derivatives are identical. The native portfolio is not needed for this closure result and remains an isolated experiment. No basis-minimality, complete finite-density amplitude or numerical-period result is asserted here.

Both arms use the same independently justified 62-source polynomial corpus: the original 52 rows, four shifted temporal-U rows, four shifted common-boost rows and two narrow raw C1 Ward rows. No reciprocal energy powers, custom elimination backend, manually chosen masters or reference values are added. The source identities and original two output sums (three input labels including the raised fixed numerator) are bound in the saved contexts. The frozen corpus is `strongest62-boost-ward/program.bin` in the sibling `2026-10-10-polynomial-raw-ward-attribution` report, SHA256 `f482d60e6aee6f3187be2eac039b26fec8fb8cc7cef185d577b5d2f8ddad14a2`.

## Complete active-closure comparison

| Measure | Unchanged native | Isolated portfolio |
| --- | ---: | ---: |
| Final basis size | 7 | 7 |
| Completed total rounds | 16 | 15 |
| Retained native rules | 3,553 | 3,324 |
| Largest completed frontier | 520 | 517 |
| Original weighted output rows audited | 2 | 2 |
| Actual basis derivative rows audited | 7 | 7 |
| Retained conditions, as stored | 677 | 598 |
| Nonzero-rational-associate condition classes | 152 | 129 |
| Initial bounded process | timeout, 300.010 s | timeout, 300.008 s |
| Authenticated continuation | pass, 364.379 s | pass, 368.694 s |
| Cumulative process wall time | 664.388 s | 668.702 s |
| Continuation peak child RSS | 50,316 KiB | 57,904 KiB |

These processes overlapped on a shared host. Wall time includes the interrupted parent's work and any repeated unfinished-round work, and excludes compilation and Nix startup. Fewer rounds or rules in the proposal do not establish a performance advantage. Build resources are recorded separately.

The common bounds are 16 total rounds, 1,024 frontier labels, 4,096 requested labels and 65,536 rules. Every queried ray and exact-point fallback uses native depth 3 and a domain budget of 1; the optional direct-zero attempt cap is 8,192 per round. The discovered direct-zero prepass contributes zero rules in these runs. Actual weighted sums, rather than individual term labels alone, determine the active frontier. Newly exposed frontier and derivative labels must receive their own native discovery requests. Retired reducible labels are recorded, not promoted to masters.

At closure, the runner rebinds exactly the discovered basis as native terminals, replays that program in its complete source context, and reduces both original weighted target sums and the actual derivative sum of every basis element. It asserts no failures and no unresolved leaves before writing `closed.bin`. [active-closure-comparison.json](active-closure-comparison.json) retains the exact common basis, all nine final rows, histories, result/program hashes and runtime scope. [compare-active-closure.py](compare-active-closure.py) verifies exact saved-row equality; it does not run another reduction or evaluate a period.

The baseline native library is SHA256 `ac85472324303a84e909023a9aac6fe77f4993934985519fd88275cae0e1aad0`. The isolated portfolio library is SHA256 `4d9423cd6fca369f184b29013998f13180556d3f6cc6767d445ed8c8ebcf6e45`; its distinct persistence tag rejects cross-loading baseline proofs. All prototype source changes, exact dependencies, generated arities, build commands and binary hashes are retained in `native-build/` and the probe build bindings. This experiment changes no production native source.

## Retained conditions

The exact condition sets differ. After only nonzero rational scaling and repetition are identified, 111 classes are common, 41 occur only in the baseline and 18 only in the proposal. No factorization, product splitting or common-root argument is used to identify further classes. Neither set is declared equivalent to, or a subset of, the other. Every original condition remains in the result and [active-closure-condition-comparison.json](active-closure-condition-comparison.json).

As a separate admissibility check, substituting exactly `epsilon = -5/4` makes every retained condition a nonzero rational multiple of `eta^k`, with nonnegative integer `k`. Thus both recorded guard sets are nonzero for every positive eta at that epsilon. This evaluates no integral and does not license eta zero or any endpoint limit.

## Earlier point controls and prototype scope

The portfolio first retains the ordinary full-corpus native candidate, then tries two bounded original-row subsets. Trial source rows retain their original ordinals and the full original source context; native source-domain validation, uniform zero projection, descent order, exact materialization, sealing and replay remain mandatory. An alternative is selected only after native certification, a strictly shorter exact RHS and the conservative retained-condition gate. The full candidate survives failed or tied trials. This is a discovery heuristic, not another reduction backend or a new physical source theorem.

Each optional arm is capped at depth 3, 2,048 attempted source rows, 512 accepted rows, 64 exact trace rows and 16,384 exact trace terms. The 19-point original-52 comparison gives 192 versus 181 RHS terms, with two strict improvements. The 19-point strongest-62 comparison gives 143 versus 127, also with two strict improvements. All 19 condition gates pass in each comparison; every point proof is independently sealed and encoded/decoded/replayed. Nineteen focused implementation controls and both cross-schema rejection controls also pass. These local results alone did not establish closure; the later resumed controls above do.

[design.md](design.md) is the preserved historical design-stage document, whose initial status statements predate the implementation and completed runs. Point summaries and hardened comparisons retain their original narrow scope. The first condition parser's ambiguity was preserved and corrected by the separately bound fail-closed `condition-polynomials-v2.py`; the hardened comparisons are the final point-condition evidence.

## Timeout, resume and archival provenance

Both first active pilots ended with explicit exit 124 before closure. They were resumed only from their last committed marker and native program: baseline round 007 and proposal round 006. The continuation retains the global budgets and cumulative requested-point history, checks the same source/target/options/library identities, decodes with full native source replay, and replays saved original targets and completed derivatives before resuming. Twenty-five metadata admission controls and an independent static review are in `active-resume/`. A failed initial continuation-wrapper compile and its source snapshot remain preserved; the successful build has a separate binding.

The raw parent checkpoint pairs and their configuration, progress, resource and run bindings are retained, as are both final `closed.bin` files and final result JSON files. All other archived payloads have deterministic gzip files, raw and compressed SHA256 values, sizes and byte-for-byte restoration checks in [archive-map.json](archive-map.json). Frozen JSON references deliberately retain their original paths; restore an archived payload to that path before rerunning a historical script. [manifest.json](manifest.json) binds the finalized report, excluding generated Python bytecode caches and itself. No input corpus outside this report is moved or rewritten.
