# Mixed-energy polynomial source controls

These validation-only controls retain the production reducer, the original polynomial completion domain, every original source row, and every original zero domain. They import source definitions from the failed production cut-3 checkpoint and discard its rules. No period values are read, no inverse-energy domain is admitted, and no production policy is changed.

The public input/geometry owner reconstructs the three-loop cut `[3]` and proves the occupied energy is exactly `E0 = F6 − F8`. Neither contributing factor alone is proportional to `E0`. The existing policy consequently emits 54 rows: the 52 legacy rows and two guarded raw Ward rows, without the shifted temporal/boost groups used for cut 0.

For each contributing ordinary completion, the experiment reindexes the complete temporal `U0` or common-boost identity by `old index = new index + 1`. It first intersects the original row with the polynomial domain, then substitutes the index in every coefficient and condition and pulls back the guard to `new index <= −1`. This is a presentation of original polynomial identities, not division by a mixed energy. The probe checks that potentially positive completion images vanish at the old exponent-zero boundary. Common-boost actions on every physical scalar factor vanish exactly. Added rows are appended to the complete original corpus.

Seven matched arms cover the original corpus, shifts in slot 6 alone, slot 8 alone, both slot orders, temporal-only shifts, and boost-only shifts.

The initial 15 actual requests all give 15 applicable native rules and 174 total RHS terms in every arm at depth 3/domain cap 1. RHS coefficients and condition multisets are identical; some condition sequences reorder two rows. This stage shows no benefit.

The later test selects 32 actual labels from completed cut-3 round-1/2 submitted requests and native uncovered frontiers. Selection is lexicographic: eight labels with negative slot 6 only, eight with negative slot 8 only, eight with both negative, then eight remaining frontier labels. Every arm uses a native requested ray at depth 3/cap 1, followed by an exact point at cap 1 only for `NoApplicableRule` or `ConditionVanished`. Every label needs the point fallback, and every arm applies at all 32 labels.

Slot-6 temporal shifts change two later RHS rows: 7 to 10 terms and 24 to 12 terms, reducing the sum from 141 to 132. Both-slot and temporal-only arms give the same two changes. Slot-8 and boost-only arms change no RHS. Changed condition lists are retained in full. The result is a modest presentation difference, including one local regression; it establishes no new applicability coverage, closed basis, or numerical improvement. No full pilot or production promotion follows from this experiment.

Both stages use native source replay and generated-program encode/decode checks. `build-binding.json`, per-arm resource/binding files, and `preserved-inputs.json` bind the source, executable, dependencies, inputs, and failed K checkpoint copies. `later/selection.json` binds the downstream selection. `summarize.py` reproduces the two summaries. The initial standalone compile failure omitted `CARGO_CRATE_NAME`; its log and binding are preserved as attempt 1, and the successful build records the required compile-time value. Timings were measured on a shared host and are not performance comparisons.

`archives.json` maps every archived proof to its original path, raw hash and deterministic gzip hash; all archives were decompressed and compared before raw removal. Small source corpora and checkpoint copies remain raw. `manifest.json` freezes the final experiment files. Historical JSON references remain unchanged and can be restored using the archive map.
