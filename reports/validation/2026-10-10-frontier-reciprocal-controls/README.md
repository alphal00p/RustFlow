# Later-frontier polynomial and reciprocal source controls

The sampled later frontier shows **no additional target coverage or smaller RHS from the reciprocal corpus**. Both corpora produce a native rule for 18 of 19 labels, return 38 RHS terms in total, and leave the same label explicitly unresolved. Seventeen RHS maps are identical. The two other maps have the same size but use different descendants or coefficients; no equivalence between those maps is asserted here.

| Corpus | Sources | Applied | Unresolved | RHS terms | Proof-source rows | Process wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Strongest polynomial | 62 | 18 | 1 | 38 | 112 | 0.8061 s |
| Widened normalized boost | 56 | 18 | 1 | 38 | 106 | 0.8066 s |

Both processes used 15,464 KiB peak child RSS. These are complete prebuilt-probe runtimes, excluding compilation. Each had an explicit 300-second timeout.

## Selection and native proof contract

The input is the completed `round-007.json` from `2026-10-10-native-source-portfolio-experiment/active-pilot/polynomial62-baseline`, whose frontier contains 290 distinct labels. Selection uses only labels. It stratifies by the exact required-cut power in slot 0 and occupied-energy power in slot 6. All nine strata are represented: two endpoints per stratum ordered by upper-occupation power, virtual-energy powers and full label, plus the lower median of the largest stratum. This covers cut powers 1–3, occupied-energy indices −3–0, and upper-occupation powers 0–2. [selection.json](selection.json) records the complete rule, populations, selected labels and frozen input hash.

The unchanged baseline generic probe reconstructs source rows, source guards, roles and certified zero domains from each explicit v2 corpus. It discards every historical reduction rule. It searches one exact point domain per label, with depth 3, seed 0 and a domain budget of 1. Every discovered program is encoded, decoded against its complete source context and replayed. Status, exact RHS coefficients and all conditions match before and after persistence. The source corpora intentionally differ in their added identities and energy-index admission; this is a comparison of the requested full corpora, not an isolated estimate of widening alone.

The polynomial corpus remains the raw, unchanged `strongest62-boost-ward/program.bin` in `2026-10-10-polynomial-raw-ward-attribution`, SHA256 `f482d60e6aee6f3187be2eac039b26fec8fb8cc7cef185d577b5d2f8ddad14a2`. The reciprocal corpus is `normalized-first/program.bin` in `2026-10-10-normalized-boost-reciprocal-experiment`; its exact hash is in [reciprocal56-binding.json](reciprocal56-binding.json). The baseline native library is SHA256 `ac85472324303a84e909023a9aac6fe77f4993934985519fd88275cae0e1aad0`. Executable, source, dependency, command and point hashes are retained in the copied build binding and per-process bindings.

## Energy images and conditions

Neither arm uses a positive occupied-energy source seed or returns a positive-energy RHS term. The reciprocal proofs contain six nonzero-coefficient positive-energy source images, all for sample point 7. Every such image has cut power zero and vanishes by the existing required-cut role. No positive-energy image survives the known zero projection. The polynomial proofs have none. [proof-energy-images.json](proof-energy-images.json) records the exact specialized source coefficients and labels; this audits the proof-source images and final one-step RHS, not every transient internal elimination row.

Raw condition lists differ at two points. Exact polynomial comparison shows that all 19 pairs define the same conditions after removing nonzero rational scaling and repetition. The original condition lists are preserved without alteration in [summary.json](summary.json) and the full results. The common unresolved label is `[1,1,0,0,1,0,0,0,0,0,0,0]`, with `NoApplicableRule` and an explicit `SearchExhausted` gap.

These controls do not prove closure, numerical values, or that reciprocal admission is unnecessary for every later label. They provide no new evidence for promoting that admission from this sample. No production code, existing frozen report, source corpus or active pilot input was modified.

## Reproduction and archives

`select-and-run.py` creates the deterministic sample and executes the two bounded controls. `analyze.py` specializes proof-source coefficients using exact rational polynomial arithmetic and compares conditions; it performs no integral elimination. Before rerunning analysis, restore any archived JSON/proof payloads using [archive-map.json](archive-map.json). Archives are deterministic gzip with original and compressed SHA256 values and byte-for-byte restoration checks. `manifest.json` binds the final report files.
