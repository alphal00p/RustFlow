# Bounded lookahead after a native direct match

This isolated search experiment has no measured benefit on the 19 finite-density
reduction probes. All 19 resulting right-hand sides and condition lists are
identical to the baseline, with 192 terms in total. Production search remains
unchanged, and no full closure pilot was run with this variant.

The proposed search retains the first direct candidate, withholds matching
direct rows from modular elimination, and permits up to 512 further source-row
attempts. It uses the existing native GPLU and exact materializer, limits an
indirect trace to 64 rows and 16,384 terms, and selects a later indirect result
only when its exact right-hand side is strictly shorter. Guard-rejected attempts
count against the additional budget. Ordinary search does not enable this
behavior; guarded proofs use a distinct experimental persistence tag.

Only one of the 19 probes actually enters this lookahead. Its one-term direct
candidate is retained after the alternative gives ten terms. The motivating
odd-integral example reveals a different limitation: its first full-corpus
winner already comes from a multi-row proof, so a policy triggered only by
direct matches does not run. A fresh virtual-source-only control still proves
that example zero, while the full corpus produces two terms. The initial
first-direct-match explanation was therefore incorrect for that example.

All 11 focused search tests pass, including exact replay, condition retention,
budget fallback, ordinary-path behavior and propagation of fatal errors.
Baseline and experimental proofs reject each other's persistence schemas.
Fallback is guaranteed only for the recorded budget/exhaustion cases; optional
arithmetic or sample failures still propagate, and selecting a candidate before
sealing does not guarantee the same conditional coverage as the baseline.

The library and test executable were compiled in a separate temporary source
copy against the hashed production dependency graph. Library compilation took
183.890 seconds and test compilation 278.351 seconds. The resource wrapper's
generic "compilation excluded" label does not apply to those two invocations;
their build bindings record the correction. The failed initial test compilation
and the missing `CARGO_CRATE_NAME` environment fix are preserved. Production
source hashes were verified unchanged. Numerical timings from this shared host
are not a controlled performance benchmark.

`negative-lookahead-summary.json` records the outcome.
`negative-lookahead-artifacts.json` freezes the original 107 evidence files;
`manifest.json` additionally binds this reader summary and that frozen manifest.
Neither the search controls nor the source-subset zero witness establish a
closed differential system or a three-loop numerical result.
