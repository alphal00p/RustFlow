# Shifted polynomial temporal-U source pilot

The fresh native 56-row pilot remains unclosed. It reaches the 1024-label frontier
cap after four completed rounds, with frontiers 44→152→528→1025 and native rule
counts 15→58→235→853. The last original output rows contain 19 and 1061 terms;
734 cumulative points have been explicitly requested and eight actual frontier
derivatives have been processed. No positive-energy completion appears in those
remainders. The process takes 163.911seconds and 79,932 KiB peak child RSS.

The matched historical 52-row production-order baseline has frontiers
54→203→633→1025, rules 15→68→298→1015, and final original rows 19/1506 terms
at the same cap. Its prebuilt runtime is 250.672 seconds. These unoptimized wrapper
runs are not a controlled machine-load benchmark. The exact count improvement
is real; it does not establish a finite closed basis or numerical flow.

Only the source presentation differs from that baseline. Four original temporal-
U source rows are appended after substituting a6→a6+1 in every coefficient and
condition, shifting every integral label+1 in slot6, and pulling the original
domain back exactly to a6<=-1. The resulting rows are twice the shell-normal
identity on their original polynomial domain. All 52 original rows, original
output coefficients, source conditions, zero domains and storage-tail bounds
remain. No inverse-energy admission, new zero theorem, target value, chosen
master or historical rule is introduced. A new measure presentation epoch binds
the added rows. The native source/replay owner is the unchanged v2 production
library, SHA ac85472324303a84e909023a9aac6fe77f4993934985519fd88275cae0e1aad0.

The active-closure scheduler, derivative sums, strict failure handling and
explicit queried-point history are identical to the baseline. Limits are 16
rounds, 1024 frontier, 4096 cumulative requests, 65536 rules, native depth 3,
one ray-domain plus one exact-point fallback per query, 8192 direct-zero attempts
per round and a 300-second timeout. The process exited 0 because the diagnostic
finished and recorded its frontier-budget failure. Successful closure would
still require the same final original-target and complete-basis-derivative
native audit. That audit was not reached.

`source-transform.patch` records the complete runner change;
`source-transform-binding.json` ties it to the unchanged historical baseline
body and frozen small-control manifest. `independent-source-audit.json` records
a separate read-only audit of substitutions, domains and geometry. The compiled
source/library/binary hashes, actual command and resources are retained.
`summary.json` records completed rounds and the exact matched baseline.

Verbose round programs/metadata are losslessly gzip-archived through
`archives.json`. Original hashes remain in the summaries, and every archive was
restored byte-for-byte before raw removal. Use `gzip -dk FILE.gz` to restore a
referenced raw path before rerunning `summarize.py`. `manifest.json` hashes the
final stored files. No active probe remains and no production source changed.
