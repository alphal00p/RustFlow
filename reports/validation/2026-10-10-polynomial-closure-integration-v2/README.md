# Optional polynomial source policy integration

The coherent release build and all seven regression groups pass: 105 tests,
including the exact 62-source public-factory comparison. The 84 earlier native
guarded tests are reused only after verifying identical native sources, features,
test executable and library. The complete massless two-loop fixed-dimension and
Laurent regressions and the massive two-loop regression also pass.
All four three-loop attempts fail during preparation; no three-loop prediction
or four-loop numerical acceptance is claimed.

The policy is `WeightedSourcePolicy::PolynomialClosure`, serialized as
`polynomial-closure-v1`. The `v2` in this report names the second integration
attempt, not a different physical source policy. Production keeps the legacy
default and the unchanged RustRed search. The policy adds exact polynomial
temporal/boost source shifts and, only under the bound singleton theorem, the
raw surface Ward pair. It widens no massless energy-power domain.

The [initial attempt](../2026-10-10-polynomial-closure-integration/polynomial-closure-attempt-summary.json)
compiled successfully and passed the exact source comparison, but two new
nonunit-routing test fixtures assigned inadmissible charge-2 edges. The
correction uses rational half-scaled loop routing, integer loop charge 2 and
unchanged unit edge charges. It asserts the actual determinant and nonunit
energy map. Only those two test bodies changed between source snapshots; no
production admission rule or generator changed.

The release capture binds 2,480 source assets and twelve build artifacts. The
source-equivalence test compares every ordered source row, label, domain,
condition, role and certified zero domain with the independently derived
experimental corpus. It discards all historical rules, permits only nonzero
rational scaling of added rows, and keeps both actual condition lists. Its
new source context roundtrips through native persistence.

The build, gate and prediction launchers in this directory preserve separate
resource records and reject overwritten attempts. Prediction launchers load
only physical definitions and actual production factories. Comparisons require
complete saved predictions and matching source, build, input and proof
identities before reading the frozen independent references. The full
three-loop fixed-dimension timeout is 1,800 seconds because the successful
standalone single-cut closure alone took about 664 seconds; the algebraic
budgets and numerical profiles retain their prior values. These shared-host
runs are not a timing comparison.

The [checkpoint summary](polynomial-closure-v2-checkpoint-summary.json) binds
the completed gates, source snapshot, build, numerical resources and failures.
The [lower-loop report](polynomial-closure-v2-lower-loop-regressions.json) records
the following saved-prediction comparisons:

| Numerical regression | Independent references | Profile refinements | Harness time |
| --- | ---: | ---: | ---: |
| Massless two-loop, D = 13/2 | 40 | 30 | 6.73 s |
| Massless two-loop, Laurent orders -2 through 0 | 30 | 24 | 138.82 s |
| Massive two-loop, D = 12/5, one profile | 10 | 0 | 47.64 s |

The massive profile additionally passes 10 historical same-profile comparisons
and two sector-assembly checks. Its single profile does not establish a new
precision-refinement gate. The independent reference uncertainty is an empirical
refinement estimate, not a rigorous numerical error interval. No supplied oracle
numerical records were compared. Resource times exclude compilation and Nix
setup, and concurrent independent runs prevent a controlled speed comparison.

Both full three-loop requests stop in cut `[0]`, with frontier sizes
`7, 15, 105, 6` over four rounds. They retain three unresolved occurrences of the
same child label under `ConditionVanished { rule: 4218, condition: 0 }`.
The first attempt uses one conditional-point refinement pass and takes 8.43 s;
the separately bound retry permits three passes and takes 8.34 s. Both actually
perform the same one point refinement before exposing the failing child. No
numerical prediction was written, and no three-loop Laurent run was started.
The unchanged native failure remains explicit; enlarging this pass budget did
not resolve it. A [separate addendum](polynomial-closure-v2-domain-retry-addendum.json)
records the next retry with three passes and a shared domain budget of 16,384.
It performs two point refinements, then fails at a different exceptional child
under rule 8961 after four rounds, with frontier sizes `7, 15, 105, 7` and three
unresolved occurrences. It takes 16.16 s and also writes no predictions. The
earlier domain-exhaustion diagnosis was inferred from search/control-flow
evidence; no direct remaining-domain counter was saved. None of these outcomes
is equivalent to the earlier source-only standalone closure.

The final [per-residual control](polynomial-closure-v2-per-residual-control-addendum.json)
keeps the three-pass, 8,192-domain settings and changes only domains per residual
from 32 to 1. It stops at the frontier cap after nine rounds in 3.59 s, with
frontier sizes `3, 15, 45, 105, 210, 378, 630, 990, 1025` and 2,875 unresolved
terms. No three-loop comparison ran because no attempt saved a complete
prediction. The two-loop results above remain valid.

The native-only conditional-chain diagnostics replay saved current requests and
construct individual replacement rules. They distinguish a guarded exceptional
point from invalid physics and document the remaining search/allocation problem;
they are neither a closed full reduction nor numerical amplitude evidence.
Their separate manifests retain exact inputs and immutable proof artifacts.

Completed verbose native payloads and numerical logs are preserved by the
lossless gzip archive map, with original and compressed hashes and restore
verification. Inputs, settings, predictions, comparisons, successful final
closed programs and the specifically reserved failed round checkpoints remain
raw. Existing report hashes refer to original uncompressed bytes; restore an
archived path from the map before rerunning a script that reads that path.
