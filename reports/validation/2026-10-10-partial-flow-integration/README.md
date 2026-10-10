# Partial occupied-flow integration — safe pause

The coherent build and all 281 tests passed (101 native and 180 RustFlow tests; one native test ignored). All existing lower-loop numerical comparisons also passed:

- Massless fixed dimension: 40 reference checks and 30 refinements; 4.25 s harness time.
- Massless Laurent: 30 reference checks and 24 refinements across five profiles; 96.74 s.
- Massive fixed dimension: 10 independent checks, 10 historical checks and 2 assembly checks; 47.82 s.

At the user's explicit pause request, the active fresh selected three-loop partial-flow producer was stopped with SIGTERM. Its wrapper recorded exit -15 (launcher exit 241), 89.6236 s and 102,576 KiB peak child RSS. It published six provisional rounds with frontiers 14, 93, 222, 590, 835 and 928. It saved no numerical profile and no closed native program. This interruption is not an algorithmic failure, a selected numerical acceptance, or a full three-loop/four-loop result. No comparator was run for it.

`safe-pause-summary.json` and `pause-integrity-audit.json` bind the completed gates, lower comparisons, termination record and all observed partial checkpoints. All 2,496 captured source assets and the runtime executable were unchanged after termination. All owned producer/wrapper processes ended. Checkpoint JSON was parsed and paired program byte lengths checked; no new proof replay or numerical computation was started for the pause audit. The unpublished next-round direct-zero/history artifacts remain identified as incomplete work, without a resume claim.

The production change integrates the private origin, source-image class, finite-label endpoint and virtual-soft boundary owners with the existing AMF pipeline. Original raw coefficient domains are retained before cancellations and support-zero classifications, and the final reduced system is bound after the audited domains are merged. Native replay, original weighted-target and basis-derivative closure, ordinary hard boundary seeds and the existing rational target/Frobenius projector remain required. The isolated predecessor's 50 passing selected tests are included in the coherent total, not counted twice.

The interrupted control used original cuts `[0,3]`, shifted slots `[1,2]`, epsilon `-5/4`, and planned four profiles `18:60:8`, `28:60:8`, `28:80:8`, `28:80:12` sharing one fresh preparation. Its 1,800-second cap was not reached. No saved rules, checkpoints or independent reference values entered the producer. Lower-loop settings and numerical tolerances remain those of the memo predecessor. The selected-sector comparator is prepared but unexecuted; its initial and hardened versions are preserved under `comparator-review/`.

Source-owner limits remain explicit: cancellation is cooperative between operations, and adapter budgets do not bound native CAS scratch or individual CAS latency. Boundary report validation uses existing Debug text plus immutable source binding, not a second theorem prover. Shared-host resource numbers are not speed comparisons.
