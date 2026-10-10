# Singleton physical-zero integration: validated corrected attempt

**Passed 218 coherent tests and all lower-loop numerical comparisons.** The gates include 91 native tests and 127 RustFlow/interface/source/capacity tests; one native test remains ignored. The build captures 13 artifacts against 2,485 source assets, including the direct massless-flow test target.

The massless two-loop fixed-D gate passed 40 independent value checks and 30 profile refinements. Its Laurent gate passed 30 independent coefficient checks and 24 refinements. The single massive profile passed 10 independent and 10 historical value checks plus two assembly checks; it is not a new precision-refinement gate. Exact source/build/executable identities and saved predictions precede every reference comparison.

The massless assembly retains all four physical cuts. Two singleton contributions use separately reported physical endpoint-zero certificates, with original coefficient conditions and finite-jet witnesses. They are not finite-eta source zeros, native closures, or transported AMF solutions. The remaining double cut retains actual native source closure, recursive boundary construction, transport and endpoint validation. Massive occupied cuts retain their actual flows. The direct PreparedOccupiedFlow path is unchanged.

The original test failure is preserved in the sibling initial-attempt report. The corrected fixture asserts the actual Cn/H0 representation and its n-1 mass jets, rather than expecting explicit upper-support indices in original input terms. The adjacent comment was corrected; the production algorithm was unchanged between attempts. Report-local metadata validators preserve historical reference definitions, numerical arithmetic and tolerances; old strict-AMF comparator files are unchanged.

`checkpoint-summary.json` records counts/resources; `final-integrity.json` verifies the final source/build state and every successful comparison input. `archive-map.json` provides deterministic gzip and restored hashes for obsolete completed proof payloads and numerical logs. Current closed/provisional proofs and all directly bound comparison inputs remain raw. Timings are shared-host measurements, not speed comparisons.

No full three-loop or four-loop numerical acceptance is claimed. The independent HC selected cut-3 AMF control and the new standalone joint three-loop attempt are separate reports.
