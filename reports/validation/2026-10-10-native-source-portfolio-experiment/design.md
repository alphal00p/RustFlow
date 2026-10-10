# Bounded original-source portfolio — design and point controls

Status: design only. The unchanged native library has replayed four source-only
controls. No portfolio implementation, whole closure or amplitude was run.
The preceding direct-hit-only experiment found no improvement on these 19 points:
the motivating odd-zero target already receives an indirect first winner.

The proposed optional search starts after the ordinary full-corpus first candidate,
whether direct or indirect. Retain that exact candidate. Search two independently
bounded original-source subsets (A, then C) with the same native seed order, integral order,
modular GPLU, exact materializer, guarded source-domain validation, zero projection
and canonicalization. Choose a trial only if its exact canonical RHS is strictly
shorter; ties retain the original. An exact zero is optimal by that local metric.
No physical value, graph name, source ID or historical winning rule selects rows.

## Three algebraic selectors

Every predicate operates on the original symbolic row before seed substitution
or zero projection, requires a nonempty row, and protects every relevant axis.
A literal constant index is rejected on a protected axis; constant zero does not
mean an unchanged symbolic index. Coefficient dependence is checked by exact
polynomial degree in the native `index_variables()` positions, including no test
on names or string representations.

| Selector | Predicate on every term | Rows | Applied points | Raw RHS terms | Baseline plus strictly shorter alternatives |
| --- | --- | ---: | ---: | ---: | ---: |
| Full baseline | All original sources | 52 | 19/19 | 192 | 192 |
| A: strict | Symbolic shift zero on every RequiredCut and Occupation | 18 | 7/19 | 24, incomplete | 160 |
| B: nonincreasing cuts | Symbolic shift zero on Occupation; symbolic shift <=0 on RequiredCut | 34 | 19/19 | 192 | 183 |
| C: index-dependent B | B, and at least one row coefficient depends on a symbolic integral-index variable | 32 | 19/19 | 199 | 181 |

C selects exactly original ordinals 20 through 51 on this corpus. Those are the
historical 32 virtual-derivative rows, but that identification describes this
fixture only. The predicate is a search heuristic, not a universal characterization
of virtual differentiation or a support theorem. A excludes 16 such rows because
inverse-propagator numerator conversion can lower a cut. B also retains two
index-independent multiplication identities, which C excludes. C can exclude
useful constant-coefficient identities on other families.

The final column is a post-hoc diagnostic using the already replayed point rules.
It is not an implemented rule selection result. Every selected alternative has a
different raw nonzero-condition list. In C's two improvements, the added conditions
are nonzero constant multiples of existing ones. A also changes the actual
exceptional locus, for example replacing epsilon conditions by eta. Thus these
term counts never imply unchanged parameter coverage. Full details, source
correspondence, gaps, proof hashes and conditions are in
`selector-controls-comparison.json`.

C gives 16→7 terms at `[2,1,1,2,1,-1,0,0,0,0,0,0]` and 2→0 at the odd target
`[1,0,1,0,1,0,0,-1,0,1,0,0]`. Replacing every baseline rule by C would instead grow
192→199 terms; retaining the baseline matters. An unimplemented four-arm minimum
is 157 terms, with the same condition caveats and greater optional search cost.

## Native implementation boundary

The production source owners remain unchanged. A coherent isolated prototype can
be small:

1. Compute a private list of allowed original source ordinals from the immutable
   guarded corpus. Keep its original relative order. Record a digest and exact
   selector version. The full role array, index-variable map, coefficients,
   source domains, nonzero conditions and measure-specific zero domains stay in
   the original context.
2. Add a separate private subset traversal to
   `crates/rustred-core/src/solver/search.rs::search_attempt_inner`. Its existing
   `source_order` is used by `SourceVisitOrder`, whose public contract requires a
   complete permutation. Do not weaken that contract. At every seed, select stored
   basis ordinals; `SeedSource::basis_row` remains the full original ordinal.
3. Reuse existing `RuleTrialLimits`, `RuleTrialStats`, native bounded trace and
   `TrialSearchError` for an optional **guarded** attempt. Existing ordinary
   `search_attempt` does not pass a `GuardedSearchScope`; the new private guarded
   wrapper must. The selected row must still call `scope.instantiate` with its
   original ordinal and domain. No rebuilt/reindexed source system belongs in
   the prototype.
4. Run the ordinary full-corpus search unchanged first in `solve_guarded_case`.
   On a nonzero candidate, run the optional subset from a fresh probe. No modular
   rows or pivot history from the first attempt are reused. Both candidates have
   the same requested case and native integral ordering.
5. For the point pilot, pre-seal the full candidate and each trial against the
   unchanged full context, original ordering and discovery domain using existing
   `seal_candidate`. Only an alternative that independently seals and covers the
   requested point/domain participates in the RHS-length comparison. A trial
   Certification failure is a recorded rejection; retain the full-source proof.
   Fatal algebra/replay errors remain explicit. If the baseline cannot be sealed
   before the existing sign split, skip optional trials and return it unchanged
   to the existing guarded path. This preserves that path rather than bypassing
   its partition logic. The selected candidate is sealed again by the current
   owner; the duplicate replay cost is explicit for this small isolated pilot.
6. Compare exact canonical RHS length, first A then C, retaining the earlier
   candidate on a tie. Return the winner with separate attempt statistics.
   Existing sign splitting,
   source replay, recentering-domain intersection, descent certification,
   exceptional-condition extraction and application remain authoritative in
   `solver/guarded/search.rs` and `solver/guarded/replay.rs`.

Pre-sealing can preserve the full-source candidate when an optional proof fails
certification or shrinks the requested point/domain. It **does not guarantee
fallback on stronger parameter conditions during later application**. Exact
replayed condition lists must remain attached to the selected rule. Do not claim
baseline parameter coverage solely because both candidates seal. A production
design promising that coverage must retain both through application and explicitly
route exceptional loci to an applicable baseline, or prove that the alternative
introduces no new locus. That is a separate change, not an inference from point
tests. Final code placement remains subject to native owner review.

The current controls use the existing public constructor to build separate subset
contexts. They retain every selected original row, domain and condition byte-for-
byte and record an explicit original-ordinal bijection, then rediscover and replay
rules and round-trip each saved native program. They validate the search idea;
they do not establish full-context traversal/recentering correctness for the
proposed implementation. `control-artifact-relocation.json` resolves historical
execution paths after the completed controls were moved out of the frozen prior
experiment directory.

## Budgets and error policy

Do not consume baseline resources on the optional trial. After retaining a full
candidate, give each trial its own depth, attempted-row and exact-trace limits.
Count a selected row before guard rejection or empty-row projection. Skipped
unselected rows are not attempted, but the selector scan is separately recorded.
An empty subset skips the trial. A no-candidate baseline retains its existing
failure; the initial portfolio does not manufacture coverage by treating an
optional-search result as a baseline fallback.

The existing 512 attempted-row suggestion cannot reach C's odd-zero winner: that
winner occurs at seed 39 of a 32-row subset, requiring between 1217 and 1248 row
attempts. Its 16→7 winner is seed 6, requiring 161–192 attempts. Strict A's odd-zero
winner is seed 103, requiring 1837–1854 attempts. A 2048-row exploratory cap covers
these observed winners, but was not enforced in the current depth-three controls.
They used depth 3 and one requested point domain without an attempted-row cap.
The proposed isolated limits are depth 3, 2048 attempted selected rows, 512
successfully instantiated rows (including empty projected rows), 64 exact-trace
rows and 16384 exact-trace terms **per alternative**. Reuse existing
`RuleTrialLimits` for depth/max_rows/trace limits, with one separate attempted-row
cap. The guarded trial must increment the attempted count before
`scope.instantiate` and increment accepted rows only after it returns Some; the
existing `account_trace=false` placement already supplies the latter count. No
new independent-row or matrix budget is proposed. Existing finite row widths
and the accepted-row bound limit matrix size, without preempting a single native
algebra operation. These controls did not enforce the new limits and do not
silently count as that bounded run.

Depth/seed exhaustion and explicit row/trace-budget exhaustion return the retained
baseline. Arithmetic, unsupported-power, unlucky-sample, malformed source,
canonicalization and exact-lift failures remain explicit according to the prior
experiment's clarified contract. Do not silently label any of them budget fallback.
Stats must distinguish the winner's source proof from total baseline+trial work.

## Required isolated gates before any recommendation

Use a unique experiment persistence identity. Test original full-source ordinals
through encode/decode and exact replay, protected literal-vs-symbolic indices,
multiple cut/occupation axes, empty selectors, untouched source guards and zero
premises, and unchanged ordinary search. Exercise baseline first-direct and
first-indirect cases, a later shorter zero, equal/larger alternatives, every budget
fallback, rejected-row accounting, native exact failure propagation, and changed
conditions/downstream sealing gaps. Then run the same complete 19-point controls
with the actual row/trace limits. Only after demonstrable benefit and accepted
coverage limits should an independently bounded active-closure pilot be considered.
No three-loop or four-loop acceptance changes follow from these controls.

## Proposed two-alternative selection on saved point rules

This table models A followed by C after the full baseline, with strict improvement
only. It is not a new native portfolio run. Both alternatives are independently
replayed point results; their conditional scope remains explicit.

| Point | Full RHS | Selected RHS | Arm | Selected conditions |
| --- | ---: | ---: | --- | --- |
| `[1,1,1,1,2,0,0,0,0,0,0,0]` | 9 | 4 | A | `eta` |
| `[2,1,1,1,2,-1,0,0,0,0,0,0]` | 22 | 4 | A | `eta` |
| `[2,1,1,2,1,-1,0,0,0,0,0,0]` | 16 | 7 | C | `2*eta`, `eta` |
| `[1,0,1,0,1,0,0,-1,0,1,0,0]` | 2 | 0 | A (C ties) | `-3+2*epsilon`, `-6+4*epsilon` |
| `[1,0,1,0,2,0,0,0,0,0,0,0]` | 2 | 1 | A | `-1+epsilon`, `-2+2*epsilon` |

The other 14 full-source rules remain selected, giving 192→157 displayed RHS
terms. `trial-cost-bounds.json` records bounds for the observed winners, including
why 512 attempted rows would miss the odd zero. No performance, active frontier,
closure or amplitude improvement is inferred.

## Reviewed condition gate for the isolated prototype

The selected alternative must first seal against the unchanged full source
context and cover the requested singleton domain. Then every actual alternative
condition must be either a nonzero constant or an exact nonzero rational multiple
of an original sealed baseline condition. Native Symbolica integer-polynomial
cross multiplication verifies this without factoring or division by a parameter.
Zero and foreign variable maps are rejected. Conditions are retained verbatim;
dropping baseline conditions is permitted. The comparison always uses the
original baseline, even if an earlier alternative has already dropped a condition.

This conservative gate rejects the new eta or epsilon loci in three of the five
saved improvements. The two surviving points are 16→7 and the odd target 2→0,
giving 192→181 terms in the saved-point selection model. The independent exact
audit is `role-subset-independent-audit/condition-check.json`. The isolated
prototype remains point-only; rays and ordinary searches retain their original
paths. Full closure, runtime benefit and physical predictions remain untested.
