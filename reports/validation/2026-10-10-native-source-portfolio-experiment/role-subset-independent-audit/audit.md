# Original-source subset selection audit

This is a search-design audit, not a closure or numerical result. `select.py`
reads the archived exact 52-row singleton corpus and independently regenerated
geometry roles. It does not choose by source ID, graph name, or known value.
`selection.json` retains the selected original rows, domains and conditions.

The strict predicate requires every protected coordinate to be a symbolic
unchanged index in every original term. It selects **18 rows**, not the 32
virtual differentiated Lorentz rows. Scalar numerator conversion introduces
cut-lowering shifts into 16 virtual Lorentz rows; these include `lorentz/2/0`,
the useful direction in the earlier odd-zero witness. Strict selection also
keeps two multiplication identities on occupation-index-one domains.

Requiring symbolic occupation shifts zero and symbolic required-cut shifts
nonpositive selects **34 rows**: all 32 virtual Lorentz rows plus those same
two multiplication identities. Source IDs describe this result only; they are
not part of the predicate. An optional coefficient-index-dependence filter is
another search heuristic, not a validity theorem, and can exclude useful
constant-coefficient identities. It needs its own matched control.

## Required contracts

* Selection gates traversal of original ordinals in the immutable full guarded
  corpus. Preserve index/variable maps, row coefficients, domains, conditions,
  role array, zero domains and original ordinal identities. A freshly renumbered
  smaller corpus cannot be sealed as if its ordinals referred to the full one.
* `Power::constant(0)` is not `Power::symbolic(0)`. The conservative generic
  selector rejects literal protected indices rather than interpreting their
  numeric value as a relative shift. Empty rows supply no candidate.
* H0 and positive-H source branches retain their original separate guards.
  The occupation-index-one multiplication rows are not valid at H0 or H2.
  Every seed image must still fit the original source guard. A selected row
  must not have its guard enlarged because it preserves occupation shifts.
* Multi-cut selection protects every RequiredCut and Occupation axis. An empty
  useful subset is an expected search outcome, not proof that an integral is
  zero. Full-corpus fallback remains required.
* Required-cut nonpositive original shifts do not by themselves guarantee that
  a final recentered recurrence never raises a cut. Recentring subtracts the
  winning target shifts, and elimination combines differently seeded rows.
  Any proposed final-rule shape guarantee must inspect the sealed candidate.
* Existing `SourceVisitOrder` requires a complete permutation. A bounded
  original-ordinal subset schedule must be a separate search facility; weakening
  ordinary permutation validation would change an unrelated API contract.
* Native role-image validation, guarded zero projection, exact lifting,
  recentering-domain intersection, condition extraction, descent certification
  and independent full-corpus replay remain authoritative. The subset is only
  a search heuristic and supplies no new equations or support theorem.

## Portfolio and budget risks

A shorter RHS can have stronger nonzero conditions, a smaller valid domain,
or fail descent/certification. To preserve baseline coverage, compare valid
sealed rules and retain the full-source candidate wherever the subset candidate
is not valid. Choosing one unsealed candidate solely by RHS length gives no
such guarantee. Fatal arithmetic/sample/exact-lift errors must remain explicit;
budget/depth/exhaustion fallback is a different matter.

Each portfolio branch needs explicit source-attempt, seed/depth, exact-trace
and domain budgets, charged before skipped seeds. Original source count is not
the number of visited attempts. A fair comparison retains the full-source
baseline budget and separately reports the optional subset cost. Selecting a
first indirect candidate can address the motivating failure that direct-only
lookahead never sees, but shorter local rules do not establish global closure
or smaller active frontiers.

No production file was changed and no native search was run for this audit.
