# Fixed-shell singleton source-order diagnostic

All four controls used the same 19 stored target points, native seed depth 3 and
one search domain per point. Every returned rule passed native source replay.
The original source guards, conditions and explicit measure-zero boxes were
retained. The subset arms remove only the source rows stated below; they are
validation experiments, not production source policies.

| Source presentation | Rows | Applicable point rules | Total RHS terms | New upper-surface terms | RHS identical to original |
|---|---:|---:|---:|---:|---:|
| Original | 52 | 19 | 192 | 2 | 19 |
| Full corpus, virtual rows first | 52 | 19 | 192 | 2 | 19 |
| Virtual rows only | 32 | 19 | 199 | 0 | 15 |
| Virtual rows plus multiplication | 36 | 19 | 192 | 0 | 17 |

The full-corpus reorder changes some exact proof traces but leaves all 19 RHSs
and nonzero-condition lists identical. It therefore gives no evidence for a
simple source-order remedy at these bounds.

Removing compact-derivative sources prevents new upper-surface terms, as
expected. With multiplication retained, two RHSs change: one grows from 24 to
33 terms and another shrinks from 16 to 7. The aggregate RHS count is unchanged.
Dropping multiplication as well changes two additional results: a one-term
surface relation grows to 10 terms, while an odd surface descendant becomes an
exact zero under its retained dimensional conditions. Restoring multiplication
selects its original two-term recurrence again. This is direct evidence that
the first discovered recurrence can obscure a later useful identity; it is not
evidence of a missing physical identity or a minimal basis.

The scalar root and its tested first derivative retain the same 9-term and
20-term recurrences in every control. These tests do not iterate derivative
closure, compute boundary values, or validate an amplitude. Individual runtimes
are preserved but were not collected as an isolated performance benchmark.

The source corpus contains the constant-medium-vector IBPs but no explicit
self-temporal Euler vector field `E_i U`. That polynomial identity is already
in the shifted-index span of the medium-vector sources. The separately owned
matched diagnostic regenerated all 52 original rows, checked that exact span,
and appended or prepended six self-Euler rows: both presentations retained
identical RHSs and conditions on all 19 points. No missing measure theorem or
useful new production source option is inferred from that control either.

These four controls use the preserved canonical native library recorded in
`priority-run-bindings.json`, before explicit measure-zero source projection.
Any comparison with that new implementation must extract the archived source
corpus, rediscover and replay the rules under the new library, and save a
separate versioned result. Old guarded proof programs are not imported into the
new guarded schema.
