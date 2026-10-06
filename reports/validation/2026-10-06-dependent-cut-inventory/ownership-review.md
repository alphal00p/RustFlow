# Scoped independent review

The parent agent approved the native-inventory/grouped-projection design before
implementation. A separate algebraic-solver agent reviewed the implementation
read-only, without edits or Cargo execution in this worktree.

The review covered original cut/virtual prescriptions, positive regulator rates,
the normalized raised-cut identity, exact coefficient guards and sums before
fitting. It found one substantive omission: original numerator denominator
guards were collected too late to survive removed-cut terms and cancellations.
The correction reuses `physical_conditions::rational_denominator_conditions`
before native rewriting, partial fractions or filtering, and validates the
captured guards at every exact epsilon sample, including empty groups. The
reviewer confirmed that this resolved the finding; regression cases cover a
forbidden epsilon and regular zero samples.

A separate native-owner audit identified a statistical-zero defect in the
pinned HEPKit affine parser for generic symbolic Add coefficients. This cut
entrypoint admits exact real rational specialization only. It now reconstructs
the native denominator rows and external Gram values from the already admitted
canonical coefficients before calling the owner. The independent reviewer
confirmed that this preserves physical slot/momentum order and ensures actual
Num coefficients reach native partial fractions. The exact `10^-1000`
raised-cut regression makes a dropped term observable. This is a scoped exact
input adapter, not a replacement partial-fraction implementation, nor a claim
that generic symbolic or Gaussian quadratic routing is supported here.

The first source-matched runtime gate exposed a separate ordering bug in domain
admission for a numerator proportional to a native cut denominator. Replacing
the tensor dimension inside a native scalar-product call before identifying
that call changed its structural identity. The correction identifies native
scalar products first, retains exact grouping on the original native atoms,
and decodes scalar coefficients afterward. The independent reviewer confirmed
that represented denominator guards and scalar-symbol admission remain intact.
The forbidden-epsilon/cut-zero regression requires this corrected path.

The review found no further blocker in that admitted scope. It does not replace
the runtime and independent numerical gates recorded in the final validation
manifest. Wider native-owner corrections and dependency updates require their
own fresh integration validation.
