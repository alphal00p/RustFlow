# Frozen per-target reconstruction protocol

No physical fit has begun. This prospective protocol selects a separate scalar candidate for each integrated target, jointly fitting that target's exact normalization/exponent channels. Distinct scalar/raised targets are not forced into one common annihilator.

Every target attempts the bounded direct theta family first. If and only if training reports no candidate within the fixed shapes, it tries the predeclared bounded recurrence family. Ambiguous/zero prefixes remain explicit. All attempts across targets share one300s wall clock. Each family keeps order/degree4, <=25unknowns and <=24shapes, for at most48shapes per target. All selected candidates are frozen before either target's holdout is opened. No holdout can trigger another fit or family change.

`run-target-fits.py EXE CONFIG NEW_OUTPUT_DIR` consumes `integrated-per-target-fit-plan.v1` with `targets:[{id,training,holdout}]`. Each split file uses the frozen fitter schema plus `target_index:id`. Targets run in numeric ID order. The final executable remains the frozen reconstruction fitter; this launcher does not edit/rebuild it.

Five synthetic checks pass in0.780s: distinct target operators remain first order; every candidate is frozen before any holdout; and c_n=1/(n+2)^4 correctly exercises theta-order4 miss followed by a degree4 recurrence candidate. Its coefficient generating function has an order4 equation with startup forcing, demonstrating why the predeclared recurrence fallback is useful. No synthetic holdout test is an identity certificate. Earlier theta-only synthetic results are retained as history.

A recurrence remains a coefficient-core relation until startup forcing and beta conjugation are handled. Separate certified companions may later form a block system; no homogeneous full-integral operator or physical period is asserted here. Physical datasets/results will live in separate directories after the producer supplies exact64-coefficient splits and source-bound grading/normalization metadata.
