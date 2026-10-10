# Conditional scalar ODE transport

This report tests numerical continuation of the scalar finite-data candidate from the sibling physical-fit report. It does not certify that equation, solve the raised target, or establish the original three-loop acceptance gate. No numerical reference was read while constructing the equation, boundary data or predictions.

`derive.py` converts the frozen shift-two / degree-four recurrence exactly. With `t=1/eta`, `beta=-1/4`, and `g=t^beta F`, it retains the forcing

    sum_j t^(2-j) p_j(theta_t-beta-j) g = -(14/117)t^beta.

The companion states are `(g, theta_t g, theta_t^2 g, theta_t^3 g, h=t^beta)`. Since `d/deta=-theta_t/eta`, the first three rows have `-1/eta` on the superdiagonal, the fourth row is `A_k/(eta*A_4)` with forcing coefficient `-B/(eta*A_4)`, and `h'=h/(4*eta)`. All 64 known scalar coefficients verify the conversion, including startup. This is a conversion check, not an independent infinite identity.

The highest eta-derivative coefficient is `-32/1053*eta^2*(eta+1)*(eta+4)`. The ordinary positive-real route avoids its finite singularities. The output `g(eta)=eta^(1/4)F(1/eta)` retains the fractional branch. The common native unit-mass hard master and two compact angular normalizations remain factored out exactly as in the producer's `series-contract.json`; no additional physical normalization is inserted here.

## Completed controls

The unchanged public RustFlow `DifferentialSystem::compile` and `CompiledSystem::transport` owners ran all twelve combinations of requested digits 30/50, boundary terms 32/48/64 and start eta 16/32. Working precision adds 20 guard digits. Taylor transport order is 64. Every native continuation call has a proposal cap at most 128 and visits a segment of eta 16,8,4,2,1. The wrapper subtracts prior accepted steps when setting each later call cap. It does not count earlier rejected proposals across calls, so an aggregate 128-proposal bound is not established; successful-call rejected counts were not retained. All 60 finite checkpoints succeeded. The largest relative difference from the finest predeclared profile is `1.016e-21`.

The existing public Frobenius owner then prepared the exact candidate matrix, evaluated orders 32/48/64, matched at eta=1/8, and selected the eta-zero limit. The exact exponents are `0,1/4,5/2,7/2,19/4`. All twelve profiles returned a finite conditional limit, with largest relative refinement `4.470e-22`. The finest normalized candidate value is

    0.0128065020700745910831852094348846972289735933415900725...

An additional finest-profile matching-point control at eta=1/16 exhausted its unchanged remaining call cap (88 proposals after subtracting 40 prior accepted steps). Its exact failure is preserved; it was not retried or used to tune the equation. Thus matching-point independence is not established by that control. The successful profiles vary precision, initial eta, boundary order and Frobenius order at the common matching point.

The finite process used 2.383 seconds and 15,436 KiB peak RSS; the endpoint process used 1.994 seconds and 15,444 KiB. Their combined 4.376 seconds remains below the fixed 300-second transport cap. The extra failed endpoint control is included in those measurements. Compile time is recorded separately and is not numerical runtime.

`physical_limit` uses the existing numerical endpoint-tolerance policy; this isolated public API does not invoke the stronger private physical admission owner. These outputs remain conditional numerical diagnostics of an unproved equation at fixed D=13/2. They are not a symbolic-epsilon, meromorphic-continuation or physical endpoint certificate. Any independent reference comparison must be recorded separately after the saved predictions.

## Reproducibility

`build.json` and `endpoint-build.json` bind the standalone sources, executables and copied exact RustFlow/serde libraries. Builds use `rustc`, do not run Cargo, and do not modify the shared target. The original syntax-only compile failures and the attempted missing-executable launch are retained under `historical/`. `finite-run.json` and `endpoint-run.json` contain process resources; their adjacent provenance files bind actual arguments. `companion.json` binds the selected recurrence and scalar coefficient input files. The raised-target holdout was never read by these scripts.
