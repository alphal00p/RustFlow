# Integrated finite-density moment-flow exploration

The polynomial boundary construction succeeds for the actual nonzero three-loop double cut. Both original targets have 64 exact integrated coefficients, generated without occupation-aware IBP or reference inputs. A reliable integrated flow for both targets has not been established: the scalar has a held-out-validated recurrence candidate, while the raised target has no candidate within the declared limits. Conditional transport and direct-series continuation are recorded separately below.

This is an isolated exploration. Production code and the previous closure work are preserved. No full three-loop, Laurent, or four-loop acceptance follows from these results.

## Input and fixed protocol

The input is `examples/finite_density/massless_three_loop_chain.json`, SHA256 `edee8ee4870bd4aafc7b56100a32557085b3ac21e0e47fb44d74560b054dd603`. The actual nonzero cut is `[0,3]`. Both the scalar and raised mixed-numerator target are retained. The shells and occupied regions stay fixed; all uncut physical slots `[1,2,4]` receive the auxiliary mass shift. The first fixed dimension is `D=13/2`.

`plan.json` was committed before physical fitting or finite-eta reference generation. Coefficient generation is limited to 64 terms per branch and 900 seconds, with a 24-term/180-second pilot. Each target reserves 48 training and 16 held-out coefficients. Recurrence shift, coefficient degree, theta order and variable degree are each at most four, with at most 25 unknowns per candidate and 300 seconds for the entire fitting stage. The prospectively frozen per-target protocol tries a homogeneous theta equation first and a recurrence only if training finds no theta candidate. All candidate hashes are frozen before any held-out file is opened. Failure does not trigger larger orders or reuse of held-out data for selection.

## Boundary construction and reusable algebra

The native region census has exactly two regions. The all-soft virtual polynomial is scaleless, with native ordinary zero evidence. The virtual-hard region leaves ordinary massive vacuum integrals multiplying compact-momentum polynomials. There are no residual compact denominators after the complete all-uncut expansion for this input. The scalar uses even half-grades; the raised terms use odd half-grades with a compensating leading offset. Their integrated branch is `eta^(D/2-3)` times an integer inverse-eta series. At `D=13/2`, this is `eta^(1/4)` times the series. The demonstrated nonresonant branch has no logarithms. A Laurent expansion in the regulator would generate logarithms from the dimension-dependent power and requires a separate calculation.

Two generic extensions make high orders affordable. Gaussian integration before expansion converts polynomial Schwinger monomials into dimension-shifted ordinary vacuum integrals. A Gram-determinant identity converts them back to the original dimension. Aggregated Wick contractions then integrate independent compact angular monomials, with radial moments and surface distributions supplied by the existing dimensional moment rules. Neither extension selects a physical invariant or inserts an evaluated formula for this graph.

The angular adapter passes a rank-128 check and an independent audit of 1,392 Cartesian moment comparisons for three and four vectors. The native Gram dimension bridge passes 16 one-to-four-loop Gaussian checks. RustRed's compact search indices stop at 63, but its existing proof-backed ordinary vacuum `Reducer` accepts wide indices; checks through 129 pass without changing RustRed's representation or authoring a recurrence. The physical producer preserves off-shell diagonal terms until the cut product identities are applied, including the raised-cut and Fermi-surface contributions.

The Gaussian producer exactly matches the native region/TensorProjector route through four coefficients of each target, and the wide-index consumer matches the original ordinary backend through 24. Its focused moment suite passes 19 tests, including 144 radial normalization cases, 144 upper-surface identities and 60 polynomial comparisons. The final 64-term generation takes 0.766 seconds; all coefficient runs, including the preserved failed compact-index attempt, total 2.0393 seconds.

## Independent validation and equation search

The complete off-shell source polynomials give a convergent normalized germ for `|1/eta| < 1/4`. Gaussian/simplex and compact-domain bounds establish this without reference values. The fractional prefactor uses the branch continued from positive real eta. This local statement is not an endpoint error bound or a four-loop continuation theorem.

Saved 64-term predictions at eta=8 agree with independent references to relative errors of `1.14e-25` for the scalar and `2.68e-24` for the raised target. The final 48-to-64 changes are `2.38e-20` and `4.23e-19`; all precision checks pass. Raw eta=2 sums fail, as expected outside the initial disk. Independent eta=8 and eta=2 quadrature refinements pass their prescribed accuracy checks; the eta=1/2 initial reference grid remains unresolved. A separate fixed-original-mass diagnostic checks bulk, surface and total raised contributions; it is a derivative diagnostic, not a replacement for the precision reference.

The scalar recurrence has shift order two and polynomial degree four. Its exact training matrix has rank 14 for 15 unknowns, using 46 training equations, and all 16 held-out equations pass. Neither homogeneous-theta nor recurrence fitting finds a raised-target candidate within the fixed shapes. That target's held-out file remains unopened by the fitter. The total physical fit takes 0.845 seconds. Finite-prefix agreement establishes a candidate only; it is not an all-index identity.

## Conditional flow and direct continuation

The scalar recurrence converts to a fourth-order forced equation, retaining the startup term `B(t)=-14/117` and the fractional exponent `beta=-1/4`. Finite candidate transport succeeds for 12 independent precision/truncation/start profiles, with maximum relative refinement `1.016e-21`. This verifies numerical transport of the candidate equation; it does not prove that equation for the integral.

The existing Frobenius owner gives a finite conditional scalar endpoint for all 12 profiles at matching point eta=1/8, with maximum relative refinement `4.470e-22`. The exact candidate exponents are `0,1/4,5/2,7/2,19/4`. An extra finest-profile matching-point control at eta=1/16 hits its unchanged step cap and is retained without retry. Thus precision and truncation stability are demonstrated at a common matching point, but that additional matching-point independence check fails. The total finite and endpoint transport runtime, including this failure, is 4.376 seconds. The public endpoint API uses its existing numerical tolerance policy; these are conditional diagnostics rather than a physical endpoint certificate.

After the conditional outputs were frozen, the independent comparison applied the native master/measure prefactor and saved physical predictions before reading references. All 24 comparisons pass: 12 at eta=2 and 12 at eta=0. Their worst relative errors are `8.815e-22` and `4.470e-22`, respectively. This adds independent numerical evidence for the scalar candidate without establishing an all-order identity, a raised-target result or a Laurent result.

The proposed generic symbolic-index certificate interface is not yet supplied by the integer-order coefficient producer, so no all-index recurrence proof is claimed. There is also a concrete obstruction to the fixed direct-eta certificate ansatz: the candidate acting on the source parameter kernel has a nonzero fourth-order pole, while the permitted certificate factors have poles of order at most two and their divergence can have order at most three. Exact substitution at an interior source point gives a nonzero leading coefficient `-7595/11501568`. This excludes that bounded ansatz only; it is not a proof that the candidate is false or that no certificate exists.

The fallback uses the same 64 coefficients, with no new integration or fitted reference data. For `w=4/(eta+4)`, it forms `H(w)=(1-w)^(1/4) G(w/[4(1-w)])` by exact binomial composition. Its 64-term sums at eta=2 agree with independent references to relative errors `2.41e-16` and `5.20e-17` for the scalar and raised targets. A source-only Cauchy argument certifies scalar truncation bounds of `3.135e-28` relative at eta=8 and `2.839e-8` at eta=2; numerical rounding is a separate issue. An independently audited raised bound sums absolute residue, original-numerator, kernel-derivative and upper-surface contributions rather than using the signed leading coefficient as a norm. It gives relative truncation bounds `6.609e-27` at eta=8 and `3.959e-7` at eta=2 on the same grid. These conservative bounds establish finite-eta control, but do not certify the stronger precision suggested by numerical agreement at eta=2.

At eta=0, the 64-term fallback errors are `5.22e-3` and `2.19e-4`. The largest tested Padé approximants have errors near `5.24e-5` for both targets, with four endpoint Padé precision checks failing. The transformed endpoint is `w=1`, where the finite-radius Cauchy bound gives no geometric tail control. Therefore neither endpoint accuracy nor a reliable direct-continuation fallback is established within this prefix budget.

## Present limits

The coefficient prototype accepts original virtual numerator degree at most one. Its first wide-artifact adapter has deliberately narrow exact family and parameter admission; other vacuum families still use the ordinary backend. Higher Gaussian tensor numerators, multi-hard-loop ultraviolet face control, and production integration remain separate work. The retained polynomial/report and native vacuum-cache budgets are explicit, but they do not bound all temporary CAS memory.

The original requirement to validate the complete nonzero three-loop case, including its Laurent result, before four-loop numerical work remains in force. No such acceptance is claimed here.

The bounded exploration is complete with a positive boundary result and an unresolved integrated-transport result. The reusable moment and ordinary-vacuum machinery is worth retaining. Replacing production with the fitted scalar equation would be premature: the raised target misses the declared reconstruction limits, the scalar identity lacks a certificate, and direct continuation fails the endpoint criterion. Occupation-aware closure remains paused, and no four-loop numerical run is justified by this exploration. A subsequent attempt needs a separately specified bounded method for certifying the integrated equation (including boundary terms) or controlling the same series at its endpoint; simply enlarging the present fit is not evidence of reliability.

## Evidence

- `../2026-10-10-integrated-moment-expansion/status.md`: frozen coefficients, region proof, dimensional moments and resource records.
- `angular-aggregate/` and `gram-shift/`: generic algebra identities, source, tests and preserved failed attempts.
- `wide-vacuum/`: existing ordinary vacuum consumer beyond the compact index range.
- `../2026-10-10-integrated-flow-reconstruction/` and `../2026-10-10-integrated-flow-target-protocol/`: fitter tests and frozen selection limits.
- `../2026-10-10-integrated-flow-physical-fit/`: physical candidate selection and held-out evidence.
- `../2026-10-10-integrated-flow-conditional-transport/`: exact recurrence conversion, conditional transport, Frobenius endpoint and preserved failed control.
- `../2026-10-10-integrated-flow-conditional-validation/`: independent comparisons of saved conditional scalar predictions at eta=2 and eta=0.
- `../2026-10-10-integrated-moment-validation/`: independent references and physical mass/surface diagnostics.
- `../2026-10-10-integrated-moment-series-comparison/`: predictions saved before numerical reference comparisons.
- `../2026-10-10-integrated-moment-continuation/`: exact same-prefix transformation, source bounds, finite-eta checks and failed endpoint fallback.
- `../2026-10-10-integrated-continuation-independent-audit/`: source, cut signs, normalization, conditional equation and separate raised-bound review; 22 exact checks pass.
- `checkpoint-verification.json`: final frozen-manifest integrity check and unchanged production-tree check.

All checkpoints use author and committer Kaapo Seppänen <seppanenkaapo@gmail.com>. No changes have been pushed.
