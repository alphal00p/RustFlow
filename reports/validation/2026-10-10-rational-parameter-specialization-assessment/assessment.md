# Exact rational parameter specialization: design assessment

This is a read-only design assessment. No native API, production source, numerical choice, or build was changed. Exact epsilon specialization before native coefficient accumulation is a plausible way to reduce bivariate rational-function work, but no performance benefit has been measured. It does not replace the current validation batch.

## Smallest useful native experiment

`GuardedProgram::apply` currently specializes integral-index variables at an integer point, checks every retained rule condition, checks each original RHS denominator before division, and only then removes zero coefficients. `reduce` repeatedly applies this owner while accumulating exact coefficients. The existing indexed-field specialization API specializes integer indices; it is not a public rational parameter-assignment API.

An additive immutable specialization view could bind an already replayed symbolic program to an exact ordered parameter assignment, initially epsilon alone while eta remains symbolic. Keep source definitions, proof objects, integral domains, order, zero domains and terminals immutable. Replay arbitrary decoded proofs symbolically exactly as today. Specialize both stored polynomial guards and coefficient numerator/denominator separately before native accumulation, and preserve their original provenance. Rational substitution can clear only known nonzero constant denominators to recover integer-coefficient polynomials. Variable-map conversion must be exact and checked; retaining an unused epsilon coordinate is also possible and may avoid a new map owner.

Reject unknown, duplicate, index-variable or eta assignments; zero rational denominators; foreign coefficient maps; and arithmetic/budget failures. Test a specialized guard for exact zero before cancellation, including zero-numerator and zero-target cases. A vanished guard or coefficient pole may use the existing later-rule fallback, but must remain an explicit failure if no admissible rule applies. A nonzero eta-dependent guard must remain a condition for the transport path; it is not automatically satisfied. Preserve all previously accumulated conditions and original input/target coefficient conditions. A specialized zero RHS is valid only after all these checks. This produces a consequence at the assigned epsilon, never a generic-epsilon proof or a new geometric zero certificate.

The first bounded diagnostic should use a copied authenticated program and the already identified expensive required label. Compare its exact result and all guard domains against exact substitution into the saved generic result, retain eta symbolically, and report decode/replay separately from application. Include exceptional samples, a guard that vanishes before cancellation, a canceled target with a pole, an eta-dependent remaining pole, and later-rule fallback. No reference values or new epsilon choices based on reference agreement are needed. This assessment does not launch that experiment.

## RustFlow integration obligations

The current occupied preparation constructs sources and closes a symbolic-epsilon reduced system before `evaluate_report` receives a sample. Specializing only inside the current numerical evaluator cannot remove the expensive symbolic closure applications. A useful integration therefore needs a distinct per-sample reduction/preparation owner or a reusable symbolic proof graph with exact sample-specific application. Symbolic source discovery/replay may be shared; the final original-target and every-basis-derivative audits must still be completed for the actual sample.

A sample can make coefficients vanish, invalidate pivots or change the reached basis. Either retain a separately proved common basis, or bind each exact sample to its own ordered basis, matrix, target maps, guards and mandatory closure audit. Do not treat one sample's closure as a generic closure or import another sample's candidate maps. The current Laurent driver prepares once and samples that preparation; a per-sample preparation design changes this contract even though interpolation of fully validated complete-amplitude values remains possible.

Boundary regions must use the actual original family and sample basis, including raised-shell and occupation labels, physical masses, moving upper surfaces and normalization. The sealed symbolic-D source/origin/contour statements are meromorphic identities, so they can be evaluated at a regular exact sample outside the initial high-D proof domain. This requires the existing finite-label admission and guard checks, actual dimension-dependent seed checks, and an admissible continuation; specializing a coefficient does not create a new physical admission. Sample-specific zeros cannot justify discarding unaudited boundary regions.

Cache and checkpoint identities must separate the generic symbolic program from its exact assignment and resulting system. At minimum bind the reduced rational numerator/denominator of epsilon, parameter names/maps, symbolic source/program identity, ordered basis, original and specialized conditions, physical admission/deformation/normalization, source policy, and any numerical precision/order data. Existing Frobenius caches compare the full matrix, but a new specialized-system API must also carry branch provenance and sample identity; equal numeric matrices alone do not prove identical physical branch selection.

## Endpoint blocker to a naive implementation

`project_limit` deliberately skips a column when its exponent has nonzero epsilon derivative. Substituting epsilon into the differential system before constructing its Frobenius basis erases this derivative. A dimension-dependent branch can then look like an ordinary rational power, including an integer or zero power. The existing `ensure_generic_indicial` also compares generic and specialized exponents to reject additional resonances; the generic information cannot simply be discarded.

Thus native exact specialization is algebraically separable from a complete AMF implementation. Numerical-flow integration needs either authenticated generic indicial/region branch information retained through a compatible specialized recurrence, or a separately justified exactly specialized endpoint construction that implements the original meromorphic limit. It must handle or explicitly reject additional resonances and preserve rational target weights before endpoint projection. This is a production blocker: generic exponent dependence must not be inferred from one or two numerical samples. A naive specialized matrix passed to today's projector is unsafe. The previously observed D=7 resonance is a concrete limitation to retain, not a sample to tune around silently.

## Source locations reviewed

- `vendor/rustred/crates/rustred-core/src/solver/guarded/lifecycle.rs`: `apply` and `reduce`.
- `vendor/rustred/crates/rustred-core/src/solver/instantiate.rs`: integral-index polynomial specialization.
- `vendor/rustred/crates/rustred-core/src/algebra/indexed/specialization.rs`: checked integer-index specialization and limits.
- `vendor/rustred/crates/rustred-core/src/algebra/coefficient/{model,context}.rs`: coefficient ring and variable-map ownership.
- `src/finite_density/flow.rs`: symbolic preparation and sample evaluation.
- `src/finite_density/flow_boundary.rs`: actual-label boundary admission and symbolic exponent matching.
- `src/finite_density/assembly.rs`: complete-amplitude sampling.
- `src/engine/connection.rs`: exact sample guard checks before empty-basis return, boundary/transport/endpoint stages.
- `src/engine/frobenius_cache.rs`, `src/engine.rs`, `src/frobenius.rs`: cached generic preparations, resonance check and epsilon-dependent physical projection.
