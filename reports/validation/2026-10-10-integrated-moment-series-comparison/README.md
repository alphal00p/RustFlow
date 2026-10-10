# Saved integrated-moment series validation

The source-derived 64-term coefficient prefix agrees with the independent finite-eta reference at eta = 8 for both original targets of the selected three-loop double cut. The scalar and raised mixed-medium 64-term relative errors are 1.14e-25 and 2.68e-24. The final 48-to-64-term changes are 2.38e-20 and 4.23e-19. These are empirical comparisons, not rigorous remainder bounds.

All 12 eta-8 reference checks, all 4 eta-8 truncation checks, and all 12 precision checks pass the prospective criterion. All 12 eta-2 raw-sum reference checks fail and remain recorded: eta = 2 lies outside the source-proved initial disk, so these sums provide no continuation evidence. The independent eta = 0.5 quadrature remains unresolved in the separate frozen reference report; no additional reference nodes were introduced here.

`series-predictions.json` was saved and hashed before `compare-saved.py` read numerical references. The evaluator only reads the producer's frozen rational coefficients and exact formal master identities. `source/normalize.rs` evaluates those masters with the existing native ordinary tadpole owner and applies the existing native measure and exact routing-determinant owners plus the producer-declared compact factor. It preserves every target phase. No coefficient, recurrence, or equation is derived from the reference representation.

`plan.json` fixes 32/48/64 terms, 30/50 digits, and eta = 8/2. Thirteen contract controls reject insufficient prefixes, missing coefficients, duplicate channels, mismatched inputs/cuts, unknown masters, vanished conditions and ambiguous leading branches. The 24-term prefix was explicitly rejected for the planned 64-term evaluation. The first adapter compile failed on two Rust borrow types; its exact source and diagnostics are retained beside the successful build.

The producer's source-based bound gives an analytic normalized germ for |1/eta| < 1/4, which covers eta = 8. `germ-independent-review.json` records the scoped audit. This bound uses the actual off-shell source polynomials, not reference values, and does not supply a numerical tail constant. The first conservative 1/16 bound is retained as historical evidence.

The separate frozen reference report checks the original physical mass derivative with the original numerator and chemical potential fixed, including its nonzero upper surface. This report compares the total raised series; it does not claim a separate producer comparison for each raised component.

There is no endpoint, Laurent, ODE-identity, complete three-loop amplitude, or four-loop acceptance claim. Production source and weighted closure remain untouched. Executable/library hashes are retained without committing the large binaries. No timing comparison is made.

Start with `summary.json`, then `comparison.json`. `producer-binding.json` binds exact source/run copies; `execution-provenance.json` binds the arithmetic runtime; `artifact-manifest.json` binds all local evidence.
